"""E4: graph neural network with interventional counterfactual extraction.

All guides are modelled jointly. A message-passing network is trained to predict the
log1p-normalised outcome vector of every cell from (i) its own features (cell type, sample,
local density, own guide, whether it has a guide call) and (ii) messages from its neighbours
on the per-sample Delaunay graph pruned at D_max = ``bins_um[-1]``. Each message is an MLP of
the neighbour's hidden state, the neighbour's own-guide one-hot and the edge distance, so the
network can represent distance-dependent spillover from any guide while adjusting for every
other guide in the neighbourhood (the "nuisance exposure" that E1 to E3 remove by restricting
to clean controls).

Estimands are then read off by interventional counterfactuals on the fitted network:

- autonomous effect of target k: for eligible control cells (``group_masks`` column ``C[:, k]``,
  so the control policy and ``clean_controls`` are respected exactly as in E2), set the cell's
  own-guide one-hot to k and ``has_call`` to 1, and take the change in its prediction;
- spillover of target k in ring b: for the same control cells, pick one random graph neighbour
  whose edge distance lies in ring b, set that neighbour's own guide to k, and take the change
  in the recipient's prediction. Recipients with no graph neighbour in ring b are dropped for
  that ring (the Delaunay graph contains only a subset of the pairs within D_max, so far rings
  are thinner than in E1 to E3).

Counterfactuals are evaluated in batches: a batch may intervene on several nodes at once, and a
recipient is kept only when exactly one intervened node lies inside its receptive field
(``n_layers`` hops, checked by label propagation on the graph), so every reported difference is
an exact single-neighbour counterfactual. Estimates are the mean difference over recipients, per
cell type and overall. Standard errors are a 200-resample bootstrap over recipients of the
per-cell differences and p-values are two-sided normal.

Caveat on uncertainty. The bootstrap captures only the sampling variability of the recipient
set given the fitted network. It does not include model uncertainty (the variability of the
fitted function across training runs, seeds or architectures), so the reported intervals are
narrower than the true uncertainty and the mandatory NTC calibration is the check that matters
for this estimator. Training is also only loosely tied to the estimand: the network is fitted
by prediction loss, not by the counterfactual contrast, so weakly supported contrasts (few
perturbed cells of k, few edges in ring b) are extrapolations of the fitted function.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from spatialspill.estimators.base import EstimateTable, Estimator
from spatialspill.estimators.groups import GroupConfig, group_masks
from spatialspill.exposure import Exposure
from spatialspill.graphs import GraphConfig, build_graph

if TYPE_CHECKING:
    import torch


def _make_model(n_in: int, n_guide: int, hidden: int, n_out: int, n_layers: int) -> torch.nn.Module:
    import torch
    from torch import nn

    class SpillConv(nn.Module):
        """One round of message passing; messages depend on [h_j, own-guide_j, d_ij]."""

        def __init__(self, hidden: int, n_guide: int) -> None:
            super().__init__()
            self.msg = nn.Sequential(
                nn.Linear(hidden + n_guide + 1, hidden), nn.ReLU(), nn.Linear(hidden, hidden)
            )
            self.upd = nn.Sequential(nn.Linear(2 * hidden, hidden), nn.ReLU())

        def forward(
            self,
            h: torch.Tensor,
            g: torch.Tensor,
            src: torch.Tensor,
            dst: torch.Tensor,
            e: torch.Tensor,
        ) -> torch.Tensor:
            m = self.msg(torch.cat([h[src], g[src], e], dim=1))
            agg = torch.zeros_like(h).index_add_(0, dst, m)
            return self.upd(torch.cat([h, agg], dim=1))

    class Net(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.enc = nn.Sequential(nn.Linear(n_in, hidden), nn.ReLU())
            self.convs = nn.ModuleList([SpillConv(hidden, n_guide) for _ in range(n_layers)])
            self.head = nn.Linear(hidden, n_out)

        def forward(
            self,
            x: torch.Tensor,
            src: torch.Tensor,
            dst: torch.Tensor,
            e: torch.Tensor,
            guide_cols: torch.Tensor,
        ) -> torch.Tensor:
            g = x[:, guide_cols]
            h = self.enc(x)
            for conv in self.convs:
                h = h + conv(h, g, src, dst, e)
            return self.head(h)

    return Net()


class E4GNN(Estimator):
    name = "E4_gnn"

    def __init__(
        self,
        hidden: int = 64,
        n_layers: int = 2,
        epochs: int = 300,
        lr: float = 1e-3,
        min_cells: int = 5,
        groups: GroupConfig | None = None,
        n_boot: int = 200,
        max_recipients: int = 5000,
        loss_subsample: int = 200_000,
        max_passes: int = 500,
        device: str = "auto",
        mps_min_nodes: int = 20_000,
        seed: int = 0,
        report_by_cell_type: bool = True,
    ) -> None:
        self.hidden = hidden
        self.n_layers = n_layers
        self.epochs = epochs
        self.lr = lr
        self.min_cells = min_cells
        self.groups = groups or GroupConfig()
        self.n_boot = n_boot
        self.max_recipients = max_recipients
        self.loss_subsample = loss_subsample
        self.max_passes = max_passes
        self.device = device
        self.mps_min_nodes = mps_min_nodes
        self.seed = seed
        self.report_by_cell_type = report_by_cell_type

    # ----------------------------------------------------------------- inputs
    def _device(self, n: int) -> Any:
        import torch

        if self.device != "auto":
            return torch.device(self.device)
        if torch.backends.mps.is_available() and n >= self.mps_min_nodes:
            return torch.device("mps")
        return torch.device("cpu")

    @staticmethod
    def _edges(adata: AnnData, d_max: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Directed edge list (both directions) of the per-sample Delaunay graph pruned at D_max."""
        tmp = AnnData(obs=adata.obs[["sample"]].copy())
        tmp.obsm["spatial"] = np.asarray(adata.obsm["spatial"], dtype=float)
        build_graph(tmp, GraphConfig(kind="delaunay", max_edge_um=d_max))
        D = tmp.obsp["spatial_distances"].tocoo()
        return D.row.astype(np.int64), D.col.astype(np.int64), D.data.astype(np.float64)

    @staticmethod
    def _features(adata: AnnData, exposure: Exposure) -> tuple[np.ndarray, np.ndarray, int]:
        """Node feature matrix, indices of the own-guide block, index of the has_call column."""
        obs = adata.obs
        ct = pd.get_dummies(obs["cell_type"].astype(str)).to_numpy(dtype=np.float32)
        smp = pd.get_dummies(obs["sample"].astype(str)).to_numpy(dtype=np.float32)
        dens = np.log1p(exposure.n_neighbors_within.astype(np.float64))
        dens = (dens - dens.mean()) / (dens.std() + 1e-8)
        G = exposure.label_matrix(exposure.labels).toarray().astype(np.float32)
        has_call = (exposure.labels != "none").astype(np.float32)
        X = np.column_stack([ct, smp, dens.astype(np.float32), G, has_call[:, None]])
        off = ct.shape[1] + smp.shape[1] + 1
        guide_cols = np.arange(off, off + G.shape[1])
        return X.astype(np.float32), guide_cols, X.shape[1] - 1

    # --------------------------------------------------------------- training
    def _train(
        self,
        model: Any,
        x: Any,
        src: Any,
        dst: Any,
        e: Any,
        y: Any,
        guide_cols: Any,
        rng: np.random.Generator,
    ) -> None:
        import torch

        n = x.shape[0]
        idx = None
        if n > self.loss_subsample:
            idx = torch.as_tensor(
                rng.choice(n, self.loss_subsample, replace=False), device=x.device
            )
        opt = torch.optim.Adam(model.parameters(), lr=self.lr)
        model.train()
        for _ in range(self.epochs):
            opt.zero_grad()
            pred = model(x, src, dst, e, guide_cols)
            if idx is not None:
                loss = torch.mean((pred[idx] - y[idx]) ** 2)
            else:
                loss = torch.mean((pred - y) ** 2)
            loss.backward()
            opt.step()
        model.eval()

    # ---------------------------------------------------------- counterfactuals
    @staticmethod
    def _propagate(lab: Any, src: Any, dst: Any, reduce: str, hops: int) -> Any:
        for _ in range(hops):
            out = lab.clone()
            out.scatter_reduce_(0, dst, lab[src], reduce=reduce, include_self=True)
            lab = out
        return lab

    def _counterfactual(
        self,
        model: Any,
        x: Any,
        src: Any,
        dst: Any,
        e: Any,
        pred0: Any,
        recipients: np.ndarray,
        intervened: np.ndarray,
        k_col: int,
        has_call_col: int,
        guide_cols: Any,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Per-recipient prediction change when ``intervened[i]`` is set to target column ``k``.

        Returns ``(diffs, used)``: differences (n_used x n_out, standardised scale) and the indices
        of the recipients that were evaluated.
        """
        import torch

        n = x.shape[0]
        dev = x.device
        remaining = rng.permutation(len(recipients))
        per_pass = max(32, n // 40)
        diffs: list[np.ndarray] = []
        used: list[np.ndarray] = []
        passes = 0
        while len(remaining) and passes < self.max_passes:
            passes += 1
            cand = remaining[:per_pass]
            rest = remaining[per_pass:]
            rec = recipients[cand]
            J = intervened[cand]
            Ju = np.unique(J)
            lab_id = torch.arange(1, len(Ju) + 1, dtype=torch.float32, device=dev)
            Ju_t = torch.as_tensor(Ju, device=dev)
            lab_max = torch.zeros(n, device=dev)
            lab_max[Ju_t] = lab_id
            lab_min = torch.full((n,), float("inf"), device=dev)
            lab_min[Ju_t] = lab_id
            lab_max = self._propagate(lab_max, src, dst, "amax", self.n_layers)
            lab_min = self._propagate(lab_min, src, dst, "amin", self.n_layers)
            rec_t = torch.as_tensor(rec, device=dev)
            ok = (lab_max[rec_t] == lab_min[rec_t]).cpu().numpy()
            if len(cand) == 1:
                ok[:] = True
            if ok.any():
                x_cf = x.clone()
                Jok = torch.as_tensor(np.unique(J[ok]), device=dev)
                x_cf[Jok[:, None], guide_cols[None, :]] = 0.0
                x_cf[Jok, guide_cols[k_col]] = 1.0
                x_cf[Jok, has_call_col] = 1.0
                with torch.no_grad():
                    pred1 = model(x_cf, src, dst, e, guide_cols)
                d = (pred1[rec_t[ok]] - pred0[rec_t[ok]]).cpu().numpy()
                diffs.append(d)
                used.append(rec[ok])
            remaining = rng.permutation(np.concatenate([rest, cand[~ok]]))
        if not diffs:
            return np.empty((0, pred0.shape[1])), np.empty(0, dtype=int)
        return np.concatenate(diffs), np.concatenate(used)

    # ------------------------------------------------------------------ output
    def _rows(
        self,
        target: str,
        kind: str,
        ring: int,
        diffs: np.ndarray,
        used: np.ndarray,
        treated_mask: np.ndarray,
        ct: np.ndarray,
        groups_out: list[str],
        outcome_names: list[str],
        rng: np.random.Generator,
    ) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        ct_used = ct[used]
        for grp in groups_out:
            sel = np.ones(len(used), dtype=bool) if grp == "all" else ct_used == grp
            n_rec = int(sel.sum())
            if n_rec < 2:
                continue
            nt = int(treated_mask.sum() if grp == "all" else (treated_mask & (ct == grp)).sum())
            d = diffs[sel]
            est = d.mean(axis=0)
            W = rng.multinomial(n_rec, np.full(n_rec, 1.0 / n_rec), size=self.n_boot)
            boot = (W @ d) / n_rec
            se = boot.std(axis=0, ddof=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                z = est / se
            pv = 2 * stats.norm.sf(np.abs(z))
            ident = n_rec >= self.min_cells and nt >= self.min_cells
            for g, oname in enumerate(outcome_names):
                ok = bool(ident and np.isfinite(pv[g]))
                rows.append(
                    {
                        "estimator": self.name,
                        "target": target,
                        "kind": kind,
                        "ring": ring,
                        "cell_type": grp,
                        "outcome": oname,
                        "estimate": float(est[g]),
                        "se": float(se[g]),
                        "ci_low": float(est[g] - 1.959964 * se[g]),
                        "ci_high": float(est[g] + 1.959964 * se[g]),
                        "pvalue": float(pv[g]) if ok else np.nan,
                        "n_treated": nt,
                        "n_control": n_rec,
                        "identified": ok,
                    }
                )
        return rows

    # --------------------------------------------------------------------- fit
    def fit(
        self, adata: AnnData, exposure: Exposure, outcomes: np.ndarray, outcome_names: list[str]
    ) -> EstimateTable:
        import torch

        torch.manual_seed(self.seed)
        rng = np.random.default_rng(self.seed)
        obs = adata.obs
        n = adata.n_obs
        bins = np.asarray(exposure.bins_um, dtype=float)
        d_max = float(bins[-1])

        src_np, dst_np, dist = self._edges(adata, d_max)
        ring_of_edge = np.searchsorted(bins, dist, side="left") - 1
        X_np, guide_cols_np, has_call_col = self._features(adata, exposure)

        Y = np.asarray(outcomes, dtype=np.float32)
        mu = Y.mean(axis=0)
        sd = Y.std(axis=0) + 1e-6
        Ys = (Y - mu) / sd

        dev = self._device(n)
        x = torch.as_tensor(X_np, device=dev)
        src = torch.as_tensor(src_np, device=dev)
        dst = torch.as_tensor(dst_np, device=dev)
        e = torch.as_tensor((dist / d_max).astype(np.float32), device=dev)[:, None]
        y = torch.as_tensor(Ys, device=dev)
        guide_cols = torch.as_tensor(guide_cols_np, device=dev)
        model = _make_model(
            X_np.shape[1], len(guide_cols_np), self.hidden, Y.shape[1], self.n_layers
        ).to(dev)
        self._train(model, x, src, dst, e, y, guide_cols, rng)
        with torch.no_grad():
            pred0 = model(x, src, dst, e, guide_cols)

        # one random graph neighbour per (cell, ring); -1 when the cell has none in that ring
        pick = np.full((exposure.n_bins, n), -1, dtype=np.int64)
        order = rng.permutation(len(src_np))
        for b in range(exposure.n_bins):
            m = order[ring_of_edge[order] == b]
            rows_b, first = np.unique(src_np[m], return_index=True)
            pick[b, rows_b] = dst_np[m][first]

        is_ntc = obs["is_ntc"].to_numpy()
        is_pert = obs["is_perturbed"].to_numpy()
        T_auto, rings, C = group_masks(exposure, is_ntc, is_pert, self.groups)
        T_auto_d = T_auto.toarray() > 0
        rings_d = [R.toarray() > 0 for R in rings]
        C_d = C.toarray() > 0
        ct = obs["cell_type"].astype(str).to_numpy()
        groups_out = [*sorted(set(ct)), "all"] if self.report_by_cell_type else ["all"]

        rows: list[dict[str, object]] = []
        for k, target in enumerate(exposure.targets):
            ctrl = np.flatnonzero(C_d[:, k])
            if len(ctrl) > self.max_recipients:
                ctrl = np.sort(rng.choice(ctrl, self.max_recipients, replace=False))
            if len(ctrl) < 2:
                continue
            diffs, used = self._counterfactual(
                model, x, src, dst, e, pred0, ctrl, ctrl, k, has_call_col, guide_cols, rng
            )
            rows += self._rows(
                target,
                "autonomous",
                -1,
                diffs * sd,
                used,
                T_auto_d[:, k],
                ct,
                groups_out,
                outcome_names,
                rng,
            )
            for b in range(exposure.n_bins):
                rec_b = ctrl[pick[b, ctrl] >= 0]
                if len(rec_b) < 2:
                    continue
                diffs, used = self._counterfactual(
                    model,
                    x,
                    src,
                    dst,
                    e,
                    pred0,
                    rec_b,
                    pick[b, rec_b],
                    k,
                    has_call_col,
                    guide_cols,
                    rng,
                )
                rows += self._rows(
                    target,
                    "spillover",
                    b,
                    diffs * sd,
                    used,
                    rings_d[b][:, k],
                    ct,
                    groups_out,
                    outcome_names,
                    rng,
                )
        tab = EstimateTable()
        tab.add(rows)
        return tab
