"""Manuscript figures and tables, all from files under results/ and reports/.

Usage: uv run python scripts/make_figures.py --fish results/<hash> --multi results/<hash> [--bench results/<hash> --gate results/<hash>]
Writes paper/figures/*.png and paper/tables/*.tex. Every number in the manuscript comes from
these outputs; the LaTeX source only includes them.
"""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "paper" / "figures"
TAB = ROOT / "paper" / "tables"


def fig_clonality() -> None:
    rows = []
    for name in [
        "perturb_fish",
        "perturb_multi",
        "perturb_map",
        "perturb_dbit",
        "spatial_perturbseq",
    ]:
        f = ROOT / "reports" / "audit" / name / "clonality.json"
        a = ROOT / "reports" / "audit" / name / "audit.json"
        if not f.exists() or not a.exists():
            continue
        c = json.loads(f.read_text())
        au = json.loads(a.read_text())
        g = au["graphs"]["delaunay_pruned"]
        rows.append(
            {
                "dataset": name,
                "clonality_z": g["clonality_z"],
                "ratio": g["same_target_edges_observed"]
                / max(g["same_target_edges_null_mean"], 1e-9),
                "frac_in_components": c["same_guide_components"][
                    "fraction_assigned_in_components_ge2"
                ],
                "misassign_cor": (c.get("misassignment_signature") or {}).get("median_cor", np.nan),
            }
        )
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].bar(df["dataset"], df["ratio"], color="#4477aa")
    ax[0].set_yscale("log")
    ax[0].set_ylabel("same-target adjacent pairs, observed / null")
    ax[0].axhline(1, color="k", lw=0.8)
    ax[1].bar(df["dataset"], df["frac_in_components"], color="#ee7733")
    ax[1].set_ylabel("assigned cells in same-guide components")
    ax[2].bar(df["dataset"], df["misassign_cor"], color="#cc3311")
    ax[2].set_ylabel("r(autonomous, unassigned-neighbour profile)")
    for a_ in ax:
        a_.tick_params(axis="x", rotation=30, labelsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_clonality.png", dpi=150)
    df.to_csv(TAB / "clonality.csv", index=False)


def fig_calibration(fish: Path, multi: Path) -> None:
    summ = pd.read_csv(ROOT / "results" / "summary" / "real_data_calibration.csv")
    summ = summ[summ["dataset"].isin(["perturb_fish", "perturb_multi"])]
    summ["label"] = (
        summ["estimator"]
        + np.where(summ["spatial_basis"] > 0, " +basis" + summ["spatial_basis"].astype(str), "")
        + np.where(summ["tile_um"] > 0, " tile" + summ["tile_um"].astype(int).astype(str), "")
    )
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    for ax, ds in zip(axes, ["perturb_fish", "perturb_multi"]):
        s = summ[summ["dataset"] == ds]
        x = np.arange(len(s))
        w = 0.27
        ax.bar(x - w, s["ntc_p05_auto"], w, label="autonomous", color="#cc3311")
        ax.bar(x, s["ntc_p05_r1"], w, label="ring 1", color="#ee7733")
        ax.bar(x + w, s["ntc_p05_r2"], w, label="ring 2", color="#4477aa")
        ax.axhline(0.05, color="k", ls="--", lw=1)
        ax.set_xticks(x)
        ax.set_xticklabels(s["label"], rotation=35, ha="right", fontsize=8)
        ax.set_title(ds)
    axes[0].set_ylabel("NTC pseudo-targets with p < 0.05")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_calibration_real.png", dpi=150)


def fig_benchmark(bench: Path | None, gate: Path | None) -> None:
    frames = []
    for p in (bench, gate):
        if p and (p / "benchmark_summary.csv").exists():
            s = pd.read_csv(p / "benchmark_summary.csv")
            s["run"] = p.name
            frames.append(s)
    if not frames:
        return
    s = pd.concat(frames)
    sp = s[(s["kind"] == "spillover") & (s["ring"] > 0)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, metric, title in zip(
        axes,
        ["ntc_fpr_p05", "power_q", "fdp_q"],
        ["NTC false positives (p<0.05)", "power at q<0.1", "false discovery proportion at q<0.1"],
    ):
        piv = sp.pivot_table(index="scenario", columns="estimator", values=metric, aggfunc="mean")
        piv.plot(kind="bar", ax=ax, width=0.85)
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=30, labelsize=8)
        ax.legend(fontsize=7)
        if metric == "ntc_fpr_p05":
            ax.axhline(0.05, color="k", ls="--", lw=1)
        if metric == "fdp_q":
            ax.axhline(0.1, color="k", ls="--", lw=1)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_benchmark.png", dpi=150)
    agg = (
        sp.groupby(["scenario", "estimator"])[["ntc_fpr_p05", "power_q", "auroc", "fdp_q"]]
        .mean()
        .round(2)
    )
    agg.to_csv(TAB / "benchmark_spillover.csv")
    tex = agg.reset_index().rename(
        columns={
            "ntc_fpr_p05": "NTC FP (p<0.05)",
            "power_q": "power (q<0.1)",
            "auroc": "AUROC",
            "fdp_q": "FDP (q<0.1)",
        }
    )
    (TAB / "benchmark_spillover.tex").write_text(
        tex.to_latex(index=False, na_rep="--", float_format="%.2f", escape=True)
    )


def fig_ranking(fish: Path) -> None:
    rk = pd.read_csv(fish / "niche_ranking.csv")
    fig, ax = plt.subplots(figsize=(7, 6))
    rk = rk.sort_values("score")
    ax.barh(
        rk["target"],
        rk["score"],
        color=np.where(rk["score"] > rk["ntc_score_q95"].iloc[0], "#cc3311", "#bbbbbb"),
    )
    ax.axvline(
        rk["ntc_score_q95"].iloc[0],
        color="k",
        ls="--",
        lw=1,
        label="NTC pseudo-target 95th percentile",
    )
    ax.set_xlabel("niche score (mean z^2 - 1, rings 1 and 2)")
    ax.legend(fontsize=8)
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_ranking.png", dpi=150)
    rk.to_csv(TAB / "niche_ranking_fish.csv", index=False)


def tables(fish: Path, multi: Path) -> None:
    rows = []
    for d in (fish, multi):
        cal = glob.glob(str(d / "*_calibration.json"))
        af = d / "artifact_fraction.json"
        pr = d / "prediction.json"
        r = json.loads(Path(cal[0]).read_text()) if cal else {}
        a = json.loads(af.read_text()) if af.exists() else {}
        p = json.loads(pr.read_text()) if pr.exists() else {}
        rows.append(
            {
                "dataset": r.get("dataset"),
                "run": d.name,
                "n_tests_all": r.get("n_tests"),
                "ntc_p05_auto": r.get("ntc_p05_autonomous"),
                "ntc_p05_r1": r.get("ntc_p05_ring1"),
                "ntc_p05_r2": r.get("ntc_p05_ring2"),
                "auto_hits": r.get("n_autonomous_hits"),
                "auto_fdp_ntc": r.get("ntc_fdp_autonomous"),
                "spill_hits": r.get("n_spillover_hits"),
                "spill_fdp_ntc": r.get("ntc_fdp_spillover"),
                "artifact_r2_ring0_minus_last": a.get("r2_ring0_minus_last"),
                "alpha_ring0": a.get("alpha_ring0"),
                **{f"pred_{k}": v.get("mean_r") for k, v in p.items()},
                **{f"prednull_{k}": v.get("null_mean_r") for k, v in p.items()},
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "real_data_summary.csv", index=False)
    (TAB / "real_data_summary.tex").write_text(df.round(3).to_latex(index=False))


def fig_power_ablation_sibling() -> None:
    """Second-pass outputs: power grid figure, ablation table, sibling benchmark table."""
    runs = sorted(RES.glob("*/"), key=lambda d: d.stat().st_mtime)
    pw = [d for d in runs if (d / "power_summary.csv").exists()]
    if pw:
        import shutil

        shutil.copy(pw[-1] / "power_grid.png", FIG / "fig6_power.png")
        shutil.copy(pw[-1] / "power_summary.csv", TAB / "power_summary.csv")
    ab = [
        d
        for d in runs
        if (
            (d / "benchmark_summary.csv").exists()
            and "benchmark_ablation" in (d / "config.yaml").read_text()
        )
        or ((d / "config.yaml").exists() and "E1_sample_strata" in (d / "config.yaml").read_text())
    ]
    if ab:
        s = pd.read_csv(ab[-1] / "benchmark_summary.csv")
        sp = s[(s["kind"] == "spillover") & (s["ring"] > 0)]
        agg = (
            sp.groupby(["scenario", "estimator"])[["ntc_fpr_p05", "power_q", "fdp_q", "n_tests"]]
            .mean()
            .round(2)
        )
        agg.to_csv(TAB / "ablation.csv")
        (TAB / "ablation.tex").write_text(
            agg.reset_index()
            .rename(
                columns={
                    "ntc_fpr_p05": "NTC FP (p<0.05)",
                    "power_q": "power (q<0.1)",
                    "fdp_q": "FDP (q<0.1)",
                    "n_tests": "tests",
                }
            )
            .to_latex(index=False, na_rep="--", float_format="%.2f", escape=True)
        )
    sb = [d for d in runs if (d / "sibling_benchmark.csv").exists()]
    if sb:
        df = pd.read_csv(sb[-1] / "sibling_benchmark.csv")
        agg = (
            df.groupby(["scenario", "detection"])[
                ["pi", "r_unc_spill", "r_cor_spill", "r_unc_auto", "r_cor_auto"]
            ]
            .mean()
            .round(2)
        )
        agg.to_csv(TAB / "sibling.csv")
        (TAB / "sibling.tex").write_text(
            agg.reset_index()
            .rename(
                columns={
                    "pi": "pi hat",
                    "r_unc_spill": "r(unc., spill)",
                    "r_cor_spill": "r(corr., spill)",
                    "r_unc_auto": "r(unc., auto)",
                    "r_cor_auto": "r(corr., auto)",
                }
            )
            .to_latex(index=False, float_format="%.2f", escape=True)
        )


def numbers(fish: Path, multi: Path, gate: Path | None) -> None:
    """LaTeX macros for every number quoted in the manuscript prose."""
    macros: dict[str, object] = {}

    def add(name: str, value: object, fmt: str = "{:.3f}") -> None:
        if value is None or (isinstance(value, float) and not np.isfinite(value)):
            macros[name] = "n/a"
        elif isinstance(value, int | np.integer):
            macros[name] = f"{int(value):,}"
        else:
            macros[name] = fmt.format(value)

    # exploratory pixel-level E2 run on Perturb-DBiT (unperturbed-pixel controls, A4)
    dbit = sorted(RES.glob("*/perturb_dbit_calibration.json"), key=lambda q: q.stat().st_mtime)
    r = json.loads(dbit[-1].read_text()) if dbit else {}
    add("DbitNTests", int(r.get("n_tests", 0)))
    add("DbitNtcP", r.get("ntc_fraction_p_below_0.05"))
    add("DbitNtcQ", r.get("ntc_fraction_q_below_fdr"))
    add("DbitAutoFdp", r.get("ntc_fdp_autonomous"), "{:.2f}")
    add("DbitSpillFdp", r.get("ntc_fdp_spillover"), "{:.2f}")

    for tag, d in (("Fish", fish), ("Multi", multi)):
        cal = glob.glob(str(d / "*_calibration.json"))
        r = json.loads(Path(cal[0]).read_text()) if cal else {}
        add(f"{tag}NtcPAuto", r.get("ntc_p05_autonomous"))
        add(f"{tag}NtcPRingZero", r.get("ntc_p05_ring0"))
        add(f"{tag}NtcPRingOne", r.get("ntc_p05_ring1"))
        add(f"{tag}NtcPRingTwo", r.get("ntc_p05_ring2"))
        add(f"{tag}NtcQAuto", r.get("ntc_q_autonomous"))
        add(f"{tag}NtcQSpill", r.get("ntc_q_spillover"))
        add(f"{tag}AutoHits", int(r.get("n_autonomous_hits", 0)))
        add(f"{tag}SpillHits", int(r.get("n_spillover_hits", 0)))
        add(f"{tag}AutoFdp", r.get("ntc_fdp_autonomous"), "{:.2f}")
        add(f"{tag}SpillFdp", r.get("ntc_fdp_spillover"), "{:.2f}")
        add(f"{tag}NTests", int(r.get("n_tests", 0)))
        af = d / "artifact_fraction.json"
        a = json.loads(af.read_text()) if af.exists() else {}
        add(f"{tag}ArtifactFrac", a.get("r2_ring0_minus_last"))
        add(f"{tag}AlphaRingZero", a.get("alpha_ring0"))
        add(f"{tag}ArtifactTargets", int(a.get("n_targets", 0)))
        pr = d / "prediction.json"
        pz = json.loads(pr.read_text()) if pr.exists() else {}
        for k, v in pz.items():
            key = (
                k.replace("replogle_k562_gwps", "Rep")
                .replace("go", "Go")
                .replace(":", "")
                .replace("autonomous", "Auto")
                .replace("spillover", "Spill")
            )
            add(f"{tag}Pred{key}", v.get("mean_r"), "{:.2f}")
            add(f"{tag}PredNull{key}", v.get("null_mean_r"), "{:.2f}")
        rk = d / "niche_ranking.csv"
        if rk.exists():
            df = pd.read_csv(rk)
            add(f"{tag}TopTarget", str(df.iloc[0]["target"]), "{}")
            add(f"{tag}TopScore", float(df.iloc[0]["score"]), "{:.2f}")
            add(f"{tag}NtcScoreQ", float(df["ntc_score_q95"].iloc[0]), "{:.2f}")
            add(f"{tag}NAboveNtc", int((df["score"] > df["ntc_score_q95"]).sum()))
            add(f"{tag}NRanked", len(df))
    # Runs quoted in the prose, pinned by hash (NOTEBOOK 2026-10-05). Picking the last row per
    # setting took a bin-sensitivity rerun instead of the pre-registered primary tile run.
    summ = pd.read_csv(ROOT / "results" / "summary" / "real_data_calibration.csv")
    pinned = {
        "EOneGlobal": "5d30bb4682",  # Perturb-FISH E1, global strata
        "EOneTile": "8c106a7091",  # Perturb-FISH E1, 250 um tiles, bins 0/15/30/60 (primary)
        "FishAnalytic": "246b4deb3e",  # Perturb-FISH E2, 40-centre basis, analytic SEs
    }
    for tag, run in pinned.items():
        r = summ[(summ["run"] == run) & (summ["dataset"] == "perturb_fish")]
        if len(r):
            r = r.iloc[0]
            add(f"{tag}NTests", int(r["n_tests"]))
            add(f"{tag}NtcP", float(r["ntc_p05"]))
            add(
                f"{tag}NtcPRingZero",
                float(r["ntc_p05_r0"])
                if "ntc_p05_r0" in r and np.isfinite(r["ntc_p05_r0"])
                else None,
            )
            add(f"{tag}AutoHits", int(r["autonomous_hits"]))
            add(f"{tag}SpillHits", int(r["spillover_hits"]))
            for col, name in (
                ("ntc_p05_auto", "Auto"),
                ("ntc_p05_r1", "RingOne"),
                ("ntc_p05_r2", "RingTwo"),
            ):
                v = float(r[col]) if np.isfinite(r[col]) else None
                add(f"{tag}NtcP{name}", v)
    m1 = summ[(summ["run"] == "53c85fb53d") & (summ["dataset"] == "perturb_multi")]  # E1 global
    if len(m1):
        r = m1.iloc[0]
        add("MultiEOneNtcP", float(r["ntc_p05"]))
        add("MultiEOneNTests", int(r["n_tests"]))
    ht = ROOT / "results" / "summary" / "heldout_technology.json"
    if ht.exists():
        h = json.loads(ht.read_text())
        add("HeldoutRho", h["spearman"], "{:.2f}")
        add("HeldoutNull", h["null_mean"], "{:.2f}")
        add("HeldoutP", h["empirical_p"], "{:.2f}")
        add("HeldoutOverlap", int(h["n_overlapping_targets"]))
    cl = TAB / "clonality.csv"
    if cl.exists():
        c = pd.read_csv(cl).set_index("dataset")
        for ds, tag in (
            ("perturb_fish", "Fish"),
            ("perturb_multi", "Multi"),
            ("perturb_map", "Map"),
            ("perturb_dbit", "Dbit"),
            ("spatial_perturbseq", "Stereo"),
        ):
            if ds in c.index:
                add(f"{tag}ClonRatio", float(c.loc[ds, "ratio"]), "{:.0f}")
                add(f"{tag}ClonZ", float(c.loc[ds, "clonality_z"]), "{:.0f}")
                add(f"{tag}MisCor", float(c.loc[ds, "misassign_cor"]), "{:.2f}")
    for tag in ("EOneTile", "ETwoPerm"):
        for sdtag in ("Half", "One", "Two"):
            for suffix in ("MinRecipients", "MaxPower", "MaxRecipients"):
                macros.setdefault(f"Power{tag}{sdtag}{suffix}", "n/a")
    pw = TAB / "power_summary.csv"
    if pw.exists():
        ps = pd.read_csv(pw)
        for est, tag in (("E1_tile", "EOneTile"), ("E2_spatial_perm", "ETwoPerm")):
            sub = ps[ps["estimator"] == est]
            for sd, sdtag in ((0.5, "Half"), (1.0, "One"), (2.0, "Two")):
                ss = sub[sub["spill_lfc_sd"] == sd]
                if len(ss):
                    ok = ss[ss["power"] >= 0.8]
                    add(
                        f"Power{tag}{sdtag}MinRecipients",
                        float(ok["ntc_recipients"].min()) if len(ok) else None,
                        "{:.0f}",
                    )
                    add(f"Power{tag}{sdtag}MaxPower", float(ss["power"].max()), "{:.2f}")
                    add(
                        f"Power{tag}{sdtag}MaxRecipients",
                        float(ss["ntc_recipients"].max()),
                        "{:.0f}",
                    )
    sbf = TAB / "sibling.csv"
    if sbf.exists():
        sb = pd.read_csv(sbf)
        row = sb[(sb["scenario"] == "clonal") & (sb["detection"] == 0.4)]
        if len(row):
            r = row.iloc[0]
            add("SibPiClonalLow", float(r["pi"]), "{:.2f}")
            add("SibRUncAuto", float(r["r_unc_auto"]), "{:.2f}")
            add("SibRCorAuto", float(r["r_cor_auto"]), "{:.2f}")
            add("SibRUncSpill", float(r["r_unc_spill"]), "{:.2f}")
            add("SibRCorSpill", float(r["r_cor_spill"]), "{:.2f}")
    sr = sorted(RES.glob("*/sibling_real.json"), key=lambda f: f.stat().st_mtime)
    if sr:
        j = json.loads(sr[-1].read_text())
        add("SibRealNTargets", int(j.get("n_targets", 0)))
        add("SibRealPiMedian", j.get("pi_median"), "{:.2f}")
        iqr = j.get("pi_iqr")
        macros["SibRealPiIqr"] = f"{iqr[0]:.2f} to {iqr[1]:.2f}" if iqr else "n/a"
        add("SibRealAgreeCor", j.get("sign_agreement_corrected"), "{:.2f}")
        add("SibRealAgreeUnc", j.get("sign_agreement_uncorrected"), "{:.2f}")
        add("SibRealRCor", j.get("r_corrected"), "{:.2f}")
    lines = [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in macros.items()]
    (TAB / "numbers.tex").write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fish", required=True)
    ap.add_argument("--multi", required=True)
    ap.add_argument("--bench", default=None)
    ap.add_argument("--gate", default=None)
    a = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    fish, multi = Path(a.fish), Path(a.multi)
    fig_clonality()
    fig_calibration(fish, multi)
    fig_benchmark(Path(a.bench) if a.bench else None, Path(a.gate) if a.gate else None)
    fig_ranking(fish)
    tables(fish, multi)
    fig_power_ablation_sibling()
    numbers(fish, multi, Path(a.gate) if a.gate else None)
    print("figures in", FIG, "tables in", TAB)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
