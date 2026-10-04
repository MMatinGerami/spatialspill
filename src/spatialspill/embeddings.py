"""Gene embeddings for held-out-gene prediction (C5).

Three sources, all reduced to a dense (genes x d) matrix keyed by gene symbol:

- ``perturbseq_embedding``: Replogle et al. 2022 pseudobulk profiles (perturbations x genes).
  Each perturbed gene is embedded by the PCA scores of its mean expression profile relative to
  control; genes not perturbed in the screen have no embedding (missing by design).
- ``go_embedding``: binary gene x GO-term matrix from a GAF file (propagated to ancestors with
  the OBO ontology when given), reduced by truncated SVD.
- ``esm2_embedding``: mean-pooled ESM2 residue embeddings of the canonical protein sequence
  (sequences fetched from UniProt REST for the requested symbols; cached on disk).

All functions return ``pandas.DataFrame`` indexed by upper-cased symbol so that human and mouse
orthologues with the same symbol align; this is a deliberate simplification recorded in the
limitations.
"""

from __future__ import annotations

import gzip
import json
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, TruncatedSVD


def perturbseq_embedding(
    bulk_h5ad: str | Path, n_components: int = 32, control_pattern: str = "non-targeting"
) -> pd.DataFrame:
    """PCA embedding of perturbation pseudobulk profiles (rows: perturbed genes)."""
    import anndata as ad

    a = ad.read_h5ad(bulk_h5ad)
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    X = np.asarray(X, dtype=np.float64)
    names = pd.Series(a.obs_names.astype(str))
    # Replogle bulk files name perturbations "<id>_<GENE>_<guides>_<ENSG>"; the normalised
    # tables are already centred on non-targeting controls (which are not rows)
    parts = names.str.split("_")
    genes = np.array(
        [p[1] if (len(p) > 1 and p[0].isdigit()) else p[0] for p in parts], dtype=object
    )
    genes = np.char.upper(genes.astype(str))
    ctrl = names.str.contains(control_pattern, case=False).to_numpy()
    if ctrl.any():
        X = X - X[ctrl].mean(0, keepdims=True)
    keep = ~ctrl
    n_comp = min(n_components, keep.sum() - 1, X.shape[1])
    Z = PCA(n_components=n_comp, random_state=0).fit_transform(X[keep])
    df = pd.DataFrame(Z, index=genes[keep], columns=[f"pc{i}" for i in range(n_comp)])
    return df[~df.index.duplicated()]


def read_gaf(path: str | Path, aspects: Iterable[str] = ("P", "F", "C")) -> pd.DataFrame:
    """Gene symbol x GO term binary matrix from a GAF 2.x file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    pairs: set[tuple[str, str]] = set()
    aspects = set(aspects)
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[8] not in aspects:
                continue
            if f[3].startswith("NOT"):
                continue
            pairs.add((f[2].upper(), f[4]))
    df = pd.DataFrame(sorted(pairs), columns=["gene", "term"])
    df["v"] = 1
    return df.pivot_table(index="gene", columns="term", values="v", fill_value=0).astype(np.int8)


def read_obo_parents(path: str | Path) -> dict[str, list[str]]:
    parents: dict[str, list[str]] = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                cur = None
            elif line.startswith("id: GO:"):
                cur = line[4:]
                parents.setdefault(cur, [])
            elif cur and line.startswith("is_a: "):
                parents[cur].append(line[6:].split(" ")[0])
            elif cur and line.startswith("relationship: part_of "):
                parents[cur].append(line.split(" ")[2])
    return parents


def propagate(matrix: pd.DataFrame, parents: dict[str, list[str]]) -> pd.DataFrame:
    """Add ancestor terms (is_a, part_of) to every annotation."""
    anc_cache: dict[str, set[str]] = {}

    def ancestors(t: str) -> set[str]:
        if t in anc_cache:
            return anc_cache[t]
        out: set[str] = set()
        stack = list(parents.get(t, []))
        while stack:
            p = stack.pop()
            if p not in out:
                out.add(p)
                stack.extend(parents.get(p, []))
        anc_cache[t] = out
        return out

    cols = list(matrix.columns)
    all_terms = set(cols)
    for t in cols:
        all_terms |= ancestors(t)
    all_terms_l = sorted(all_terms)
    tidx = {t: i for i, t in enumerate(all_terms_l)}
    M = np.zeros((matrix.shape[0], len(all_terms_l)), dtype=np.int8)
    X = matrix.to_numpy()
    for j, t in enumerate(cols):
        rows = np.flatnonzero(X[:, j])
        M[rows, tidx[t]] = 1
        for a in ancestors(t):
            M[rows, tidx[a]] = 1
    return pd.DataFrame(M, index=matrix.index, columns=all_terms_l)


def go_embedding(
    gaf_path: str | Path,
    obo_path: str | Path | None = None,
    n_components: int = 32,
    min_genes_per_term: int = 5,
) -> pd.DataFrame:
    M = read_gaf(gaf_path)
    if obo_path is not None:
        M = propagate(M, read_obo_parents(obo_path))
    keep_terms = np.asarray(M.to_numpy().sum(axis=0)) >= min_genes_per_term
    M = M.iloc[:, np.flatnonzero(keep_terms)]
    n_comp = min(n_components, M.shape[1] - 1)
    Z = TruncatedSVD(n_components=n_comp, random_state=0).fit_transform(
        M.to_numpy(dtype=np.float32)
    )
    return pd.DataFrame(Z, index=M.index, columns=[f"svd{i}" for i in range(n_comp)])


def fetch_uniprot_sequences(
    symbols: Iterable[str], organism_id: int, cache: str | Path
) -> dict[str, str]:
    """Canonical reviewed sequences from UniProt REST, cached as JSON. Network access required."""
    import requests

    cache = Path(cache)
    seqs: dict[str, str] = json.loads(cache.read_text()) if cache.exists() else {}
    for s in symbols:
        key = f"{organism_id}:{s.upper()}"
        if key in seqs:
            continue
        q = f"gene_exact:{s} AND organism_id:{organism_id} AND reviewed:true"
        params: dict[str, str] = {"query": q, "fields": "sequence", "format": "json", "size": "1"}
        r = requests.get("https://rest.uniprot.org/uniprotkb/search", params=params, timeout=60)
        r.raise_for_status()
        res = r.json().get("results", [])
        seqs[key] = res[0]["sequence"]["value"] if res else ""
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(seqs))
    return {k.split(":", 1)[1]: v for k, v in seqs.items() if k.startswith(f"{organism_id}:")}


def esm2_embedding(
    sequences: dict[str, str], model_dir: str | Path, max_len: int = 1022, device: str = "cpu"
) -> pd.DataFrame:
    """Mean-pooled last-layer ESM2 embeddings for each sequence (keyed by symbol)."""
    import torch
    from transformers import AutoModel, AutoTokenizer

    tok = AutoTokenizer.from_pretrained(str(model_dir))
    model = AutoModel.from_pretrained(str(model_dir)).to(device).eval()
    rows, idx = [], []
    with torch.no_grad():
        for sym, seq in sequences.items():
            if not seq:
                continue
            enc = tok(seq[:max_len], return_tensors="pt").to(device)
            h = model(**enc).last_hidden_state[0, 1:-1]
            rows.append(h.mean(0).cpu().numpy())
            idx.append(sym.upper())
    return (
        pd.DataFrame(
            np.stack(rows), index=idx, columns=[f"esm{i}" for i in range(rows[0].shape[0])]
        )
        if rows
        else pd.DataFrame()
    )
