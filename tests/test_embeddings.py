import gzip

import numpy as np

from spatialspill.embeddings import go_embedding, propagate, read_gaf, read_obo_parents


def test_gaf_obo_propagation_and_svd(tmp_path):
    gaf = tmp_path / "a.gaf.gz"
    lines = ["!gaf-version: 2.2"]
    genes = [f"G{i}" for i in range(12)]
    for i, g in enumerate(genes):
        term = "GO:0000002" if i % 2 else "GO:0000003"
        lines.append(
            "\t".join(
                [
                    "UniProtKB",
                    f"P{i}",
                    g,
                    "",
                    term,
                    "PMID:1",
                    "IDA",
                    "",
                    "P",
                    "",
                    "",
                    "protein",
                    "taxon:9606",
                    "20200101",
                    "UniProt",
                ]
            )
        )
    lines.append(
        "\t".join(
            [
                "UniProtKB",
                "P99",
                "G0",
                "NOT",
                "GO:0000003",
                "PMID:1",
                "IDA",
                "",
                "P",
                "",
                "",
                "protein",
                "taxon:9606",
                "20200101",
                "UniProt",
            ]
        )
    )
    with gzip.open(gaf, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    obo = tmp_path / "go.obo"
    obo.write_text(
        "[Term]\nid: GO:0000001\nname: root\n\n[Term]\nid: GO:0000002\nis_a: GO:0000001 ! root\n\n[Term]\nid: GO:0000003\nrelationship: part_of GO:0000001 ! root\n"
    )
    M = read_gaf(gaf)
    assert M.shape == (12, 2) and M.loc["G0", "GO:0000003"] == 1
    P = propagate(M, read_obo_parents(obo))
    assert "GO:0000001" in P.columns and (P["GO:0000001"] == 1).all()
    E = go_embedding(gaf, obo, n_components=2, min_genes_per_term=3)
    assert E.shape[0] == 12 and np.isfinite(E.to_numpy()).all()
