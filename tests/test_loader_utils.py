import gzip

import numpy as np
import pandas as pd

from spatialspill.loaders._util import dominant_guide, guide_to_target, is_ntc_guide
from spatialspill.loaders.perturb_dbit import _parse_pixel, read_expression, read_sg_table


def test_guide_to_target_formats():
    assert guide_to_target("Bmi1.g2") == "Bmi1"
    assert guide_to_target("Rps19_sg1") == "Rps19"
    assert guide_to_target("sg54894_SLC22A16") == "SLC22A16"
    assert guide_to_target("sgrna12550_Pax8") == "Pax8"
    assert guide_to_target("NTC2") == "NTC"
    assert guide_to_target("control_gBC_00009") == "NTC"
    assert is_ntc_guide("Non-targeting_1") and not is_ntc_guide("Nsd1.g2")


def test_dominant_guide_rules():
    df = pd.DataFrame(
        {
            "pixel": ["1x1", "1x1", "1x1", "2x2", "3x3", "3x3"],
            "grna": ["A.g1", "A.g1", "B.g1", "C.g1", "D.g1", "E.g1"],
            "count": [3, 2, 1, 1, 1, 1],
        }
    )
    out = dominant_guide(df, "pixel", "grna", "count", min_count=1, min_frac=0.5)
    assert out.loc["1x1", "guide"] == "A.g1"
    assert np.isclose(out.loc["1x1", "guide_confidence"], 5 / 6)
    assert out.loc["2x2", "guide"] == "C.g1"
    assert out.loc["3x3", "guide"] == "none"  # tie at 50% is accepted only if >= min_frac
    assert out.loc["3x3", "n_guides"] == 2
    out2 = dominant_guide(df, "pixel", "grna", "count", min_count=2)
    assert out2.loc["2x2", "guide"] == "none"


def test_read_sg_table_three_formats(tmp_path):
    p1 = tmp_path / "a.txt.gz"
    with gzip.open(p1, "wt") as fh:
        fh.write("BA\tBB\tgrna\tcount\n1\t2\tX.g1\t4\n")
    p2 = tmp_path / "b.txt.gz"
    with gzip.open(p2, "wt") as fh:
        fh.write('"","barcode","grna","count"\n"1","10x12","sgrna1_Y",1\n')
    p3 = tmp_path / "c.txt.gz"
    with gzip.open(p3, "wt") as fh:
        fh.write("bc10x14\tZ_sg1\t1\t\n")
    for p, pix, g in [(p1, "1x2", "X.g1"), (p2, "10x12", "sgrna1_Y"), (p3, "10x14", "Z_sg1")]:
        df = read_sg_table(p)
        assert list(df.columns) == ["pixel", "grna", "count"]
        assert df.iloc[0]["pixel"] == pix and df.iloc[0]["grna"] == g


def test_read_expression_both_orientations(tmp_path):
    a = tmp_path / "pix_by_gene.tsv.gz"
    with gzip.open(a, "wt") as fh:
        fh.write("\tGeneA\tGeneB\n1x1\t1\t0\n2x1\t0\t3\n")
    b = tmp_path / "gene_by_pix.tsv.gz"
    with gzip.open(b, "wt") as fh:
        fh.write("\t1x1\t2x1\nX4933401J01Rik\t1\t0\nGeneB\t0\t3\n")
    da = read_expression(a)
    db = read_expression(b)
    assert list(da.index) == ["1x1", "2x1"] and list(db.index) == ["1x1", "2x1"]
    assert "4933401J01Rik" in db.columns  # leading X stripped from numeric-start gene names
    assert _parse_pixel("bc10x14") == (10, 14) and _parse_pixel("GeneA") is None
