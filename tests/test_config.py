from omegaconf import OmegaConf

from spatialspill.config import config_hash, load_config, results_dir


def test_hash_deterministic_and_order_invariant():
    a = {"x": 1, "y": {"z": 2}}
    b = {"y": {"z": 2}, "x": 1}
    assert config_hash(a) == config_hash(b)
    assert config_hash(a) != config_hash({"x": 2, "y": {"z": 2}})
    assert config_hash(OmegaConf.create(a)) == config_hash(a)


def test_load_and_results_dir(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("seed: 0\nmodel:\n  k: 3\n")
    cfg = load_config(cfg_file, ["model.k=4"])
    assert cfg.model.k == 4
    d = results_dir(cfg, root=tmp_path / "results")
    assert (d / "config.yaml").exists()
    assert d.name == config_hash(cfg)
