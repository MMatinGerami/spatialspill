"""Configuration hashing so results are keyed by the exact config that produced them."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import yaml
from omegaconf import DictConfig, OmegaConf


def to_plain(cfg: DictConfig | dict[str, Any]) -> dict[str, Any]:
    if isinstance(cfg, DictConfig):
        out = OmegaConf.to_container(cfg, resolve=True)
        assert isinstance(out, dict)
        return {str(k): v for k, v in out.items()}
    return dict(cfg)


def config_hash(cfg: DictConfig | dict[str, Any], length: int = 10) -> str:
    """Deterministic short hash of a config (sorted-key JSON, sha256)."""
    plain = to_plain(cfg)
    s = json.dumps(plain, sort_keys=True, default=str)
    return hashlib.sha256(s.encode()).hexdigest()[:length]


def load_config(path: str | Path, overrides: list[str] | None = None) -> DictConfig:
    base = OmegaConf.load(path)
    assert isinstance(base, DictConfig)
    if overrides:
        base = OmegaConf.merge(base, OmegaConf.from_dotlist(overrides))
        assert isinstance(base, DictConfig)
    return base


def results_dir(cfg: DictConfig | dict[str, Any], root: str | Path = "results") -> Path:
    """``results/<hash>/`` for this config; writes ``config.yaml`` and ``code_version.txt`` (git revision) into it.

    The hash covers the configuration only, so rerunning the same config after a code change
    overwrites the directory; ``code_version.txt`` records which code produced the files."""
    h = config_hash(cfg)
    d = Path(root) / h
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "config.yaml", "w") as fh:
        yaml.safe_dump(to_plain(cfg), fh, sort_keys=True)
    try:
        import subprocess

        rev = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "src", "scripts"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip()
        (d / "code_version.txt").write_text(f"{rev}{' dirty' if dirty else ''}\n")
    except OSError:
        pass
    return d
