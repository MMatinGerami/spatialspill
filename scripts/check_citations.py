"""Verify that every DOI mentioned in the repository resolves.

Usage: uv run python scripts/check_citations.py [paths...]
Scans markdown, tex and bib files for ``doi:`` / ``10.xxxx/...`` patterns and queries
https://api.crossref.org/works/<doi>. Dataset DOIs (Zenodo, figshare) are registered with
DataCite, not Crossref, so a Crossref miss falls back to the DataCite API. Exits non-zero if
any DOI fails both.
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import requests

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>\]\)},;]+")
DEFAULT_PATHS = [
    "docs",
    "paper",
    "README.md",
    "NOTEBOOK.md",
    "DECISIONS.md",
    "CRITIQUE.md",
    "data/README.md",
]


def collect(paths: list[str]) -> set[str]:
    dois: set[str] = set()
    for p in paths:
        pp = Path(p)
        files = (
            [pp]
            if pp.is_file()
            else list(pp.rglob("*.md")) + list(pp.rglob("*.tex")) + list(pp.rglob("*.bib"))
        )
        for f in files:
            if not f.exists():
                continue
            for m in DOI_RE.finditer(f.read_text(errors="ignore")):
                dois.add(m.group(0).rstrip(".").rstrip(","))
    return dois


HEADERS = {"User-Agent": "spatialspill-citation-check (mailto:mmatin.gerami@gmail.com)"}


def check(doi: str) -> tuple[bool, str]:
    r = requests.get(f"https://api.crossref.org/works/{doi}", timeout=30, headers=HEADERS)
    if r.status_code == 200:
        title = (r.json()["message"].get("title") or [""])[0]
        return True, title[:90]
    d = requests.get(f"https://api.datacite.org/dois/{doi}", timeout=30, headers=HEADERS)
    if d.status_code == 200:
        titles = d.json()["data"]["attributes"].get("titles") or [{}]
        return True, "[DataCite] " + titles[0].get("title", "")[:80]
    return False, f"HTTP {r.status_code} (Crossref), {d.status_code} (DataCite)"


def main() -> int:
    paths = sys.argv[1:] or DEFAULT_PATHS
    dois = sorted(collect(paths))
    bad = 0
    for d in dois:
        ok, info = check(d)
        print(f"{'OK ' if ok else 'BAD'} {d}  {info}")
        bad += not ok
        time.sleep(0.2)
    print(f"{len(dois)} DOIs checked, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
