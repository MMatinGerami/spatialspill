#!/usr/bin/env bash
# Perturb-FISH (Binan et al., Cell 2025) processed tables from the Brain Image Library
# (doi:10.35077/ace-gem-get, CC BY 4.0). Only final tables and cell metadata are fetched;
# raw images (about 26 TB) are not.
set -euo pipefail
B=https://download.brainimagelibrary.org/0c/bd/0cbd479c521afff9
OUT="$(dirname "$0")/../data/raw/perturb_fish"
mkdir -p "$OUT"
cd "$OUT"
curl -sS -L -C - -o README.txt "$B/astrocytes/README.txt"
for sub in extras/tumors/processed/finaltables extras/THP1/processed/finaltables astrocytes/processed/finaltables; do
  wget -q -r -np -nH --cut-dirs=3 -R "index.html*" -c "$B/$sub/" || echo "wget failed for $sub"
done
mkdir -p tumors/cell_metadata astrocytes/cell_metadata
curl -sS -L -C - -o tumors/cell_metadata/cell_metadata.csv "$B/extras/tumors/processed/sample1/MERFISH/DeriveEntityMetadataTask/region_0/cell_metadata.csv"
curl -sS -L -C - -o astrocytes/cell_metadata/cell_metadata_region0.csv "$B/astrocytes/processed/sample1and2/MERFISH/DeriveEntityMetadataTask/region_0/cell_metadata.csv"
curl -sS -L -C - -o astrocytes/cell_metadata/cell_metadata_region1.csv "$B/astrocytes/processed/sample1and2/MERFISH/DeriveEntityMetadataTask/region_1/cell_metadata.csv"
find . -type f | sort | xargs ls -la
echo DONE_FISH
