#!/usr/bin/env bash
# Supporting resources for Phase 4 and 5 (see data/README.md): Replogle 2022 pseudobulk tables
# (figshare, CC BY 4.0), GO annotations and ontology (CC BY 4.0), OmniPath ligand-receptor
# tables, CellPhoneDB v5 interaction inputs, ESM2 t12 35M checkpoint (MIT, via Hugging Face).
set -uo pipefail
OUT="$(cd "$(dirname "$0")/.." && pwd)/data/raw/supporting"
mkdir -p "$OUT"/replogle "$OUT"/go "$OUT"/lr "$OUT"/esm2
cd "$OUT"
curl -sS -L -C - -o replogle/K562_gwps_normalized_bulk_01.h5ad https://ndownloader.figshare.com/files/35773217
curl -sS -L -C - -o replogle/K562_essential_normalized_bulk_01.h5ad https://ndownloader.figshare.com/files/35780870
curl -sS -L -C - -o replogle/rpe1_normalized_bulk_01.h5ad https://ndownloader.figshare.com/files/35775512
curl -sS -L -C - -o go/goa_human.gaf.gz https://current.geneontology.org/annotations/goa_human.gaf.gz
curl -sS -L -C - -o go/goa_mouse.gaf.gz https://current.geneontology.org/annotations/mgi.gaf.gz
curl -sS -L -C - -o go/go-basic.obo https://current.geneontology.org/ontology/go-basic.obo
curl -sS -L -o lr/omnipath_intercell.tsv "https://omnipathdb.org/intercell?resources=CellPhoneDB,CellChatDB,ICELLNET,Guide2Pharma,connectomeDB2020&format=tsv"
curl -sS -L -o lr/omnipath_ligrec.tsv "https://omnipathdb.org/interactions?datasets=ligrecextra&format=tsv&fields=sources,references,curation_effort&genesymbols=1"
curl -sS -L -o lr/cellphonedb_interaction_input.csv https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv
curl -sS -L -o lr/cellphonedb_gene_input.csv https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/gene_input.csv
curl -sS -L -o lr/cellphonedb_complex_input.csv https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/complex_input.csv
for f in config.json model.safetensors tokenizer_config.json vocab.txt special_tokens_map.json; do
  curl -sS -L -C - -o esm2/$f https://huggingface.co/facebook/esm2_t12_35M_UR50D/resolve/main/$f
done
ls -la replogle go lr esm2
echo DONE_SUPPORT
