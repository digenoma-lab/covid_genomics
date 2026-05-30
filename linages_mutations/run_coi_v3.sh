python call_coinfections_v3.py \
  --variants-dir iSNVs \
  --lineage-defs lineage_defs_v3.tsv \
  --glob "*.tsv" \
  --outdir coinfection_calls_v3_group \
  --label-mode group \
  --min-total-dp 100 \
  --min-alt-dp 50 \
  --major-af-min 0.70 \
  --major-min-supported-sites 6 \
  --major-min-frac 0.30 \
  --minor-af-min 0.08 \
  --minor-af-max 0.40 \
  --min-minor-support-sites 4 \
  --min-minor-median-af 0.10 \
  --max-lineages-per-site 5  



python call_coinfections_v3.py \
  --variants-dir iSNVs \
  --lineage-defs lineage_defs_v3.tsv \
  --glob "*.tsv" \
  --outdir coinfection_calls_lineage_v3 \
  --label-mode lineage \
  --min-total-dp 100 \
  --min-alt-dp 50 \
  --major-af-min 0.70 \
  --major-min-supported-sites 8 \
  --major-min-frac 0.40 \
  --minor-af-min 0.10 \
  --minor-af-max 0.40 \
  --min-minor-support-sites 5 \
  --min-minor-median-af 0.12 \
  --max-lineages-per-site 5
