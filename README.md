# covid_genomics

Analysis repository for SARS-CoV-2 genomes from Chile, with emphasis on
consensus SNVs, intra-host SNVs (iSNVs), lineage/variant comparisons, and
coinfection detection.

## Main analysis files

- `iSNV_COI_V3.Rmd`: main current iSNV and coinfection analysis.
  This notebook loads iSNV and SNV calls, applies primary and sensitivity iSNV
  filters, joins Pangolin variant calls and coinfection annotations, builds QC
  plots, compares iSNV vs SNV burden, fits regression models, and runs dN/dS
  analyses.
- `analysis/MFig_Covid.Rmd`: manuscript-style figure analysis. It contains the
  original map, demographics, GISAID context, phylogeny, coverage, SNV/iSNV,
  regression, statistical testing, and dN/dS figure workflows.
- `linages_mutations/call_coinfections_v3.py`: coinfection caller based on
  lineage-defining mutations. It scores major and minor lineage support per
  sample and writes cohort-level and per-sample outputs.


## Key inputs

The current iSNV workflow expects these inputs to exist relative to the
repository root:

- `iSNVs/*_variants.tsv`
- `SNVs/*.variants.tsv`
- `demographics/all_samples.txt`
- `analysis/all_results.csv`
- `analysis/resume_samtools_coverage.txt`
- `linages_mutations/coinfection_calls_lineage_v3/cohort_coinfection_summary.tsv`

`analysis/all_results.csv` provides the Pangolin/variant label used for each
sample. The coinfection summary is used only to flag coinfection status and
attach coinfection-support metadata.

## iSNV analysis

Run or knit:

```r
rmarkdown::render("iSNV_COI_V3.Rmd")
```

The notebook defines two iSNV call sets:

- `isnv_primary`: stricter call set for downstream interpretation.
- `isnv_sensitivity`: more permissive call set used to assess robustness.

The workflow includes:

- raw iSNV loading with an RDS cache at `analysis/raw_isnv_cache.rds`;
- filtering by allele frequency, quality, strand support, Ct, genome ends,
  consensus SNVs, and recurrent suspicious hotspots;
- comparison of `isnv_primary` vs `isnv_sensitivity`;
- sample-level iSNV burden summaries;
- SNV vs iSNV burden over time;
- regression models and non-parametric tests;
- variant/coinfection group comparisons;
- dN/dS analyses.

## Coinfection calling

The coinfection caller is:

```bash
python linages_mutations/call_coinfections_v3.py \
  --variants-dir <per-sample-variant-tables> \
  --lineage-defs <lineage-definitions.tsv> \
  --outdir linages_mutations/coinfection_calls_lineage_v3 \
  --label-mode lineage
```

Important options include:

- `--label-mode lineage` or `--label-mode group`
- `--min-total-dp`
- `--min-alt-dp`
- `--major-af-min`
- `--major-min-supported-sites`
- `--minor-af-min`
- `--minor-af-max`
- `--min-minor-support-sites`
- `--min-minor-median-af`

Main outputs:

- `cohort_coinfection_summary.tsv`
- `*.major_scores.tsv`
- `*.minor_scores.tsv`
- `*.filtered_variants.tsv`
- `lineage_defs.filtered.tsv`

## Notes

- Paths in the Rmds are mostly relative to the repository root or to the
  `analysis/` directory, depending on the notebook. Run notebooks from their
  expected working directory.
- If iSNV loading is slow, keep `analysis/raw_isnv_cache.rds`; it avoids
  repeating the expensive raw iSNV parsing and strand-bias calculations unless
  the source TSV files change.
- Some older exploratory analyses remain in `analysis/`; the current iSNV/COI
  workflow should be taken from `iSNV_COI_V3.Rmd`.
