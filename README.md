# Cas12a Editing Efficiency Processing

This repository contains the minimal processing code needed to convert raw flow-cytometry files, amplicon-sequencing variant tables, and genome-wide long variant tables into editing-efficiency/status tables.

## Workflows

- `scripts/run_fcs_efficiency.R`
  Converts raw `.fcs` files into GFP-positive fractions and GFP editing efficiencies.
- `scripts/amplicon_preprocessing/run_endogenous_amplicon_preprocessing.sh`
  Runs endogenous amplicon FASTQ preprocessing: fastp merge, BBMap mapping,
  optional BAM subsampling, FreeBayes variant calling, and VCF parsing into a
  long variant table.
- `scripts/run_endogenous_amplicon_vcf_to_efficiency.R`
  Converts a concatenated, long FreeBayes variant table for the endogenous amplicon panel into annotated variants, editing-window frequencies, and a Supplementary Table S6-compatible efficiency table.
- `scripts/run_amplicon_vcf_to_efficiency.R`
  Compatibility alias for the endogenous amplicon workflow above.
- `scripts/run_amplicon_variant_efficiency.R`
  Compatibility workflow for tables where variant frequencies and HDR labels have already been calculated.
- `scripts/run_donor_guide_variant_efficiency.R`
  Converts donor-position and shorter-guide long amplicon variant tables into designed-HDR efficiencies.
- `scripts/run_genome_wide_colony_status.R`
  Converts genome-wide long variant tables into colony-level and design-level editing status tables.
- `scripts/export_figure_tables.R`
  Collects only the rows used in manuscript Figures 1-4 into one source-data table per figure.

Detailed commands are in `README_raw_to_efficiency.md`.

After rendering the figure notebooks, export publication source-data tables with:

```sh
Rscript scripts/export_figure_tables.R
```

This writes `results/publication_figure_tables/figure_1_source_data.csv` through `figure_4_source_data.csv`, plus `source_data_summary.csv`.

For the manuscript Figure 1 source-data table, set `CAS12A_PAPER_NOTEBOOKS_ROOT` to the finalized paper-notebook folder before running the exporter. When that folder contains the 2026-05 processed files, Figure 1 uses the May GFP ON cognate/bridge-normalized rows; otherwise the exporter falls back to the rendered `results/04_cas_variants_overview/` tables.

## Requirements

Core table processing uses `dplyr`, `readr`, `stringr`, and `tidyr`.

Raw FCS processing additionally requires Bioconductor flow-cytometry packages:

```r
install.packages("BiocManager")
BiocManager::install(c("flowCore", "flowWorkspace", "openCyto"))
```

## Path Configuration

Raw data are not stored in this repository. Machine-specific input paths are centralized in `config/paths.R`; set the corresponding `CAS12A_*` environment variables or pass script-level path overrides before running the workflows on another machine.

The event-level FCS files used in the paper are distributed as versioned GitHub
Release assets. See [`data/publication_fcs/README.md`](data/publication_fcs/README.md)
for archive contents, checksums, extraction, and analysis commands.
