# Cas12a precision editing

This repository contains the scripts used to calculate editing efficiencies
for the Cas12a precision-editing manuscript.

## Analysis workflows

- `scripts/run_fcs_efficiency.R`: raw FCS files to GFP ON/OFF efficiency.
- `scripts/amplicon_preprocessing/run_endogenous_amplicon_preprocessing.sh`:
  paired FASTQ files to one long FreeBayes variant table.
- `scripts/run_endogenous_amplicon_vcf_to_efficiency.R`: endogenous-panel
  variant table to HDR, reference, and non-HDR frequencies.
- `scripts/run_donor_guide_variant_efficiency.R`: donor-position and
  shorter-guide variant tables to HDR efficiency.
- `scripts/run_genome_wide_colony_status.R`: genome-wide variant tables to
  colony and design editing status.

Commands and input-column requirements are described in
`README_raw_to_efficiency.md`. A longer description of the functions is in
`docs/repository_walkthrough.md`.

## Publication data

Small inputs used to test the workflows are in `data/publication_inputs/`.
The event-level FCS files are available from the
[`fcs-data-v1` release](https://github.com/vsbatagui/cas12a_precision_editing/releases/tag/fcs-data-v1).

## R packages

The table-processing scripts use `dplyr`, `readr`, `stringr`, and `tidyr`.
Raw FCS processing also requires:

```r
install.packages("BiocManager")
BiocManager::install(c("flowCore", "flowWorkspace", "openCyto"))
```
