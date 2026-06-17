# Cas12a Editing Efficiency Processing

This repository contains the minimal processing code needed to convert raw flow-cytometry files, amplicon-sequencing variant tables, and genome-wide long variant tables into editing-efficiency/status tables for the Cas12a manuscript.

## Workflows

- `scripts/run_fcs_efficiency.R`
  Converts raw `.fcs` files into GFP-positive fractions and GFP editing efficiencies.
- `scripts/run_amplicon_variant_efficiency.R`
  Converts long amplicon variant tables into HDR, reference, and non-HDR editing efficiencies for the amplicon panel.
- `scripts/run_donor_guide_variant_efficiency.R`
  Converts donor-position and shorter-guide long amplicon variant tables into designed-HDR efficiencies.
- `scripts/run_genome_wide_colony_status.R`
  Converts genome-wide long variant tables into colony-level and design-level editing status tables.

Detailed commands are in `README_raw_to_efficiency.md`.

## Requirements

Core table processing uses `dplyr`, `readr`, `readxl`, `stringr`, and `tidyr`.

Raw FCS processing additionally requires Bioconductor flow-cytometry packages:

```r
install.packages("BiocManager")
BiocManager::install(c("flowCore", "flowWorkspace", "openCyto"))
```

## Path Configuration

Raw data are not stored in this repository. Machine-specific input paths are centralized in `config/paths.R`; set the corresponding `CAS12A_*` environment variables or pass script-level path overrides before running the workflows on another machine.
