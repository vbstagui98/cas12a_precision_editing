# Cas12a Figure Notebooks

This workspace contains cleaned, shareable replacements for the following original analyses:

- `gfp_lib_20251118_transfer_counts_standardized.Rmd`
- `enas_ampli_colony_count_20260204.R`
- `donor_location_20251121.Rmd`
- `cas_vars_viability_20260208.R`
- `Cas_vars_GFP_ON_20251030.qmd`

The replacements are built as Quarto notebooks and all write their outputs under `results/` in this workspace instead of back into the original project folders.

There are now two notebook styles:

- Full workflow notebooks in `notebooks/`
  These validate inputs, generate cleaned tables, and make plots.
- Plot-only notebooks in `notebooks/plots/`
  These only read the exported tables and make figures. These are the easier files to share with collaborators.

## What changed

- Machine-specific paths are centralized in `config/paths.R`.
- Reusable data loading, validation, plotting, and export helpers live in `R/`.
- Each notebook states its inputs, assumptions, and outputs.
- Plot exports are written to notebook-specific folders under `results/`.
- A render script is included so a colleague can rebuild everything from one command.

## Setup

1. Edit `config/paths.R` if the input data live somewhere else on another machine.
2. Confirm the required R packages are installed.
3. Render either one notebook or all notebooks.

## Install missing packages

```bash
Rscript scripts/install_packages.R
```

## Render one notebook

```bash
quarto render notebooks/01_gfp_library_transfer_counts.qmd
```

## Render all notebooks

```bash
Rscript scripts/render_all.R
```

## Render plot-only notebooks

Run this after the tables have been generated once by the full workflow notebooks.

```bash
Rscript scripts/render_plot_notebooks.R
```

## Notebook map

- `notebooks/01_gfp_library_transfer_counts.qmd`
  Standardizes transfer-count GFP library datasets, produces heatmaps, kinetics plots, and viability summaries.
- `notebooks/02_enas_amplicon_colony_counts.qmd`
  Merges plate-count annotations with amplicon HDR data and exports the correlation and composition plots.
- `notebooks/03_donor_location.qmd`
  Rebuilds donor-position and guide-length figures plus the related colony-count overlays.
- `notebooks/04_cas_variants_overview.qmd`
  Rebuilds the Cas variant GFP ON/OFF processing workflow and exports the main time-course panels.
- `notebooks/05_cas_variants_viability.qmd`
  Joins Cas-variant editing measurements to colony counts and exports viability-versus-efficiency plots.
- `notebooks/06_genome_wide_design_rules.qmd`
  Rebuilds the genome-wide picked-colony editing screen, integrates guide and locus annotations, and exports candidate figure panels plus a cassette-design rule summary.

## Plot-only notebook map

- `notebooks/plots/01_plot_gfp_library_transfer_counts.qmd`
- `notebooks/plots/02_plot_enas_amplicon_colony_counts.qmd`
- `notebooks/plots/03_plot_donor_location.qmd`
- `notebooks/plots/04_plot_cas_variants_overview.qmd`
- `notebooks/plots/05_plot_cas_variants_viability.qmd`

## Notes for collaborators

- The notebooks are designed to fail early with clear messages when an input file or required column is missing.
- Output file names are deterministic so figure links in a manuscript can be updated once and kept stable.
- The helper functions favor cleaned snake_case column names even when the raw files do not.

## Yeast Variant Annotation Workflow

Two scripts were added for yeast variant annotation from a table with `CHR,pos,var,ref,alt`.

- `scripts/build_scepd_promoters.R`
  Converts scEPDnew TSS entries into a promoter BED file using an explicit upstream/downstream window.
- `scripts/annotate_yeast_variants.R`
  Annotates each variant for SpCas9 NGG targetability, CDS/promoter/intergenic class, and SnpEff consequences.

### Inputs

- Variant table: headered `CHR,pos,var,ref,alt` or headerless first five columns `chr pos var ref alt`
- Reference genome FASTA for the same yeast assembly
- GFF3 annotation for the same assembly
- scEPDnew promoter/TSS export for *S. cerevisiae*
- SnpEff database for the same assembly

All references must use the same chromosome naming and assembly coordinates.

### 1. Build a promoter BED from scEPDnew

Download the *Saccharomyces cerevisiae* scEPDnew entries, ideally one representative promoter per gene, then convert them to BED.

```bash
Rscript scripts/build_scepd_promoters.R \
  --scepd path/to/scEPD_scer.tsv \
  --output results/reference/scepd_promoters_500up_100down.bed \
  --upstream 500 \
  --downstream 100
```

The promoter interval is defined relative to the scEPDnew TSS:

- `+` strand: `[TSS - upstream, TSS + downstream]`
- `-` strand: `[TSS - downstream, TSS + upstream]`

### 2. Prepare SnpEff

If `snpEff` is not already on your `PATH`, download the jar and list databases:

```bash
java -jar /path/to/snpEff.jar databases | grep -Ei 'sac|cer|R64'
```

Then download the matching yeast database:

```bash
java -jar /path/to/snpEff.jar download <YEAST_DB_NAME>
```

The exact database name can differ by installation, so inspect the `databases` output first.

### 3. Run the annotation workflow

```bash
Rscript scripts/annotate_yeast_variants.R \
  --variants "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Selected_812/vars_20241106.txt" \
  --genome-fasta path/to/sacCer3.fa \
  --gff path/to/sacCer3.gff3 \
  --promoter-bed results/reference/scepd_promoters_500up_100down.bed \
  --output results/variants/vars_20241106.annotated.tsv \
  --snpeff-jar /path/to/snpEff.jar \
  --snpeff-db <YEAST_DB_NAME>
```

If you only want the Cas9 and region-class annotations first, omit SnpEff with:

```bash
Rscript scripts/annotate_yeast_variants.R \
  --variants "/Users/u0174312/Verstrepen.lab Dropbox/Vlad Batagui/Mac (2)/Documents/PhD/Cas12a_genome_wide/Selected_812/vars_20241106.txt" \
  --genome-fasta path/to/sacCer3.fa \
  --gff path/to/sacCer3.gff3 \
  --promoter-bed results/reference/scepd_promoters_500up_100down.bed \
  --output results/variants/vars_20241106.annotated.tsv \
  --skip-snpeff
```

### Output columns

The main output keeps the original input columns and adds:

- `spcas9_targetable`
- `spcas9_ngg_pam_count`
- `spcas9_ngg_pam_count_plus`
- `spcas9_ngg_pam_count_minus`
- `overlaps_cds`
- `cds_count`
- `cds_ids`
- `overlaps_promoter`
- `promoter_count`
- `promoter_ids`
- `region_class`
- `snpeff_ann_count`
- `snpeff_primary_effect`
- `snpeff_primary_impact`
- `snpeff_primary_gene`
- `snpeff_primary_feature_id`
- `snpeff_primary_hgvs_c`
- `snpeff_primary_hgvs_p`
- `snpeff_all_effects`
- `snpeff_all_impacts`
- `snpeff_all_genes`

`region_class` is mutually exclusive with precedence `coding > promoter > intergenic`.
