# Publication source data

This folder contains the figure-source workbook and clean long variant tables committed for the Cas12a precision-editing publication repository.

## Files

- `publication_figure_source_tables_clean_20260617.xlsx`: one worksheet per final figure source-data table.
- `raw_variant_tables/donor_position_all_detected_variants.csv`: all detected variants from the donor-position amplicon sequencing runs after applying the same donor sample annotations used for the efficiency summaries.
- `raw_variant_tables/short_guide_all_detected_variants.csv`: all detected variants from the short-guide amplicon sequencing run after applying the same guide-length annotations used for the efficiency summaries.
- `raw_variant_tables/genomic_panel_all_detected_variants.csv`: all detected variants from `genomic_ampli_enas_fn_harmonized_20260318.csv`, the long table used for the genomic amplicon panel figure.
- `manifest.csv`: row counts and source paths used for the export.

Regenerate these files from the project root with:

```bash
CAS12A_PAPER_NOTEBOOKS_ROOT="/path/to/cas12a_paper_notebooks_20260312" \
CAS12A_PRESENTATIONS_AMPLICON_ROOT="/path/to/9-month_figs/amplicon_seq" \
Rscript scripts/export_publication_source_data.R
```
