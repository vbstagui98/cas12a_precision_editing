# Repository Walkthrough For Publication Reproducibility

This document explains the public repository structure and the assumptions used
by the reproducibility scripts. The source of truth for amplicon-sequencing
assumptions is the original `scripts_amplicons` workflow, especially
`04_parse_vcf.ipynb` and `var_freq_R_final_20250328.Rmd`.

## Top-Level Structure

- `README.md`
  Short repository overview and dependency notes.

- `README_raw_to_efficiency.md`
  Main runbook for converting raw or near-raw data into efficiency/status
  tables.

- `config/paths.R`
  Local path defaults and `CAS12A_*` environment-variable overrides.

- `R/`
  Reusable functions. Public command-line scripts should source functions from
  here rather than implementing analysis logic inline.

- `scripts/`
  Command-line entry points and helper/export scripts.

- `data/publication_inputs/`
  Small, path-independent CSV inputs that reviewers can run without access to
  local Dropbox paths.

- `tests/`
  Regression tests for amplicon and FCS table behavior.

## Public Workflows

### FCS Workflow

Entry point:

```bash
Rscript scripts/run_fcs_efficiency.R \
  --fcs-dir path/to/FACS_DATA \
  --platemap path/to/platemap.csv \
  --output results/fcs_efficiency
```

Purpose:

- read raw `.fcs` files;
- remove samples with too few events;
- apply a common logicle transform;
- gate GFP-positive cells on `BL1.A`;
- convert GFP-positive fraction into editing efficiency;
- optionally attach guide annotations and colony-count viability.

Reviewer shortcut:

- `data/publication_inputs/fcs_population_frequencies_publication.csv`
  contains precomputed GFP population frequencies, so reviewers can test
  metadata joins and GFP-to-efficiency conversion without raw FCS files.

### Endogenous Amplicon Workflow

Entry point:

```bash
Rscript scripts/run_endogenous_amplicon_vcf_to_efficiency.R \
  --variants data/publication_inputs/endogenous_freebayes_variants_publication.csv \
  --designs data/publication_inputs/endogenous_intended_loci.csv \
  --guide-features data/publication_inputs/endogenous_guide_features.csv \
  --output-dir results/endogenous_amplicon \
  --min-dp 0
```

Generic raw-parser input assumptions:

1. merge/trim reads with `fastp`;
2. map with BBMap;
3. call variants with FreeBayes:

   ```bash
   freebayes \
     --min-base-quality 3 \
     --haplotype-length 30 \
     -C 1 \
     --min-coverage 50 \
     --min-alternate-fraction 0.001 \
     --pooled-continuous
   ```

4. parse split variants from VCF;
5. calculate `frc_alt = AO / DP * 100`;
6. calculate `frc_ref = RO / DP * 100`;
7. filter `DP > 10000`;
8. derive `mismatches`, `positions`, and `pos_mismatch`;
9. assign HDR in R using the guide-design `match` column.

Important distinction:

- The script default is `--min-dp 10000`, matching the original parser.
- The committed publication CSV is already curated from final manuscript data
  and includes a retained correction row, so the reviewer command uses
  `--min-dp 0`.

Amplicon HDR rule:

- Build `match = paste(promoter, Guide, CHROM, pos_mismatch, mismatches, sep = "_")`.
- A row is HDR if:
  - `TYPE` is `snp`, `mnp`, or `complex`;
  - `nchar(ALT) == nchar(REF)`;
  - `match` or `intended_variant_id` matches the design table.
- There is no amplicon `AO >= 2` threshold.
- A matched intended row with `AO = 1` is HDR if it survived the parser/depth
  filtering.
- Non-intended rows keep their FreeBayes `TYPE` in `MUTATION` and are marked
  `other_variant` in `variant_annotation`.
- The genome-wide `AF >= 50%` and `DP >= 4` non-target rule is not used here.

Editing-window rule:

- For samples with detected HDR, keep every variant row with the same
  `Sample`, `CHROM`, `POS`, and `REF` as the matched HDR row.
- Add one synthetic `REF` row with `AO = first(RO)`.
- Recompute frequencies as `AO / sum(AO) * 100` within each sample.
- For samples with no matched HDR row, add an explicit synthetic HDR row at 0%
  and a synthetic REF row at 100%.

### Donor-Position And Shorter-Guide Amplicon Workflow

Entry point:

```bash
Rscript scripts/run_donor_guide_variant_efficiency.R \
  --output-dir results/donor_guide_variant_efficiency
```

This workflow is separate from the endogenous panel. It uses assay-specific
expected SNP positions and the historical `Part*` metadata columns from those
amplicon experiments.

### Genome-Wide Colony-Status Workflow

Entry point:

```bash
Rscript scripts/run_genome_wide_colony_status.R --help
```

Genome-wide rules:

- intended edit: intended HDR allele has `AO >= 2`;
- non-target variants: `frc_alt >= 50` and `DP >= 4`;
- no `QUAL` threshold.

These rules must stay in the genome-wide workflow and must not be imported into
amplicon-sequencing workflows.

## Main Functions

### `R/fcs_efficiency.R`

- `check_flow_packages()`
  Confirms that `flowCore`, `flowWorkspace`, and `openCyto` are installed for
  raw event-level FCS processing.

- `read_csv_keep_names(path)`
  Reads comma- or semicolon-delimited metadata while preserving original column
  names.

- `repair_blank_index_columns(data)`
  Names blank columns as `source_index_*`.

- `normalise_well_name(x)`
  Converts wells such as `A1` to `A01`.

- `read_platemap_for_fcs(path)`
  Reads and validates the FCS plate map. Requires a `name` column.

- `read_optional_table(path)`
  Reads optional metadata/guide/count tables or returns `NULL`.

- `copy_alias_if_missing(data, canonical, aliases)`
  Copies the first available alias into a canonical column.

- `harmonise_publication_join_keys(data)`
  Adds common join aliases such as `guide_id` and `WellPosition`.

- `make_common_logicle_transform(fs, channels)`
  Estimates one common transform for retained FCS samples.

- `gate_gfp_population(fstrans, gfp_channel, gate_range)`
  Gates GFP-positive cells and returns population frequencies.

- `read_gate_fcs(...)`
  Raw FCS pipeline: read files, filter event counts, transform, gate GFP, and
  return population stats.

- `standardise_pop_stats(pop_stats)`
  Converts gate output to `name`, `population`, `parent`,
  `gfp_positive_fraction`, and `parent_frequency`.

- `join_optional_by_first_key(data, extra, keys)`
  Joins optional metadata using the first shared key.

- `editing_efficiency_from_gfp(gfp_fraction, assay)`
  Converts GFP fraction to editing efficiency: GFP ON uses the fraction, GFP
  OFF uses `1 - fraction`.

- `normalise_control_count(colony_table, control_well = "A09")`
  Calculates relative viability from `true_colony_count / control_colony_count`
  within matched timepoint/replicate/dataset groups.

- `annotate_fcs_efficiency(...)`
  Joins metadata, converts GFP fraction to editing efficiency, and adds
  publication columns.

- `run_fcs_efficiency(...)`
  Top-level writer for `pop_stats.csv`, `event_counts.csv`, and
  `fcs_efficiency.csv`.

### `R/amplicon_variant_efficiency.R`

- `normalize_column_key(x)`
  Normalizes column names for alias matching.

- `copy_column_alias(data, canonical, aliases)`
  Adds canonical columns from accepted aliases.

- `standardize_amplicon_columns(data)`
  Harmonizes variant table names such as `sample`/`Sample`,
  `crRNA_id`/`Guide`, and `direct_repeat`/`DR_cognate`.

- `read_variant_table(path)`
  Reads CSV/TSV/TXT variant tables.

- `write_csv_mkdir(data, path)`
  Writes CSV after creating the output directory.

- `check_required_files(paths)`
  Stops if required input files are missing.

- `ensure_column()` and `ensure_columns()`
  Add missing optional columns with `NA`.

- `standardize_guide_features(guide_features)`
  Harmonizes guide/design feature columns and preserves intended-locus aliases.

- `amplicon_efficiency_from_long_table(...)`
  Compatibility function for already processed tables with `MUTATION`,
  `frc_alt`, and `frc_ref`.

- `normalize_chromosome(x)`
  Makes chromosome names comparable.

- `contains_integer(text, expected)`
  Checks whether a position string contains the expected integer.

- `matches_intended_alias(candidate, intended_aliases)`
  Tests pipe-separated intended variant aliases.

- `derive_mismatch_fields(ref, alt, pos)`
  Reconstructs `mismatches`, `positions`, and `pos_mismatch` when the parser did
  not provide them.

- `standardize_endogenous_designs(design_table)`
  Builds exact intended-locus columns from `match`/`intended_match` or explicit
  chromosome/position/allele columns.

- `merge_endogenous_design_features(design_assignments, guide_features)`
  Adds guide/donor/PAM/DeepCpf1 metadata to exact intended-locus assignments.

- `annotate_endogenous_amplicon_variants(...)`
  Implements the original amplicon HDR assignment: matching substitution rows
  are HDR; all other non-reference rows remain ordinary variants. No AO
  threshold is applied.

- `build_endogenous_editing_window(annotated_variants)`
  Builds the intended-locus editing window, adds synthetic REF/HDR placeholders,
  and recomputes percentages from counts.

- `endogenous_hdr_efficiency_table(annotated_variants)`
  Builds one sample-level HDR table for regression against manuscript values.

- `endogenous_amplicon_efficiency_table(editing_window)`
  Formats the editing-window rows into publication-style output columns.

- `simplify_variant_table(df)`
  Simplifies complex VCF alleles for donor-position/shorter-guide workflows.

- `complement_nuc(nuc)`
  Returns the nucleotide complement.

- `summarise_donor_position_hdr(...)`
  Matches expected donor-position SNPs and summarizes HDR percentages.

- `summarise_fixed_hdr(...)`
  Matches the fixed shorter-guide expected HDR SNP.

- `prepare_fn_donor_variants(variants)`
  Applies historical Fn donor-position cleanup.

- `prepare_enas_donor_variants(variants)`
  Keeps enAs donor rows from the mixed input table.

- `prepare_short_guide_variants(variants)`
  Extracts guide length, strain, replicate, and timepoint from historical
  columns.

- `prepare_donor_and_guide_length_outputs(paths, output_dir)`
  Top-level donor-position/shorter-guide writer.

### `R/genome_wide_colony_status.R`

- `make_colony_assay_table(...)`
  Builds colony/design assay rows and applies the genome-wide intended
  `AO >= 2` rule.

- `make_non_target_tables(...)`
  Applies the genome-wide non-target `frc_alt >= 50` and `DP >= 4` rule.

- `build_genome_wide_colony_status(...)`
  Orchestrates colony/design/non-target status construction.

- `run_genome_wide_colony_status(input_paths, output_dir)`
  Writes the genome-wide colony, design, and non-target output CSVs.

## Tests To Read First

- `tests/test_amplicon_variant_efficiency.R`
  Verifies that amplicon HDR is detected by matched row existence, not by an AO
  threshold; verifies zero-HDR placeholders, synthetic REF rows, count-based
  frequency recomputation, and absence of `unintended_AF50_DP4`.

- `tests/test_fcs_table_harmonization.R`
  Verifies guide aliases, GFP ON/OFF conversion, and colony-count normalization.

## Validation Commands

```bash
Rscript tests/test_amplicon_variant_efficiency.R
Rscript tests/test_fcs_table_harmonization.R
Rscript scripts/run_endogenous_amplicon_vcf_to_efficiency.R \
  --variants data/publication_inputs/endogenous_freebayes_variants_publication.csv \
  --designs data/publication_inputs/endogenous_intended_loci.csv \
  --guide-features data/publication_inputs/endogenous_guide_features.csv \
  --output-dir /tmp/cas12a_endogenous_amplicon \
  --min-dp 0
```
