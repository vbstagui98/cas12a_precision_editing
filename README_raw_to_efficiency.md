# Raw Inputs To Editing Efficiencies

This document covers the processing paths needed for the publication upload:

- raw `.fcs` files to GFP editing efficiencies;
- endogenous-panel long amplicon-sequencing variant tables to HDR/ref/non-HDR efficiencies;
- donor-position and shorter-guide amplicon variant tables to designed-HDR efficiencies;
- genome-wide long variant tables to colony-level editing status.

All command-line table inputs are CSV, TSV, or TXT files.

## FCS To GFP Editing Efficiency

Download and extract the relevant archive from the
[`fcs-data-v1` release](https://github.com/vsbatagui/cas12a_precision_editing/releases/tag/fcs-data-v1).
The committed release manifest and archive checksums are documented in
`data/publication_fcs/README.md`.

```bash
Rscript scripts/run_fcs_efficiency.R \
  --fcs-dir path/to/FACS_DATA \
  --platemap path/to/platemap.csv \
  --output results/fcs_efficiency \
  --guide-features path/to/guides_donor_selected_info.csv
```

The FCS workflow follows the original Transfer_counts notebooks:

- read raw `.fcs` files;
- filter wells with `<= 500` events;
- estimate one logicle transform from retained samples;
- gate GFP-positive cells on `BL1.A`;
- report `gfp_positive_fraction` and `editing_efficiency`.

Guide annotation CSVs may use `crRNA_id`, `guide_id`, `Guide`, or `guide`.
The output retains the original analysis fields and also provides the
publication names `efficiency_pct`, `colony_count`, `control_colony_count`,
`relative_viability`, and `viability_normalization_sample`.

Editing direction:

- `GFP_OFF`: `editing_efficiency = 1 - gfp_positive_fraction`
- `GFP_ON`: `editing_efficiency = gfp_positive_fraction`

## Endogenous Concatenated FreeBayes Variants To Amplicon Efficiency

The primary publication input is the long table produced after parsing and
concatenating the FreeBayes VCF files from one sequencing run. It should follow
the original `scripts_amplicons` preprocessing chain: merge/trim reads with
`fastp`, map with BBMap, call variants with FreeBayes using
`--min-base-quality 3 --haplotype-length 30 -C 1 --min-coverage 50
--min-alternate-fraction 0.001 --pooled-continuous`, parse/split the VCF
records, calculate `frc_alt = AO / DP * 100` and
`frc_ref = RO / DP * 100`, filter `DP > 10000`, and derive `mismatches`,
`positions`, and `pos_mismatch`.

A runnable, path-independent implementation of this preprocessing chain is in
`scripts/amplicon_preprocessing/`. For example:

```bash
bash scripts/amplicon_preprocessing/run_endogenous_amplicon_preprocessing.sh \
  --fastq-dir path/to/fastq_merged_lanes \
  --analysis-dir results/amplicon_seq \
  --reference path/to/saccharomyces_cerevisiae_sequence.fasta \
  --sample-glob '*VB_Fn*_R1.fastq.gz'
```

```bash
Rscript scripts/run_endogenous_amplicon_vcf_to_efficiency.R \
  --variants path/to/concatenated_freebayes_variants.csv \
  --designs path/to/intended_loci.csv \
  --guide-features path/to/endogenous_guide_features.csv \
  --output-dir results/amplicon_efficiency
```

### Variant input

Each row represents one alternate allele from one VCF sample. Multi-allelic
VCF records should be split so that each alternate allele has its own row.

| Meaning | Accepted examples |
|---|---|
| Sample ID | `Sample`, `sample` |
| Guide ID | `Guide`, `crRNA_id`, `guide_id` |
| Chromosome | `CHROM`, `chromosome`, `chrom`, `chr` |
| VCF position | `POS`, `position`, `pos` |
| Reference allele | `REF`, `ref` |
| Alternate allele | `ALT`, `alt` |
| FreeBayes variant type | `TYPE`, `type`, `variant_type` |
| Alternate read count | `AO`, `ao`, `alt_depth` |
| Reference read count | `RO`, `ro`, `ref_depth` |
| Read depth | `DP`, `dp`, `depth` |
| Parser allele frequency | `frc_alt`, `allele_frequency_pct`, optional |
| Parser reference frequency | `frc_ref`, `reference_pct`, optional |

Sample metadata such as `Experiment`, `Timepoint`, `Replicate`, `promoter`,
`recruit`, `cas_variant`, and `DR_cognate` are retained. For the historical
`E1`-`E4` design, promoter and donor-recruitment values are derived from
`Experiment` when they are absent.

`QUAL` is optional and is never filtered. `mismatches`, `positions`, and
`pos_mismatch` are calculated from `REF`, `ALT`, and `POS` if the VCF parser
did not already produce them.

### Design input

The design CSV needs `crRNA_id` (or `Guide`) and an exact intended-locus
definition. Either provide:

- the original `match` column, formatted as
  `promoter_Guide_chromosome_position_alternate`; or
- `intended_chromosome`, `intended_position`, and `alternate_allele`.

The original manuscript assignment file also works directly because the
aliases `promoter`, `Guide`, `chr`, `position`, `ALT`, and `match` are
recognized.

`--guide-features` is optional. It accepts a CSV with the Supplementary Table
S5 feature names and copies those guide/donor annotations into the output.
Its `variant_position_annotation` field is not used to assign HDR:
in the current publication table that field records the PAM coordinate, not
the edited nucleotide coordinate, and `distance_from_pam` does not retain the
direction needed to reconstruct every intended locus.

### Calculation

- editing window = all haplotypes that co-ocur with the HDR haplotype. Thus HDR +
  WT allele + other alleles = 100
- `HDR`: a detected row exists after parser/depth filtering where
  `match = paste(promoter, Guide, CHROM, pos_mismatch, mismatches, sep = "_")`
  is present in the design table, `TYPE` is one of `snp`, `mnp`, or `complex`,
  and `nchar(ALT) == nchar(REF)`.
- Other non-reference, non-intended variants are retained
- For samples with detected HDR, the editing window contains all variant rows
  with the same `Sample`, `CHROM`, `POS`, and `REF` as the intended call.
- One reference row is added using the intended call's `RO`.
- Frequencies are calculated as
  `count / (reference RO + sum of alternate AO in the editing window) * 100`.
- Samples without a matched intended row receive an explicit HDR 0% placeholder
  and reference 100%.

The script defaults to `--min-dp 10000`, matching the original parser's
`DP > 10000` filter. Use `--min-dp 0` only when rerunning already curated final
manuscript count tables that should not receive an additional depth filter.

Outputs:

- `amplicon_variants_annotated.csv`: all input variants plus design
  annotations, count-derived raw frequencies, and assignment columns.
- `amplicon_editing_window.csv`: HDR, reference, and other alleles
  at the intended locus with frequencies recomputed from counts.
- `amplicon_efficiency.csv`: a compact table using the publication
  names from Supplementary Table S6 where applicable.

The result was regression-tested against the historical
`var_file_20250328.csv` workflow. The recalculated FnCas12a HDR and reference
percentages match the old processed output to floating-point precision.

### Compatibility Input

`scripts/run_amplicon_variant_efficiency.R` remains available for previously
processed tables where `MUTATION`/`editing_class`, `frc_alt`/`efficiency_pct`,
and `frc_ref`/`reference_pct` already exist. It is not the primary raw-data
workflow.

`scripts/run_amplicon_vcf_to_efficiency.R` remains as a compatibility alias for
`scripts/run_endogenous_amplicon_vcf_to_efficiency.R`.

## Donor-Position And Shorter-Guide Amplicon Variant Tables

This is intentionally separate from the endogenous-panel workflow. It follows
the Figure 3 donor-position and shorter-guide analyses, where the designed HDR
allele is assigned from the expected edited SNP for that assay-specific design.
It does not use the endogenous guide/promoter intended-locus table.

```bash
Rscript scripts/run_donor_guide_variant_efficiency.R \
  --output-dir results/donor_guide_variant_efficiency
```

Outputs:

- `donor_position_hdr_from_variants.csv`
- `short_guide_hdr_from_variants.csv`

Inputs are the long variant tables generated for those experiments
(`master_df_filtered.csv` in the original analysis folders). Outputs are one
designed-HDR efficiency table for donor position and one for guide length.

## Genome-Wide Variant Tables To Colony Status

```bash
Rscript scripts/run_genome_wide_colony_status.R \
  --output-dir results/genome_wide_colony_status \
  --picked-colonies-round-1 path/to/picked_colonies_round_1.csv \
  --picked-colonies-round-2 path/to/picked_colonies_round_2.csv \
  --target-coverage-round-1 path/to/target_coverage_round_1.tsv \
  --target-coverage-round-2 path/to/target_coverage_round_2.tsv \
  --mean-coverage-round-1 path/to/mean_coverage_round_1.csv \
  --mean-coverage-round-2 path/to/mean_coverage_round_2.csv \
  --normalized-variants-round-1 path/to/variants_norm_round_1.csv \
  --normalized-variants-round-2 path/to/variants_norm_round_2.csv \
  --design-annotations path/to/designs_812_sv_my_variants_annotated.tsv \
  --deepcpf1-scores path/to/input_deepcpf1_20260319.scored.csv \
  --dr-scores path/to/guides_dr_scores_20260325.tsv \
  --normalized-indels path/to/normalized_indel_designs.vcf
```

The genome-wide workflow uses the manuscript rule set:

- intended editing: intended HDR allele has `AO >= 2`;
- non-target variants: `frc_alt >= 50` and `DP >= 4`;
- no `QUAL` threshold.

Outputs:

- `genome_wide_colony_editing_status.csv`
- `genome_wide_design_editing_status.csv`
- `genome_wide_non_target_af50_dp4_variants.csv`

Optional path overrides are listed in:

```bash
Rscript scripts/run_genome_wide_colony_status.R --help
```

## Schema Checks

The command-line scripts use delimited text, not the supplementary workbook.
The harmonization layer was checked against the workbook column names:

- Supplementary Tables S1/S3: FCS outputs include the publication efficiency
  and viability names.
- Supplementary Table S4/S6: raw amplicon outputs use `editing_class`,
  `efficiency_pct`, `reference_pct`, `crRNA_id`, `crRNA_promoter`,
  `direct_repeat`, `donor_recruitment`, and `generations`.
- Supplementary Table S5: guide and donor feature names are accepted. Exact
  HDR assignment additionally requires `match` or explicit intended-locus
  columns because the published position field records the PAM coordinate.
- Supplementary Tables S7/S8: the genome-wide workflow retains its dedicated
  colony/design input files and applies the genome-wide manuscript rules:
  intended `AO >= 2`, non-target `AF >= 50%` and `DP >= 4`, no `QUAL` filter.

Run the table-level regression tests with:

```bash
Rscript tests/test_amplicon_variant_efficiency.R
Rscript tests/test_fcs_table_harmonization.R
```
