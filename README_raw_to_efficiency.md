# Raw Inputs To Editing Efficiencies

Inputs and outputs 

- raw `.fcs` files to GFP editing efficiencies;
- long amplicon-sequencing variant tables to HDR/ref/non-HDR efficiencies;
- genome-wide long variant tables to colony-level editing status.

## FCS To GFP Editing Efficiency

```bash
Rscript scripts/run_fcs_efficiency.R \
  --fcs-dir path/to/FACS_DATA \
  --platemap path/to/platemap.csv \
  --output results/fcs_efficiency \
  --guide-features path/to/guides_donor_selected_info.csv
```

The FCS workflow follows the original analysis workflow:

- read raw `.fcs` files;
- filter wells with `<= 500` events;
- estimate one logicle transform from retained samples;
- gate GFP-positive cells on `BL1.A`;
- report `gfp_positive_fraction` and `editing_efficiency`.

Editing efficiency calculated as a function of the GFP assay (total cleavage = GFP_OFF or HDR = GFP_ON)
- `GFP_OFF`: `editing_efficiency = 1 - gfp_positive_fraction`
- `GFP_ON`: `editing_efficiency = gfp_positive_fraction`

## Amplicon Variant Tables To Editing Efficiency

For the amplicon panel, the following workflow uses pre-comptued allele frequencies derived from vcf files:

```bash
Rscript scripts/run_amplicon_variant_efficiency.R \
  --variants path/to/genomic_ampli_enas_fn_harmonized_20260318.csv \
  --output results/amplicon_panel_efficiency.csv \
  --guide-features results/panel_sequence_manifest/combined_sequence_manifest.csv
```

The input must be a long table where intended edits have already been labelled in `MUTATION` and allele frequencies are present in `frc_alt` and `frc_ref`.

For donor-position and shorter-guide amplicon sequencing:

```bash
Rscript scripts/run_donor_guide_variant_efficiency.R \
  --output-dir results/donor_guide_variant_efficiency
```

Outputs:

- `donor_position_hdr_from_variants.csv`
- `short_guide_hdr_from_variants.csv`

The donor-position and shorter-guide helper assigns designed HDR from the raw long variant tables by matching the expected designed SNP.

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

Outputs:

- `genome_wide_colony_editing_status.csv`
- `genome_wide_design_editing_status.csv`
- `genome_wide_non_target_af50_dp4_variants.csv`

Optional path overrides are listed in:

```bash
Rscript scripts/run_genome_wide_colony_status.R --help
```
