# Publication analysis inputs

These files are path-independent inputs for the public analysis scripts. They
contain only samples represented in the final manuscript source tables.

## Endogenous amplicon sequencing

Run:

```bash
Rscript scripts/run_endogenous_amplicon_vcf_to_efficiency.R \
  --variants data/publication_inputs/endogenous_freebayes_variants_publication.csv \
  --designs data/publication_inputs/endogenous_intended_loci.csv \
  --guide-features data/publication_inputs/endogenous_guide_features.csv \
  --output-dir results/endogenous_amplicon \
  --min-dp 0
```

Inputs:

- `endogenous_freebayes_variants_publication.csv`: one count-level variant
  record per detected variant or explicit zero-call placeholder. It retains
  `AO`, `RO`, `DP`, and the pre-filter allele-count denominator
  `total_count`; `frc_alt` and `frc_ref` are deliberately omitted.
- `endogenous_intended_loci.csv`: exact HDR assignments by crRNA and crRNA
  promoter. Alternative equivalent FreeBayes representations are separated by
  `|` in `intended_match`.
- `endogenous_guide_features.csv`: crRNA, donor, PAM, 5' PAM nucleotide,
  distance from PAM, strand, GC content, and DeepCpf1 score.
- `endogenous_sample_metadata.csv`: one row per Figure 4 sample.

The generic raw-parser workflow defaults to the original `DP > 10000` filter
from `scripts_amplicons/04_parse_vcf.ipynb`. This committed CSV is a curated
final manuscript input that already reflects the plotted data and retains a
correction row, so the reproducibility command above uses `--min-dp 0`.

The workflow uses no QUAL filter. Intended amplicon HDR is assigned when a
matched intended substitution row exists after the parser/depth-filtering
stage; no separate `AO >= 2` threshold is applied to amplicon data.
Non-intended variants are retained as ordinary amplicon variants; the AF50/DP4
unwanted-edit rule is used only by the genome-wide colony-status workflow.

## Donor position and shorter guides (Figure 3)

The original analysis folder calls this workflow `amplicon_seq`. Run the
committed data directly with:

```bash
Rscript scripts/run_donor_guide_variant_efficiency.R \
  --fn-variants data/publication_inputs/ampliconseq_fn_donor_variants.csv \
  --enas-variants data/publication_inputs/ampliconseq_enas_donor_short_guide_variants.csv \
  --output-dir results/donor_guide_variant_efficiency
```

- `ampliconseq_fn_donor_variants.csv` contains the FnCas12a donor-position
  samples and the corrected repeat-sample assignments used by the notebook.
- `ampliconseq_enas_donor_short_guide_variants.csv` is the mixed sequencing-run
  table containing the enAsCas12a donor-position samples and shorter-guide
  samples for both Cas12a variants.

The script matches the assay-specific intended SNP and reports the detected
row's `frc_alt` as HDR efficiency. It does not apply the endogenous-panel
editing-window calculation or the genome-wide AF50/DP4 rule.
The output `figure_3_ampliconseq_efficiency.csv` applies the manuscript
timepoint selection and contains the 124 Figure 3 rows.

## Flow cytometry

The raw workflow requires the original `.fcs` files:

```bash
Rscript scripts/run_fcs_efficiency.R \
  --fcs-dir path/to/FACS_DATA \
  --platemap data/publication_inputs/fcs_publication_platemap.csv \
  --output results/fcs_efficiency
```

`fcs_publication_platemap.csv` has one row per unique FCS file used in the
Figure 1 or Figure 2 analyses. Filter it to the relevant `experiment_key`
before processing an individual FCS directory.

For reviewers who do not download the event-level FCS archive,
`fcs_population_frequencies_publication.csv` contains the GFP-gate output.
The annotation and GFP-to-editing conversion can be checked with:

```r
source("R/fcs_efficiency.R")

pop_stats <- readr::read_csv(
  "data/publication_inputs/fcs_population_frequencies_publication.csv",
  show_col_types = FALSE
)

efficiency <- annotate_fcs_efficiency(
  pop_stats,
  "data/publication_inputs/fcs_publication_platemap.csv"
)
```

Compare the result with `fcs_expected_efficiencies_publication.csv`.

## Raw data availability

CSV files cannot replace event-level FCS measurements. The `.fcs` files are
distributed as GitHub Release assets; see
[`../publication_fcs/README.md`](../publication_fcs/README.md). FASTQ and BAM
files should remain linked through ENA/SRA rather than committed to GitHub.

`manifest.csv` records dimensions and descriptions for the CSV inputs.
