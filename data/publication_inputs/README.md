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
  --output-dir results/endogenous_amplicon
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

The workflow uses no QUAL filter. Intended editing requires at least two
alternate reads. Non-intended variants are retained as ordinary amplicon
variants; the AF50/DP4 unwanted-edit rule is used only by the genome-wide
colony-status workflow.

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

CSV files cannot replace event-level FCS measurements. The `.fcs` files should
be deposited in Zenodo, FlowRepository, or another archival repository and
linked from the main project README. FASTQ and BAM files should likewise be
linked through ENA/SRA rather than committed to GitHub.

`manifest.csv` records dimensions, descriptions, and SHA256 checksums for all
CSV inputs.
