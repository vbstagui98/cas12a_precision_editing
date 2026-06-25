# Endogenous Amplicon-Seq Preprocessing

This folder contains a path-independent version of the original
`scripts_amplicon` preprocessing workflow used for the endogenous/Fn-style
amplicon panel.

The public downstream R workflow starts from the long FreeBayes table produced
here:

```bash
Rscript scripts/run_endogenous_amplicon_vcf_to_efficiency.R \
  --variants results/amplicon_seq/var_file.csv \
  --designs data/publication_inputs/endogenous_intended_loci.csv \
  --guide-features data/publication_inputs/endogenous_guide_features.csv \
  --output-dir results/endogenous_amplicon
```

## Master Script

Run all preprocessing stages:

```bash
bash scripts/amplicon_preprocessing/run_endogenous_amplicon_preprocessing.sh \
  --fastq-dir path/to/fastq_merged_lanes \
  --analysis-dir results/amplicon_seq \
  --reference path/to/saccharomyces_cerevisiae_sequence.fasta \
  --sample-glob '*VB_Fn*_R1.fastq.gz'
```

Run only a subset of samples:

```bash
bash scripts/amplicon_preprocessing/run_endogenous_amplicon_preprocessing.sh \
  --fastq-dir path/to/fastq_merged_lanes \
  --analysis-dir results/amplicon_seq \
  --reference path/to/saccharomyces_cerevisiae_sequence.fasta \
  --sample-list sample_names.txt
```

`sample_names.txt` should contain one sample per line. The script accepts names
with or without the trailing Illumina `_S...` index.

Run only selected stages:

```bash
bash scripts/amplicon_preprocessing/run_endogenous_amplicon_preprocessing.sh \
  --fastq-dir path/to/fastq_merged_lanes \
  --analysis-dir results/amplicon_seq \
  --reference path/to/saccharomyces_cerevisiae_sequence.fasta \
  --steps freebayes,parse
```

Valid stages are `fastp`, `map`, `subsample`, `freebayes`, and `parse`.
For `--steps parse`, only `--analysis-dir` is required; the script reads VCFs
from `<analysis-dir>/vcf`.

## Preserved Processing Assumptions

The master script preserves the original command choices:

- `fastp` merges paired reads with `-m`, trims the same Illumina adapter
  sequences, and keeps reads with `--length_required=100`.
- `bbmap.sh` maps merged reads with `slow k=12 minratio=0.9`.
- `samtools` optionally subsamples sorted BAM files to 100,000 reads by
  default, using the same fraction calculation style as the original script.
- `freebayes` is called with:

```bash
freebayes \
  -f reference.fasta \
  sample.subset.bam \
  --min-base-quality 3 \
  --haplotype-length 30 \
  -C 1 \
  --min-coverage 50 \
  --min-alternate-fraction 0.001 \
  --pooled-continuous
```

- The parser keeps `QUAL` but never filters on it.
- The parser calculates `frc_alt = AO / DP * 100` and
  `frc_ref = RO / DP * 100`.
- The parser keeps rows with `DP > 10000` by default.
- The parser derives `mismatches`, `positions`, and `pos_mismatch` exactly for
  the downstream HDR matching step.

## Outputs

The master script writes:

- `merged_fastp/`: merged FASTQ files from `fastp`;
- `html_fastp/`: `fastp` HTML and JSON reports;
- `Bam/`: SAM, BAM, sorted BAM, and BAM index files from mapping;
- `Bam_stats/`: BBMap scaffold statistics;
- `bam_subset/`: BAM files used for FreeBayes;
- `vcf/`: one FreeBayes VCF per sample;
- `var_file.csv`: concatenated long variant table.

`var_file.csv` is the intended input to
`scripts/run_endogenous_amplicon_vcf_to_efficiency.R`.

## Dependencies

Install these tools outside the repository:

- `fastp`
- BBMap, providing `bbmap.sh`
- `samtools`
- `freebayes`
- `python3`

The parser uses only the Python standard library.
