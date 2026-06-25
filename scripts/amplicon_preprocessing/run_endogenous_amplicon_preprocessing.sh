#!/usr/bin/env bash

set -euo pipefail

## required packages:
## fastp, bbmap/bbmap.sh, samtools, freebayes, python3
##
## This script takes in a directory containing paired fastq reads and runs the
## endogenous amplicon-seq preprocessing workflow used before the public R
## efficiency scripts:
##   1. merge reads and trim adapters with fastp
##   2. map merged reads to a reference genome with bbmap
##   3. optionally subsample sorted BAM files
##   4. call variants with freebayes
##   5. parse the VCF files into one long variant table

usage() {
    cat <<'EOF'
Usage:
  bash scripts/amplicon_preprocessing/run_endogenous_amplicon_preprocessing.sh \
    --fastq-dir path/to/fastq_merged_lanes \
    --analysis-dir path/to/amplicon_seq_analysis \
    --reference path/to/saccharomyces_cerevisiae_sequence.fasta

Required:
  --analysis-dir    Directory where merged FASTQ, BAM, VCF, and parsed CSV files are written.
  --fastq-dir       Directory with paired FASTQ files. Required for the fastp stage.
  --reference       FASTA reference used by BBMap and FreeBayes. Required for map/freebayes.

Optional:
  --sample-glob     R1 FASTQ glob inside --fastq-dir. Default: *_R1*.fastq.gz
                    Example: '*VB_Fn*_R1.fastq.gz'
  --sample-list     Plain text file with sample names to keep, one per line.
                    Names may include or omit the trailing _S sample index.
  --sample-regex    Bash/Python regex used to keep matching samples.
  --steps           Comma-separated stages. Default: fastp,map,subsample,freebayes,parse
                    Valid stages: fastp,map,subsample,freebayes,parse
  --threads         Threads for fastp and BBMap where supported. Default: 4
  --subsample-reads Reads per BAM for the subsample stage. Default: 100000
                    Use 0 to keep all reads.
  --min-dp          Parser depth filter. Default: 10000, matching the original notebook.
  --output-csv      Parsed long variant table. Default: <analysis-dir>/var_file.csv
  --force           Recompute outputs even if stage output files already exist.
  --help            Show this message.

Outputs:
  <analysis-dir>/merged_fastp/
  <analysis-dir>/html_fastp/
  <analysis-dir>/Bam/
  <analysis-dir>/Bam_stats/
  <analysis-dir>/bam_subset/
  <analysis-dir>/vcf/
  <analysis-dir>/var_file.csv
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

fastq_dir=""
analysis_dir=""
reference=""
sample_glob="*_R1*.fastq.gz"
sample_list=""
sample_regex=""
steps="fastp,map,subsample,freebayes,parse"
threads=4
subsample_reads=100000
min_dp=10000
output_csv=""
force=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq-dir)
            fastq_dir="$2"
            shift 2
            ;;
        --analysis-dir)
            analysis_dir="$2"
            shift 2
            ;;
        --reference)
            reference="$2"
            shift 2
            ;;
        --sample-glob)
            sample_glob="$2"
            shift 2
            ;;
        --sample-list)
            sample_list="$2"
            shift 2
            ;;
        --sample-regex)
            sample_regex="$2"
            shift 2
            ;;
        --steps)
            steps="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --subsample-reads)
            subsample_reads="$2"
            shift 2
            ;;
        --min-dp)
            min_dp="$2"
            shift 2
            ;;
        --output-csv)
            output_csv="$2"
            shift 2
            ;;
        --force)
            force=true
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

step_enabled() {
    case ",$steps," in
        *",$1,"*) return 0 ;;
        *) return 1 ;;
    esac
}

if [[ -z "$analysis_dir" ]]; then
    usage >&2
    exit 1
fi

if step_enabled fastp && [[ -z "$fastq_dir" ]]; then
    echo "--fastq-dir is required when the fastp stage is enabled." >&2
    exit 1
fi

if { step_enabled map || step_enabled freebayes; } && [[ -z "$reference" ]]; then
    echo "--reference is required when map or freebayes stages are enabled." >&2
    exit 1
fi

if [[ -n "$fastq_dir" && ! -d "$fastq_dir" ]]; then
    echo "FASTQ directory does not exist: $fastq_dir" >&2
    exit 1
fi

if [[ -n "$reference" && ! -f "$reference" ]]; then
    echo "Reference FASTA does not exist: $reference" >&2
    exit 1
fi

if [[ -n "$sample_list" && ! -f "$sample_list" ]]; then
    echo "Sample list does not exist: $sample_list" >&2
    exit 1
fi

if [[ ! -d "$analysis_dir" ]]; then
    mkdir -p "$analysis_dir"
fi
cd "$analysis_dir"

if step_enabled fastp; then
    mkdir -p merged_fastp
    mkdir -p html_fastp
fi
if step_enabled map; then
    mkdir -p Bam
    mkdir -p Bam_stats
fi
if step_enabled subsample; then
    mkdir -p bam_subset
fi
if step_enabled freebayes; then
    mkdir -p vcf
fi

if [[ -z "$output_csv" ]]; then
    output_csv="$analysis_dir/var_file.csv"
fi

sample_without_sequencing_index() {
    local sample="$1"
    printf '%s\n' "${sample%%_S[0-9]*}"
}

sample_from_read1() {
    local read1="$1"
    local sample
    sample="$(basename "$read1")"
    sample="${sample%.fastq.gz}"
    sample="${sample%.fq.gz}"
    sample="${sample%.fastq}"
    sample="${sample%.fq}"
    sample="${sample%_R1_001}"
    sample="${sample%_R1}"
    sample="${sample%_1}"
    sample="${sample%_L???}"
    printf '%s\n' "$sample"
}

sample_allowed() {
    local sample="$1"
    local short_sample
    short_sample="$(sample_without_sequencing_index "$sample")"

    if [[ -n "$sample_regex" ]]; then
        if [[ ! "$sample" =~ $sample_regex && ! "$short_sample" =~ $sample_regex ]]; then
            return 1
        fi
    fi

    if [[ -n "$sample_list" ]]; then
        if ! grep -Fxq "$sample" "$sample_list" && ! grep -Fxq "$short_sample" "$sample_list"; then
            return 1
        fi
    fi

    return 0
}

paired_read2() {
    local read1="$1"
    local read2="$read1"

    read2="${read2/_R1_/_R2_}"
    if [[ "$read2" == "$read1" ]]; then
        read2="${read2/_R1/_R2}"
    fi
    if [[ "$read2" == "$read1" ]]; then
        read2="${read2/_1.fastq/_2.fastq}"
    fi
    if [[ "$read2" == "$read1" ]]; then
        read2="${read2/_1.fq/_2.fq}"
    fi

    printf '%s\n' "$read2"
}

collect_read1_files() {
    local read1
    read1_files=()

    for read1 in "$fastq_dir"/$sample_glob; do
        [[ -e "$read1" ]] || continue
        read1_files+=("$read1")
    done
}

if step_enabled fastp; then
    ## This stage follows the original 01_fastp_merge.sh command and keeps the
    ## same adapter sequences and length cutoff.
    collect_read1_files
    file_count=0
    total_files="${#read1_files[@]}"

    if [[ "$total_files" -eq 0 ]]; then
        echo "No R1 FASTQ files matched $fastq_dir/$sample_glob" >&2
        exit 1
    fi

    for read1 in "${read1_files[@]}"; do
        sample="$(sample_from_read1 "$read1")"
        if ! sample_allowed "$sample"; then
            continue
        fi

        read2="$(paired_read2 "$read1")"
        if [[ ! -f "$read2" ]]; then
            echo "Missing paired R2 read for $read1; expected $read2" >&2
            exit 1
        fi

        file_count=$((file_count+1))
        merged="merged_fastp/${sample}_merged.fastq"
        html_file="html_fastp/${sample}.html"
        json_file="html_fastp/${sample}.json"

        echo "File $file_count/$total_files merging with fastp: $sample"
        if [[ -f "$merged" && "$force" != true ]]; then
            echo "Skipping existing merged FASTQ: $merged"
            continue
        fi

        fastp \
            --in1="$read1" \
            --in2="$read2" \
            -m \
            --merged_out "$merged" \
            --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --length_required=100 \
            --thread "$threads" \
            --html="$html_file" \
            --json="$json_file"
    done
fi

if step_enabled map; then
    ## Map merged reads using bbmap. The core bbmap options match 02_map.sh.
    merged_files=(merged_fastp/*_merged.fastq)
    file_count=0
    total_files=0
    for file in "${merged_files[@]}"; do
        [[ -e "$file" ]] || continue
        total_files=$((total_files+1))
    done

    if [[ "$total_files" -eq 0 ]]; then
        echo "No merged FASTQ files found in $analysis_dir/merged_fastp" >&2
        exit 1
    fi

    for file in "${merged_files[@]}"; do
        [[ -e "$file" ]] || continue
        filename="$(basename "$file" _merged.fastq)"
        if ! sample_allowed "$filename"; then
            continue
        fi

        file_count=$((file_count+1))
        sorted_bam="Bam/${filename}_sorted.bam"
        echo "File $file_count/$total_files aligning: $filename"
        if [[ -f "$sorted_bam" && "$force" != true ]]; then
            echo "Skipping existing sorted BAM: $sorted_bam"
            continue
        fi

        bbmap.sh \
            in="$file" \
            out="Bam/${filename}.sam" \
            slow \
            k=12 \
            threads="$threads" \
            ref="$reference" \
            scafstats="Bam_stats/${filename}.stats" \
            minratio=0.9

        samtools view -b -S "Bam/${filename}.sam" > "Bam/${filename}.bam"
        samtools sort "Bam/${filename}.sam" > "$sorted_bam"
        samtools index "$sorted_bam"
        rm "Bam/${filename}.sam"
    done
fi

if step_enabled subsample; then
    ## Subsample reads to speed up variant calling. Set --subsample-reads 0 to
    ## keep every read while still copying files into bam_subset/.
    for bam in Bam/*_sorted.bam; do
        [[ -e "$bam" ]] || continue
        filename="$(basename "$bam" _sorted.bam)"
        if ! sample_allowed "$filename"; then
            continue
        fi

        subset="bam_subset/${filename}.subset.bam"
        if [[ -f "$subset" && "$force" != true ]]; then
            echo "Skipping existing subset BAM: $subset"
            continue
        fi

        echo "Processing $bam"
        count="$(samtools view -c "$bam")"

        if [[ "$subsample_reads" -le 0 || "$count" -le "$subsample_reads" ]]; then
            echo "$bam has $count reads. Keeping all reads."
            cp "$bam" "$subset"
        else
            echo "$bam has $count reads. Subsampling to $subsample_reads reads."
            fraction="$(awk "BEGIN {printf \"%.6f\", $subsample_reads/$count}")"
            samtools view -b -s "$fraction" "$bam" > "$subset"
        fi

        samtools index "$subset"
    done
fi

if step_enabled freebayes; then
    ## Call variants using the exact FreeBayes options from
    ## 03_freebayes_var_call.sh. --pooled-continuous is essential for amplicon
    ## sequencing from heterogeneous populations.
    bam_files=(bam_subset/*.subset.bam)
    if [[ ! -e "${bam_files[0]}" ]]; then
        bam_files=(Bam/*_sorted.bam)
    fi

    for bam in "${bam_files[@]}"; do
        [[ -e "$bam" ]] || continue
        filename="$(basename "$bam")"
        filename="${filename%.subset.bam}"
        filename="${filename%_sorted.bam}"
        if ! sample_allowed "$filename"; then
            continue
        fi

        vcf_file="vcf/${filename}.vcf"
        echo "$filename processing with freebayes"
        if [[ -f "$vcf_file" && "$force" != true ]]; then
            echo "Skipping existing VCF: $vcf_file"
            continue
        fi

        freebayes \
            -f "$reference" \
            "$bam" \
            --min-base-quality 3 \
            --haplotype-length 30 \
            -C 1 \
            --min-coverage 50 \
            --min-alternate-fraction 0.001 \
            --pooled-continuous \
            > "$vcf_file"
    done
fi

if step_enabled parse; then
    ## Parse and concatenate FreeBayes VCF files into the long table consumed by
    ## scripts/run_endogenous_amplicon_vcf_to_efficiency.R.
    parse_args=(
        --vcf-dir "$analysis_dir/vcf"
        --output "$output_csv"
        --min-dp "$min_dp"
    )
    if [[ -n "$sample_list" ]]; then
        parse_args+=(--sample-list "$sample_list")
    fi
    if [[ -n "$sample_regex" ]]; then
        parse_args+=(--sample-regex "$sample_regex")
    fi

    python3 "$script_dir/parse_freebayes_vcf.py" "${parse_args[@]}"
fi

echo "Amplicon preprocessing complete."
echo "Parsed variant table: $output_csv"
