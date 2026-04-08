#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/score_cas12a_dr_perturbation_mxfold2.sh INPUT.csv OUTPUT.tsv [DR_LEN]

Input assumptions:
  - CSV with a header row
  - column 1 = guide_id
  - column 2 = DR+spacer sequence

Output:
  - tab-delimited table with one row per guide
  - score = 1 only if the DR is still structured and no DR nucleotide pairs to
    the spacer in the MXfold2 prediction

Optional environment variables:
  - MXFOLD2_BIN: explicit path to the mxfold2 executable
  - MXFOLD2_MODEL: passed as '--model ...' to mxfold2 predict
  - MXFOLD2_PARAM: passed as '--param ...' to mxfold2 predict
  - MIN_DR_PAIRED_NT: minimum number of DR nucleotides paired within the DR to
    call the DR "structured" (default: 2)
  - EXPECTED_DR_PREFIX: canonical DR dot-bracket motif to enforce
    (default: .....(((((....)))))
EOF
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 1
fi

input_csv=$1
output_tsv=$2
dr_len=${3:-20}
min_dr_paired_nt=${MIN_DR_PAIRED_NT:-2}
expected_dr_prefix=${EXPECTED_DR_PREFIX:-.....(((((....)))))}

if [[ ! -f "$input_csv" ]]; then
  echo "Input file not found: $input_csv" >&2
  exit 1
fi

if ! [[ "$dr_len" =~ ^[0-9]+$ ]] || (( dr_len < 1 )); then
  echo "DR_LEN must be a positive integer." >&2
  exit 1
fi

if ! [[ "$min_dr_paired_nt" =~ ^[0-9]+$ ]]; then
  echo "MIN_DR_PAIRED_NT must be a non-negative integer." >&2
  exit 1
fi

if [[ ! "$expected_dr_prefix" =~ ^[().]+$ ]]; then
  echo "EXPECTED_DR_PREFIX must contain only dot-bracket characters: . ( )" >&2
  exit 1
fi

mxfold2_bin=${MXFOLD2_BIN:-}
if [[ -z "$mxfold2_bin" ]]; then
  if command -v mxfold2 >/dev/null 2>&1; then
    mxfold2_bin=$(command -v mxfold2)
  elif [[ -x "/Users/u0174312/miniconda3/envs/mxfold2/bin/mxfold2" ]]; then
    mxfold2_bin="/Users/u0174312/miniconda3/envs/mxfold2/bin/mxfold2"
  else
    echo "mxfold2 not found. Set MXFOLD2_BIN or activate the mxfold2 env first." >&2
    exit 1
  fi
fi

score_structure() {
  local structure=$1
  local dr_len_local=$2

  awk -v structure="$structure" -v dr_len="$dr_len_local" '
    BEGIN {
      n = split(structure, chars, "")
      top = 0

      for (i = 1; i <= n; i++) {
        if (chars[i] == "(") {
          stack[++top] = i
        } else if (chars[i] == ")") {
          if (top > 0) {
            j = stack[top--]
            partner[i] = j
            partner[j] = i
          }
        }
      }

      dr_to_dr = 0
      dr_to_spacer = 0
      dr_unpaired = 0

      for (i = 1; i <= dr_len; i++) {
        if (!(i in partner)) {
          dr_unpaired++
        } else if (partner[i] <= dr_len) {
          dr_to_dr++
        } else {
          dr_to_spacer++
        }
      }

      printf "%d\t%d\t%d", dr_to_dr, dr_to_spacer, dr_unpaired
    }
  '
}

build_expected_structure() {
  local seq_len=$1
  local prefix_len=${#expected_dr_prefix}
  local suffix_len

  if (( seq_len < prefix_len )); then
    printf 'NA'
    return
  fi

  suffix_len=$((seq_len - prefix_len))
  printf '%s' "$expected_dr_prefix"
  awk -v n="$suffix_len" 'BEGIN { for (i = 0; i < n; i++) printf "." }'
}

tmp_fasta=$(mktemp)
tmp_meta=$(mktemp)
tmp_pred=$(mktemp)
cleanup() {
  rm -f "$tmp_fasta" "$tmp_meta" "$tmp_pred"
}
trap cleanup EXIT

line_no=0
entry_count=0
while IFS=, read -r raw_id raw_seq _ || [[ -n "${raw_id:-}" || -n "${raw_seq:-}" ]]; do
  ((line_no += 1))

  raw_id=${raw_id//$'\r'/}
  raw_seq=${raw_seq//$'\r'/}
  raw_id_lc=$(printf '%s' "$raw_id" | tr '[:upper:]' '[:lower:]')

  if (( line_no == 1 )) && [[ "$raw_id_lc" == "guide_id" ]]; then
    continue
  fi

  guide_id=${raw_id#\"}
  guide_id=${guide_id%\"}
  guide_seq=${raw_seq#\"}
  guide_seq=${guide_seq%\"}
  guide_seq=${guide_seq// /}

  if [[ -z "$guide_id" && -z "$guide_seq" ]]; then
    continue
  fi

  if [[ -z "$guide_id" || -z "$guide_seq" ]]; then
    echo "Skipping malformed line $line_no in $input_csv" >&2
    continue
  fi

  guide_seq=$(printf '%s' "$guide_seq" | tr '[:lower:]' '[:upper:]')
  guide_seq_rna=$(printf '%s' "$guide_seq" | tr 'T' 'U')

  if (( ${#guide_seq_rna} < dr_len )); then
    echo "Skipping $guide_id: sequence shorter than DR_LEN=$dr_len" >&2
    continue
  fi

  if [[ ! "$guide_seq_rna" =~ ^[ACGUN]+$ ]]; then
    echo "Skipping $guide_id: sequence contains characters outside A/C/G/T/U/N" >&2
    continue
  fi

  ((entry_count += 1))
  fasta_id=$(printf 'seq_%06d' "$entry_count")
  printf '>%s\n%s\n' "$fasta_id" "$guide_seq_rna" >> "$tmp_fasta"
  printf '%s\t%s\t%s\n' "$guide_id" "$guide_seq" "$guide_seq_rna" >> "$tmp_meta"
done < "$input_csv"

mkdir -p "$(dirname "$output_tsv")"
printf '%s\n' \
  "guide_id	guide_seq_input	guide_seq_rna	model	expected_structure	structure_energy	structure	canonical_structure_intact	canonical_structure_disturbed	dr_len	dr_nt_paired_to_dr	dr_nt_paired_to_spacer	dr_nt_unpaired	dr_structure_detected	dr_perturbed	score" \
  > "$output_tsv"

if (( entry_count == 0 )); then
  echo "No valid guides found in $input_csv" >&2
  exit 0
fi

cmd=("$mxfold2_bin" "predict")
if [[ -n "${MXFOLD2_MODEL:-}" ]]; then
  cmd+=("--model" "$MXFOLD2_MODEL")
fi
if [[ -n "${MXFOLD2_PARAM:-}" ]]; then
  cmd+=("--param" "$MXFOLD2_PARAM")
fi
cmd+=("$tmp_fasta")

"${cmd[@]}" | awk '
  /^>/ {
    if (state != 0) {
      exit 1
    }
    state = 1
    next
  }
  state == 1 {
    seq = $0
    state = 2
    next
  }
  state == 2 {
    structure = $1
    energy = $0
    sub(/^[^[:space:]]+[[:space:]]+\(/, "", energy)
    sub(/\)[[:space:]]*$/, "", energy)
    print seq "\t" structure "\t" energy
    state = 0
  }
' > "$tmp_pred"

meta_lines=$(wc -l < "$tmp_meta" | tr -d '[:space:]')
pred_lines=$(wc -l < "$tmp_pred" | tr -d '[:space:]')
if [[ "$meta_lines" != "$pred_lines" ]]; then
  echo "Could not parse mxfold2 output cleanly: expected $meta_lines records, got $pred_lines" >&2
  exit 1
fi

exec 3< "$tmp_meta"
exec 4< "$tmp_pred"

while IFS=$'\t' read -r guide_id guide_seq guide_seq_rna <&3 && IFS=$'\t' read -r pred_seq structure energy <&4; do
  if [[ "$guide_seq_rna" != "$pred_seq" ]]; then
    echo "Sequence mismatch while parsing mxfold2 output for $guide_id" >&2
    exit 1
  fi

  expected_structure=$(build_expected_structure "${#guide_seq_rna}")

  score_fields=$(score_structure "$structure" "$dr_len")
  IFS=$'\t' read -r dr_to_dr dr_to_spacer dr_unpaired <<< "$score_fields"

  if (( dr_to_dr >= min_dr_paired_nt )); then
    dr_structure_detected=yes
  else
    dr_structure_detected=no
  fi

  if (( dr_to_spacer > 0 )); then
    dr_perturbed=yes
  else
    dr_perturbed=no
  fi

  if [[ "$expected_structure" == "NA" ]]; then
    canonical_structure_intact=na
    canonical_structure_disturbed=na
  elif [[ "$structure" == "$expected_structure" ]]; then
    canonical_structure_intact=yes
    canonical_structure_disturbed=no
  else
    canonical_structure_intact=no
    canonical_structure_disturbed=yes
  fi

  if (( dr_to_spacer == 0 )) && [[ "$dr_structure_detected" == "yes" ]]; then
    score=1
  else
    score=0
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$guide_id" \
    "$guide_seq" \
    "$guide_seq_rna" \
    "${MXFOLD2_MODEL:-default}" \
    "$expected_structure" \
    "$energy" \
    "$structure" \
    "$canonical_structure_intact" \
    "$canonical_structure_disturbed" \
    "$dr_len" \
    "$dr_to_dr" \
    "$dr_to_spacer" \
    "$dr_unpaired" \
    "$dr_structure_detected" \
    "$dr_perturbed" \
    "$score" \
    >> "$output_tsv"
done

exec 3<&-
exec 4<&-

echo "Wrote $output_tsv" >&2
