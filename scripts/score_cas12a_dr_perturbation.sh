#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/score_cas12a_dr_perturbation.sh INPUT.csv OUTPUT.tsv [DR_LEN]

Input assumptions:
  - CSV with a header row
  - column 1 = guide_id
  - column 2 = DR+spacer sequence

What this script folds:
  leader transcript + (DR+spacer guide sequence)

Output:
  - tab-delimited table with one row per guide
  - dr_perturbed_by_spacer = whether any DR nucleotide pairs into the spacer
  - canonical_dr_intact = whether the DR window exactly matches the expected
    dot-bracket motif
  - score_generic = 1 if the DR is internally structured and not paired to the spacer
  - score_canonical = 1 if the canonical DR motif is intact and not paired to the spacer

Optional environment variables:
  - LEADER_SEQUENCE: DNA/RNA sequence prepended before the guide
  - EXPECTED_DR_STRUCTURE: exact DR dot-bracket motif to enforce
  - EXPECTED_DR_PREFIX: backward-compatible alias for EXPECTED_DR_STRUCTURE
  - MIN_DR_PAIRED_NT: minimum number of DR nucleotides paired within the DR to
    call the DR "structured" (default: 2)
  - If DR_LEN is omitted, it defaults to the length of EXPECTED_DR_STRUCTURE
EOF
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 1
fi

input_csv=$1
output_tsv=$2

leader_seq_raw=${LEADER_SEQUENCE:-agattttgtagtgccctcttgggctagcggtaaaggtgcgcattttttcacaccctacaatgttctgttcaaaagattttggtcaaacgctgtagaagtgaaagttggtgcgcatgtttcggcgttcgaaacttctccgcagtgaaagataaatgatc}
expected_dr_structure=${EXPECTED_DR_STRUCTURE:-${EXPECTED_DR_PREFIX:-.....(((((....)))))}}
dr_len=${3:-${#expected_dr_structure}}
min_dr_paired_nt=${MIN_DR_PAIRED_NT:-2}

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

if [[ ! "$expected_dr_structure" =~ ^[().]+$ ]]; then
  echo "EXPECTED_DR_STRUCTURE must contain only dot-bracket characters: . ( )" >&2
  exit 1
fi

leader_seq=$(printf '%s' "$leader_seq_raw" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')
leader_seq_rna=$(printf '%s' "$leader_seq" | tr 'T' 'U')

if [[ ! "$leader_seq_rna" =~ ^[ACGUN]+$ ]]; then
  echo "LEADER_SEQUENCE contains characters outside A/C/G/T/U/N" >&2
  exit 1
fi

leader_len=${#leader_seq_rna}

if ! command -v RNAfold >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1 && [[ -z "${SCORE_DR_CONDA_RERUN:-}" ]]; then
    exec env \
      SCORE_DR_CONDA_RERUN=1 \
      LEADER_SEQUENCE="$leader_seq_raw" \
      EXPECTED_DR_STRUCTURE="$expected_dr_structure" \
      MIN_DR_PAIRED_NT="$min_dr_paired_nt" \
      conda run -n vienna bash "$0" "$@"
  fi
  echo "RNAfold not found in PATH. Activate the 'vienna' conda env first." >&2
  exit 1
fi

score_dr_window() {
  local structure=$1
  local dr_start=$2
  local dr_end=$3

  awk -v structure="$structure" -v dr_start="$dr_start" -v dr_end="$dr_end" '
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

      dr_to_leader = 0
      dr_to_dr = 0
      dr_to_spacer = 0
      dr_unpaired = 0

      for (i = dr_start; i <= dr_end; i++) {
        if (!(i in partner)) {
          dr_unpaired++
        } else if (partner[i] < dr_start) {
          dr_to_leader++
        } else if (partner[i] <= dr_end) {
          dr_to_dr++
        } else {
          dr_to_spacer++
        }
      }

      printf "%d\t%d\t%d\t%d", dr_to_leader, dr_to_dr, dr_to_spacer, dr_unpaired
    }
  '
}

mkdir -p "$(dirname "$output_tsv")"

printf '%s\n' \
  "guide_id	guide_seq_input	guide_seq_rna	leader_seq_rna	full_prec_rna	expected_dr_structure	mfe_kcal_mol	full_structure	dr_region_structure	leader_len	dr_start	dr_end	dr_nt_paired_to_leader	dr_nt_paired_to_dr	dr_nt_paired_to_spacer	dr_nt_unpaired	dr_structure_detected	canonical_dr_intact	canonical_dr_disturbed	dr_perturbed_by_spacer	score_generic	score_canonical" \
  > "$output_tsv"

line_no=0
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

  full_prec_rna="${leader_seq_rna}${guide_seq_rna}"
  dr_start=$((leader_len + 1))
  dr_end=$((leader_len + dr_len))

  fold_output=$(printf '%s\n' "$full_prec_rna" | RNAfold --noPS 2>/dev/null)
  structure_line=$(printf '%s\n' "$fold_output" | sed -n '2p')
  full_structure=$(printf '%s\n' "$structure_line" | awk '{print $1}')
  mfe=$(printf '%s\n' "$structure_line" | sed -E 's/.*\(([[:space:]]*[-+]?[0-9]+(\.[0-9]+)?)[[:space:]]*\).*/\1/' | tr -d '[:space:]')

  if [[ -z "$full_structure" || -z "$mfe" ]]; then
    echo "Skipping $guide_id: could not parse RNAfold output" >&2
    continue
  fi

  dr_region_structure=${full_structure:$((dr_start - 1)):dr_len}
  score_fields=$(score_dr_window "$full_structure" "$dr_start" "$dr_end")
  IFS=$'\t' read -r dr_to_leader dr_to_dr dr_to_spacer dr_unpaired <<< "$score_fields"

  if (( dr_to_dr >= min_dr_paired_nt )); then
    dr_structure_detected=yes
  else
    dr_structure_detected=no
  fi

  if (( dr_to_spacer > 0 )); then
    dr_perturbed_by_spacer=yes
  else
    dr_perturbed_by_spacer=no
  fi

  if (( dr_len != ${#expected_dr_structure} )); then
    canonical_dr_intact=na
    canonical_dr_disturbed=na
  elif [[ "$dr_region_structure" == "$expected_dr_structure" ]]; then
    canonical_dr_intact=yes
    canonical_dr_disturbed=no
  else
    canonical_dr_intact=no
    canonical_dr_disturbed=yes
  fi

  if (( dr_to_spacer == 0 )) && [[ "$dr_structure_detected" == "yes" ]]; then
    score_generic=1
  else
    score_generic=0
  fi

  if (( dr_to_spacer == 0 )) && [[ "$canonical_dr_intact" == "yes" ]]; then
    score_canonical=1
  else
    score_canonical=0
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$guide_id" \
    "$guide_seq" \
    "$guide_seq_rna" \
    "$leader_seq_rna" \
    "$full_prec_rna" \
    "$expected_dr_structure" \
    "$mfe" \
    "$full_structure" \
    "$dr_region_structure" \
    "$leader_len" \
    "$dr_start" \
    "$dr_end" \
    "$dr_to_leader" \
    "$dr_to_dr" \
    "$dr_to_spacer" \
    "$dr_unpaired" \
    "$dr_structure_detected" \
    "$canonical_dr_intact" \
    "$canonical_dr_disturbed" \
    "$dr_perturbed_by_spacer" \
    "$score_generic" \
    "$score_canonical" \
    >> "$output_tsv"
done < "$input_csv"

echo "Wrote $output_tsv" >&2
