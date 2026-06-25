#!/usr/bin/env python3

"""Parse concatenated FreeBayes amplicon VCF calls into a long CSV table.

This script is a command-line replacement for the original
scripts_amplicon/04_parse_vcf.ipynb notebook. It preserves the relevant
notebook behavior:

- split multi-allelic VCF records into one row per ALT allele;
- retain QUAL but do not filter on it;
- calculate frc_alt = AO / DP * 100 and frc_ref = RO / DP * 100;
- keep rows with DP > min_dp, default 10000;
- derive mismatches, positions, and pos_mismatch;
- derive Guide, Experiment, Timepoint, and Replicate from the sample name when
  the original 20250311 naming scheme is used.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


DROP_INFO_COLUMNS = {
    "PAIRED",
    "PAIREDR",
    "PAO",
    "PQA",
    "PRO",
    "RUN",
    "AB",
    "ABP",
    "AC",
    "AF",
    "AN",
    "CIGAR",
    "DPB",
    "DPRA",
    "EPP",
    "EPPR",
    "GTI",
    "LEN",
    "MEANALT",
    "MQM",
    "RPL",
    "RPP",
    "RPPR",
    "RPR",
    "SAF",
    "SAP",
    "SAR",
    "SRF",
    "SRP",
    "SRR",
    "PQR",
    "QA",
    "QR",
    "MQMR",
    "NS",
}


PREFERRED_COLUMNS = [
    "AO",
    "DP",
    "NUMALT",
    "ODDS",
    "RO",
    "TYPE",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUAL",
    "frc_alt",
    "frc_ref",
    "Sample",
    "Guide",
    "Experiment",
    "Timepoint",
    "Replicate",
    "mismatches",
    "positions",
    "pos_mismatch",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parse FreeBayes VCF files into one long amplicon variant table."
    )
    parser.add_argument("--vcf-dir", required=True, help="Directory containing .vcf or .vcf.gz files.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument(
        "--min-dp",
        type=float,
        default=10000,
        help="Keep variants with DP > min_dp. Default: 10000.",
    )
    parser.add_argument(
        "--sample-list",
        default="",
        help="Optional text file with sample names to keep, one per line.",
    )
    parser.add_argument(
        "--sample-regex",
        default="",
        help="Optional regex used to keep matching samples.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def sample_from_vcf_path(path: Path) -> str:
    name = path.name
    if name.endswith(".vcf.gz"):
        name = name[:-7]
    elif name.endswith(".vcf"):
        name = name[:-4]
    return name.split("_S")[0]


def short_sample_name(sample: str) -> str:
    return re.split(r"_S[0-9].*", sample)[0]


def read_sample_list(path: str) -> Optional[Set[str]]:
    if not path:
        return None

    keep: Set[str] = set()
    with open(path, "r") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                keep.add(value)
    return keep


def sample_allowed(sample: str, sample_list: Optional[Set[str]], sample_regex: str) -> bool:
    short_sample = short_sample_name(sample)

    if sample_regex:
        if not re.search(sample_regex, sample) and not re.search(sample_regex, short_sample):
            return False

    if sample_list is not None:
        if sample not in sample_list and short_sample not in sample_list:
            return False

    return True


def parse_info(info: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    if not info or info == ".":
        return parsed

    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = "True"
    return parsed


def allele_value(value: str, allele_index: int, allele_count: int) -> str:
    parts = value.split(",")
    if len(parts) == allele_count:
        return parts[allele_index]
    return value


def as_float(value: object) -> float:
    if value is None:
        return math.nan
    text = str(value)
    if text in {"", "."}:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def format_number(value: float) -> str:
    if math.isnan(value):
        return ""
    if value.is_integer():
        return str(int(value))
    return repr(value)


def find_mismatches(ref: str, alt: str) -> Tuple[str, List[int]]:
    mismatches: List[str] = []
    positions: List[int] = []

    for index, (ref_base, alt_base) in enumerate(zip(ref, alt)):
        if ref_base != alt_base:
            mismatches.append(alt_base)
            positions.append(index)

    return "".join(mismatches), positions


def add_sample_metadata(row: Dict[str, str], sample: str) -> None:
    row["Sample"] = sample

    split_sample = sample.split("_")
    metadata = split_sample[5:]
    keys = ["Guide", "Experiment", "Timepoint", "Replicate"]
    for index, key in enumerate(keys):
        row[key] = metadata[index] if index < len(metadata) else ""


def parse_vcf(path: Path, sample: str) -> Iterable[Dict[str, str]]:
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, _variant_id, ref, alt_field, qual, _filter_value, info = fields[:8]
            alt_values = alt_field.split(",")
            info_values = parse_info(info)

            for allele_index, alt in enumerate(alt_values):
                row: Dict[str, str] = {}
                for key, value in info_values.items():
                    if key in DROP_INFO_COLUMNS:
                        continue
                    row[key] = allele_value(value, allele_index, len(alt_values))

                row["CHROM"] = chrom
                row["POS"] = pos
                row["REF"] = ref
                row["ALT"] = alt
                row["QUAL"] = qual
                add_sample_metadata(row, sample)
                yield row


def finalize_row(row: Dict[str, str], min_dp: float) -> Optional[Dict[str, str]]:
    ao = as_float(row.get("AO"))
    ro = as_float(row.get("RO"))
    dp = as_float(row.get("DP"))
    qual = as_float(row.get("QUAL"))

    if math.isnan(dp) or dp <= min_dp:
        return None

    frc_alt = ao / dp * 100 if not math.isnan(ao) and dp else math.nan
    frc_ref = ro / dp * 100 if not math.isnan(ro) and dp else math.nan
    mismatches, positions = find_mismatches(row.get("REF", ""), row.get("ALT", ""))

    row["AO"] = format_number(ao)
    row["RO"] = format_number(ro)
    row["DP"] = format_number(dp)
    row["QUAL"] = format_number(qual)
    row["frc_alt"] = format_number(frc_alt)
    row["frc_ref"] = format_number(frc_ref)
    row["mismatches"] = mismatches
    row["positions"] = str(positions)

    pos = as_float(row.get("POS"))
    first_position = positions[0] if positions else 0
    row["pos_mismatch"] = "" if math.isnan(pos) else str(int(pos) + first_position)

    return row


def ordered_columns(rows: List[Dict[str, str]]) -> List[str]:
    discovered: List[str] = []
    for row in rows:
        for key in row:
            if key not in discovered:
                discovered.append(key)

    ordered = [key for key in PREFERRED_COLUMNS if key in discovered]
    ordered.extend(key for key in discovered if key not in ordered)
    return ordered


def main() -> None:
    args = parse_args()
    vcf_dir = Path(args.vcf_dir)
    output = Path(args.output)
    sample_list = read_sample_list(args.sample_list)

    if not vcf_dir.is_dir():
        raise SystemExit(f"VCF directory does not exist: {vcf_dir}")

    rows: List[Dict[str, str]] = []
    vcf_paths = sorted(list(vcf_dir.glob("*.vcf")) + list(vcf_dir.glob("*.vcf.gz")))
    for path in vcf_paths:
        sample = sample_from_vcf_path(path)
        if not sample_allowed(sample, sample_list, args.sample_regex):
            continue

        print(path.name)
        for raw_row in parse_vcf(path, sample):
            row = finalize_row(raw_row, args.min_dp)
            if row is not None:
                rows.append(row)

    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ordered_columns(rows)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output}")


if __name__ == "__main__":
    main()
