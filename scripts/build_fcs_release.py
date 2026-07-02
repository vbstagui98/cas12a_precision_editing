#!/usr/bin/env python3
"""Build publication-scoped FCS release archives from the committed plate map."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
import sys
import zipfile
from collections import defaultdict
from pathlib import Path


MANIFEST_NAME = "publication_fcs_manifest.csv"
ARCHIVE_NAMES = {
    "Figure 1": "cas12a_precision_editing_figure1_fcs.zip",
    "Figure 2": "cas12a_precision_editing_figure2_fcs.zip",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build checksummed Figure 1 and Figure 2 FCS release archives."
    )
    parser.add_argument(
        "--platemap",
        default="data/publication_inputs/fcs_publication_platemap.csv",
        type=Path,
        help="Publication FCS plate map (default: %(default)s).",
    )
    parser.add_argument(
        "--search-root",
        action="append",
        required=True,
        type=Path,
        help="Root to search recursively for FCS files; may be repeated.",
    )
    parser.add_argument(
        "--output-dir",
        default="build/fcs_release",
        type=Path,
        help="Directory for manifests and ZIP archives (default: %(default)s).",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate files without writing release assets.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_platemap(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = list(reader.fieldnames or [])

    required = {"name", "figure", "paper_panel"}
    missing = sorted(required.difference(columns))
    if missing:
        raise ValueError(f"Plate map is missing columns: {', '.join(missing)}")
    if not rows:
        raise ValueError("Plate map contains no samples.")

    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        counts[row["name"]] += 1
    duplicates = sorted(name for name, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(
            "Plate map must contain one row per unique FCS basename; duplicated: "
            + ", ".join(duplicates[:10])
        )

    unsupported = sorted({row["figure"] for row in rows}.difference(ARCHIVE_NAMES))
    if unsupported:
        raise ValueError("Unsupported figure labels: " + ", ".join(unsupported))
    return rows, columns


def index_fcs_files(roots: list[Path]) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = defaultdict(list)
    for root in roots:
        if not root.is_dir():
            raise FileNotFoundError(f"Search root does not exist: {root}")
        for current, _, files in os.walk(root):
            for filename in files:
                if filename.lower().endswith(".fcs"):
                    index[filename].append(Path(current, filename))
    return index


def resolve_files(
    rows: list[dict[str, str]], index: dict[str, list[Path]]
) -> tuple[dict[str, Path], dict[str, tuple[int, str]]]:
    selected: dict[str, Path] = {}
    identity: dict[str, tuple[int, str]] = {}
    unresolved: list[str] = []
    conflicting: list[str] = []

    for row in rows:
        name = row["name"]
        candidates = [path for path in index.get(name, []) if path.stat().st_size > 0]
        if not candidates:
            unresolved.append(name)
            continue

        candidates.sort(key=lambda path: str(path))
        identities: dict[tuple[int, str], list[Path]] = defaultdict(list)
        for path in candidates:
            key = (path.stat().st_size, sha256(path))
            identities[key].append(path)

        if len(identities) != 1:
            details = "; ".join(
                f"{size} bytes, sha256={digest}: {paths[0]}"
                for (size, digest), paths in identities.items()
            )
            conflicting.append(f"{name} ({details})")
            continue

        file_identity = next(iter(identities))
        selected[name] = candidates[0]
        identity[name] = file_identity

    if unresolved:
        preview = "\n  ".join(unresolved[:20])
        suffix = "\n  ..." if len(unresolved) > 20 else ""
        raise FileNotFoundError(
            f"{len(unresolved)} required FCS files are missing or zero-byte placeholders:\n  "
            f"{preview}{suffix}"
        )
    if conflicting:
        preview = "\n  ".join(conflicting[:10])
        suffix = "\n  ..." if len(conflicting) > 10 else ""
        raise ValueError(
            f"{len(conflicting)} FCS basenames resolve to conflicting non-empty files:\n  "
            f"{preview}{suffix}"
        )
    return selected, identity


def build_manifest(
    rows: list[dict[str, str]],
    columns: list[str],
    identity: dict[str, tuple[int, str]],
) -> tuple[list[dict[str, str]], list[str]]:
    metadata_columns = [
        column for column in columns if column != "raw_fcs_available_locally"
    ]
    output_columns = [
        "archive",
        "file_bytes",
        "sha256",
        "included_in_release",
        *metadata_columns,
    ]
    output_rows = []
    for row in rows:
        size, digest = identity[row["name"]]
        output_rows.append(
            {
                "archive": ARCHIVE_NAMES[row["figure"]],
                "file_bytes": str(size),
                "sha256": digest,
                "included_in_release": "True",
                **{column: row[column] for column in metadata_columns},
            }
        )
    return output_rows, output_columns


def csv_bytes(rows: list[dict[str, str]], columns: list[str]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def write_archives(
    output_dir: Path,
    manifest_rows: list[dict[str, str]],
    manifest_columns: list[str],
    selected: dict[str, Path],
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / MANIFEST_NAME
    manifest_path.write_bytes(csv_bytes(manifest_rows, manifest_columns))

    outputs = [manifest_path]
    for figure, archive_name in ARCHIVE_NAMES.items():
        archive_rows = [row for row in manifest_rows if row["figure"] == figure]
        archive_path = output_dir / archive_name
        root = archive_name.removesuffix(".zip")
        with zipfile.ZipFile(
            archive_path, mode="w", compression=zipfile.ZIP_STORED, allowZip64=True
        ) as archive:
            archive.writestr(
                f"{root}/{MANIFEST_NAME}",
                csv_bytes(archive_rows, manifest_columns),
            )
            archive.writestr(
                f"{root}/README.txt",
                (
                    f"{figure} publication FCS files\n\n"
                    "FCS_DATA contains only event-level files referenced by the "
                    "publication plate map. Verify each file against "
                    f"{MANIFEST_NAME} before analysis.\n"
                ),
            )
            for row in archive_rows:
                archive.write(
                    selected[row["name"]],
                    arcname=f"{root}/FCS_DATA/{row['name']}",
                )
        outputs.append(archive_path)

    checksums = output_dir / "release_checksums.sha256"
    with checksums.open("w", encoding="ascii") as handle:
        for path in outputs:
            handle.write(f"{sha256(path)}  {path.name}\n")
    outputs.append(checksums)
    return outputs


def main() -> int:
    args = parse_args()
    try:
        rows, columns = read_platemap(args.platemap)
        index = index_fcs_files(args.search_root)
        selected, identity = resolve_files(rows, index)
        manifest_rows, manifest_columns = build_manifest(rows, columns, identity)
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    total_bytes = sum(size for size, _ in identity.values())
    print(f"Resolved {len(selected)} publication FCS files.")
    print(f"Total uncompressed size: {total_bytes} bytes ({total_bytes / 1024**3:.3f} GiB).")
    for figure in ARCHIVE_NAMES:
        count = sum(row["figure"] == figure for row in rows)
        print(f"{figure}: {count} files")

    if args.check_only:
        return 0

    outputs = write_archives(
        args.output_dir, manifest_rows, manifest_columns, selected
    )
    for path in outputs:
        print(f"Wrote {path} ({path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
