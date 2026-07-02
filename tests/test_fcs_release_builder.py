#!/usr/bin/env python3

import csv
import importlib.util
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_fcs_release", ROOT / "scripts" / "build_fcs_release.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FcsReleaseBuilderTest(unittest.TestCase):
    def write_platemap(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["name", "figure", "paper_panel"]
            )
            writer.writeheader()
            writer.writerows(rows)

    def test_builds_figure_archives_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            (source / "figure1.fcs").write_bytes(b"FCS3.1-figure-1")
            (source / "figure2.fcs").write_bytes(b"FCS3.1-figure-2")
            plate = root / "plate.csv"
            self.write_platemap(
                plate,
                [
                    {"name": "figure1.fcs", "figure": "Figure 1", "paper_panel": "A"},
                    {"name": "figure2.fcs", "figure": "Figure 2", "paper_panel": "B"},
                ],
            )

            rows, columns = MODULE.read_platemap(plate)
            selected, identity = MODULE.resolve_files(
                rows, MODULE.index_fcs_files([source])
            )
            manifest, manifest_columns = MODULE.build_manifest(
                rows, columns, identity
            )
            outputs = MODULE.write_archives(
                root / "output", manifest, manifest_columns, selected
            )

            self.assertEqual(4, len(outputs))
            figure1 = root / "output" / MODULE.ARCHIVE_NAMES["Figure 1"]
            with zipfile.ZipFile(figure1) as archive:
                names = archive.namelist()
            self.assertTrue(any(name.endswith("/FCS_DATA/figure1.fcs") for name in names))
            self.assertFalse(any(name.endswith("/FCS_DATA/figure2.fcs") for name in names))

    def test_rejects_zero_byte_placeholder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "empty.fcs").touch()
            rows = [{"name": "empty.fcs"}]
            with self.assertRaisesRegex(FileNotFoundError, "zero-byte"):
                MODULE.resolve_files(rows, MODULE.index_fcs_files([root]))

    def test_rejects_conflicting_duplicate_basename(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "first"
            second = root / "second"
            first.mkdir()
            second.mkdir()
            (first / "sample.fcs").write_bytes(b"one")
            (second / "sample.fcs").write_bytes(b"two")
            rows = [{"name": "sample.fcs"}]
            with self.assertRaisesRegex(ValueError, "conflicting"):
                MODULE.resolve_files(rows, MODULE.index_fcs_files([first, second]))


if __name__ == "__main__":
    unittest.main()
