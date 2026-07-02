# Publication FCS Data

The event-level FCS files used for Figures 1 and 2 are distributed as GitHub
Release assets because the complete package is larger than 1 GiB. The Git
repository contains the code and path-independent sample metadata; it does not
duplicate the binary FCS files in Git history.

Download the assets from the
[`fcs-data-v1` release](https://github.com/vbstagui98/cas12a_precision_editing/releases/tag/fcs-data-v1).

## Release contents

- `cas12a_precision_editing_figure1_fcs.zip`: Cas-variant, guide-promoter, and
  direct-repeat GFP assays used in Figure 1 and associated supplementary plots.
- `cas12a_precision_editing_figure2_fcs.zip`: GFP guide-panel assays used in
  Figure 2 and associated supplementary plots.
- `publication_fcs_manifest.csv`: one row per released FCS file, including the
  publication metadata, byte count, archive name, and SHA-256 checksum.
- `release_checksums.sha256`: SHA-256 checksums for the release assets.

The release contains 1,272 files totaling 1,241,415,160 bytes (1.156 GiB)
before ZIP container overhead: 446 Figure 1 files and 826 Figure 2 files.

Each ZIP extracts to a separate directory containing `FCS_DATA/` and an
archive-specific copy of the manifest. Files are selected by the exact basenames
in `data/publication_inputs/fcs_publication_platemap.csv`; superseded and unused
wells are excluded.

## Run the analysis

After extracting one archive, run:

```bash
Rscript scripts/run_fcs_efficiency.R \
  --fcs-dir path/to/extracted_archive/FCS_DATA \
  --platemap data/publication_inputs/fcs_publication_platemap.csv \
  --output results/fcs_efficiency
```

The full plate map can be used with either archive. Rows without a matching FCS
basename are ignored because the workflow starts from the files present in the
selected `FCS_DATA/` directory.

## Rebuild the release assets

The release builder searches one or more local archives, selects only committed
publication basenames, verifies duplicate copies, computes checksums, and fails
if any required file is absent or is an online-only zero-byte placeholder.

```bash
python3 scripts/build_fcs_release.py \
  --search-root path/to/local/fcs_archive_1 \
  --search-root path/to/local/fcs_archive_2 \
  --output-dir build/fcs_release
```

Use `--check-only` to validate completeness without writing the ZIP files.
