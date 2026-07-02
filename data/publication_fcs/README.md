# Publication FCS files

The event-level FCS files used for Figures 1 and 2 are distributed as two
GitHub release archives:

- `cas12a_precision_editing_figure1_fcs.zip`: 446 files used for the
  Cas-variant, guide-promoter, and direct-repeat experiments.
- `cas12a_precision_editing_figure2_fcs.zip`: 826 files used for the GFP guide
  panel.

Download both files from the
[`fcs-data-v1` release](https://github.com/vsbatagui/cas12a_precision_editing/releases/tag/fcs-data-v1).
The sample annotation is in
`data/publication_inputs/fcs_publication_platemap.csv`.

After extracting an archive, run:

```bash
Rscript scripts/run_fcs_efficiency.R \
  --fcs-dir path/to/FCS_DATA \
  --platemap data/publication_inputs/fcs_publication_platemap.csv \
  --output results/fcs_efficiency
```
