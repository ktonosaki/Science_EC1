# Extended Data Fig. 9

This folder contains the data and analysis script for the individual NM-line analysis corresponding to the Extended Data Fig. 9 description in the bioRxiv article. The count table is stored as a CSV and is read by the analysis script so the data are not duplicated in code.

## Files

- `data/extdata_fig9_seed_by_line.csv`: raw seed set counts for the three independent `Sup35NM` replacement lines and the `ec1 null` control.
- `scripts/extdata_fig9_seed_by_line.R`: plot generation and significance testing script for the individual-line comparison.

## Notes

The dataset in `extdata_fig9_seed_by_line.csv` is the source-of-truth data used by the analysis script. This preserves the figure-specific data and analysis logic in a transparent, reproducible form for repository storage.
