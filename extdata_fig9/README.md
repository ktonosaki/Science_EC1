# Extended Data Fig. 9

This folder contains the data and analysis scripts for the Extended Data Fig. 9 panel-wise analyses. The count tables are stored as CSV files and read by scripts so the data are not duplicated in code.

## Files

- `data/extdata_fig9_seed_by_line.csv`: raw seed-set counts for panel a (`ec1` vs independent `NM` lines).
- `data/extdata_fig9_abnormal_by_line_counts.csv`: per-pistil `normal`/`abnormal` counts for panel b (`ec1` vs independent `NM` lines).
- `scripts/extdata_fig9_seed_by_line.R`: panel a plot/statistics script using pairwise logistic LRT (+ BH correction).
- `scripts/extdata_fig9_abnormal_by_line.R`: panel b plot/statistics script using pairwise logistic LRT.
- `plots/extdata_fig9_seed_by_line_summary.pdf`: panel a summary plot.
- `plots/extdata_fig9_abnormal_by_line_summary.pdf`: panel b summary plot.
- `results/extdata_fig9_seed_by_line_summary.csv`: panel a quartile/median summary statistics.
- `results/extdata_fig9_seed_by_line_pvalues.csv`: panel a pairwise p-values for annotations.
- `results/extdata_fig9_abnormal_by_line_summary.csv`: panel b quartile/median summary statistics.
- `results/extdata_fig9_abnormal_by_line_pvalues.csv`: panel b pairwise p-values for annotations.

## Notes

The datasets in `extdata_fig9_seed_by_line.csv` and `extdata_fig9_abnormal_by_line_counts.csv` are source-of-truth data for panel a and panel b respectively. This preserves panel-specific data and analysis logic in a transparent, reproducible form for repository storage.
