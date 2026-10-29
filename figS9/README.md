# Fig. S9

This folder contains the data and analysis scripts for the Extended Data Fig. 9 panel-wise analyses. The count tables are stored as CSV files and read by scripts so the data are not duplicated in code.

## Files

- `data/figS9a_overpenetration_by_line_counts.csv`: per-pistil `normal`/`abnormal` counts for panel a over-penetration (`ec1` vs independent `NM` lines).
- `data/figS9b_developing_seed_by_line.csv`: raw seed-set counts for panel b developing-seed analysis (`ec1` vs independent `NM` lines).
- `scripts/figS9a_overpenetration_by_line.R`: panel a plot/statistics script using one-model logistic pairwise contrasts with BH correction across all six pairwise comparisons.
- `scripts/figS9b_developing_seed_by_line.R`: panel b plot/statistics script using one-model logistic pairwise contrasts with BH correction across all six pairwise comparisons.
- `plots/figS9a_overpenetration_by_line.pdf`: panel a summary plot.
- `plots/figS9b_developing_seed_by_line.pdf`: panel b summary plot.
- `results/figS9a_overpenetration_by_line_summary.csv`: panel a quartile/median summary statistics.
- `results/figS9a_overpenetration_by_line_pvalues.csv`: panel a pairwise p-values for annotations.
- `results/figS9b_developing_seed_by_line_summary.csv`: panel b quartile/median summary statistics.
- `results/figS9b_developing_seed_by_line_pvalues.csv`: panel b pairwise p-values for annotations.

## Notes

The datasets in `figS9a_overpenetration_by_line_counts.csv` and `figS9b_developing_seed_by_line.csv` are source-of-truth data for panel a and panel b respectively. This preserves panel-specific data and analysis logic in a transparent, reproducible form for repository storage.
