# Figure 3

This folder contains the raw data and analysis scripts used for Figure 3. The count table is stored as a CSV and is read by the analysis script so the data are not duplicated in code.

## Files

- `data/fig3_abnormal_ovule.csv`: counts for abnormal ovule phenotypes reported in the figure analysis.
- `scripts/fig3_logistic.R`: logistic model script used for the main Figure 3 analysis.
- `results/fig3_abnormal_ovule_summary.csv`: quartile and median summary statistics used for the panel overlays.
- `results/fig3_abnormal_ovule_pvalues.csv`: model-based and Fisher exact test p-values used for the figure annotation labels.

## Data provenance

The observation counts in `fig3_abnormal_ovule.csv` are the source-of-truth values used by the analysis script.
