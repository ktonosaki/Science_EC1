# Figure 3

This folder contains the raw data and analysis scripts used for Figure 3. The count table is stored as a CSV and is read by the analysis script so the data are not duplicated in code.

## Files

- `data/fig3f_abnormal_ovule.csv`: counts for abnormal ovule phenotypes reported in Fig. 3f.
- `data/fig3g_overpenetration_counts.csv`: per-pistil total/abnormal counts used for Fig. 3g over-penetration analysis.
- `scripts/fig3f_abnormal_ovule.R`: logistic model script for Fig. 3f (Middle of rSY and Proximal to EC).
- `scripts/fig3g_overpenetration.R`: Fig. 3g script using binomial logistic regression (LRT) on count data.
- `plots/fig3f_abnormal_ovule.pdf`: publication-style summary plot for Fig. 3f.
- `plots/fig3g_overpenetration.pdf`: publication-style Fig. 3g over-penetration plot.
- `results/fig3f_abnormal_ovule_summary.csv`: quartile and median summary statistics used for the Fig. 3f panel overlays.
- `results/fig3f_abnormal_ovule_pvalues.csv`: model-based and Fisher exact test p-values used for Fig. 3f annotation labels.
- `results/fig3g_overpenetration_summary.csv`: quartile and median summary statistics for Fig. 3g.
- `results/fig3g_overpenetration_pvalues.csv`: p-value table for Fig. 3g.

## Data provenance

The observation counts in `fig3f_abnormal_ovule.csv` are the source-of-truth values used by the Fig. 3f logistic analysis. The Fig. 3g script uses per-pistil count data from `fig3g_overpenetration_counts.csv`, and p-values are reported from the logistic LRT output.
