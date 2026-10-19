# Figure 5

This folder contains the raw data and analysis script used for the Fig. 5 seed-normalization analysis. The count table is stored as a CSV and is read by the analysis script so the data are not duplicated in code.

## Files

- `data/fig5f_normal_seed.csv`: raw counts used to quantify normal seed formation across genotypes (Fig. 5f).
- `data/fig5e_abnormal_ovule_counts.csv`: per-pistil `normal`/`abnormal` counts used for Fig. 5e.
- `scripts/fig5f_normal_seed.R`: plotting and summary script for Fig. 5f genotype comparisons and significance annotations.
- `scripts/fig5e_abnormal_ovule.R`: abnormal-ovule plotting script corresponding to Fig. 5e.
- `plots/fig5f_normal_seed.pdf`: publication-style Fig. 5f summary plot.
- `plots/fig5e_abnormal_ovule.pdf`: Fig. 5e abnormal-ovule summary plot.
- `results/fig5f_normal_seed_summary.csv`: quartile and median summary statistics used for the Fig. 5f overlays.
- `results/fig5f_normal_seed_pvalues.csv`: pairwise logistic LRT results used for Fig. 5f annotations.
- `results/fig5e_abnormal_ovule_summary.csv`: quartile and median summary statistics for Fig. 5e.
- `results/fig5e_abnormal_ovule_pvalues.csv`: pairwise logistic-model contrast results (BH-adjusted) used for Fig. 5e annotations.

## Notes

The dataset in `fig5f_normal_seed.csv` is the source-of-truth data for Fig. 5f seed-set analysis, and `fig5e_abnormal_ovule_counts.csv` is the source-of-truth data used by the Fig. 5e abnormal-ovule script. This archive keeps the figure-specific raw values and analysis logic in one place for transparent reproducibility.
