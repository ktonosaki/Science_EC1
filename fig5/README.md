# Figure 5

This folder contains the raw data and analysis script used for the Fig. 5 seed-normalization analysis. The count table is stored as a CSV and is read by the analysis script so the data are not duplicated in code.

## Files

- `data/fig5_normal_seed.csv`: raw counts used to quantify normal seed formation across genotypes.
- `scripts/fig5_normal_seed.R`: plotting and summary script for genotype comparisons and significance annotations.

## Notes

The dataset in `fig5_normal_seed.csv` is the source-of-truth data used by the analysis script. This archive keeps the figure-specific raw values and analysis logic in one place for transparent reproducibility.
