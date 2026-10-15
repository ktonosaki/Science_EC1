# Figure 3

This folder contains the raw data and analysis scripts used for Figure 3, including the logistic-regression summary and the reconstructed panels f/g.

## Files

- `data/fig3_abnormal_ovule.csv`: counts for abnormal ovule phenotypes reported in the figure analysis.
- `data/fig3_overpenetration.csv`: reconstructed values for the figure f/g over-penetration summary.
- `scripts/fig3_logistic.R`: logistic model script used for the main Figure 3 analysis.
- `scripts/fig3_fg_overpenetration.R`: script used to generate the panel-style figure f/g summary.

## Data provenance

The observation counts in `fig3_abnormal_ovule.csv` are derived from the values embedded in the original analysis script. The over-penetration values in `fig3_overpenetration.csv` were reconstructed from the plotted points in the publication figure image, so they should be treated as a transparent, reproducible approximation until the original raw measurements are available.
