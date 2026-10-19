# Science_EC1 figure archive

This repository stores raw data and analysis scripts by figure.

## Structure

```text
fig3/
  README.md
  data/
    fig3_abnormal_ovule.csv
    fig3_g_overpenetration_counts.csv
  scripts/
    fig3_logistic.R
    fig3_g_overpenetration.R
  plots/
    fig3_abnormal_ovule_summary.pdf
    fig3_g_overpenetration_summary.pdf
  results/
    fig3_abnormal_ovule_pvalues.csv
    fig3_abnormal_ovule_summary.csv
    fig3_g_overpenetration_pvalues.csv
    fig3_g_overpenetration_summary.csv

fig5/
  README.md
  data/
    fig5_abnormal_ovule_counts.csv
    fig5_normal_seed.csv
  scripts/
    fig5_e_abnormal_ovule.R
    fig5_normal_seed.R
  plots/
    fig5_e_abnormal_ovule_plot.pdf
    fig5_normal_seed_summary.pdf
  results/
    fig5_e_abnormal_ovule_pvalues.csv
    fig5_e_abnormal_ovule_summary.csv
    fig5_normal_seed_pvalues.csv
    fig5_normal_seed_summary.csv

extdata_fig9/
  README.md
  data/
    extdata_fig9_seed_by_line.csv
  scripts/
    extdata_fig9_seed_by_line.R
  plots/
    extdata_fig9_seed_by_line_summary.pdf
  results/
    extdata_fig9_seed_by_line_pvalues.csv
    extdata_fig9_seed_by_line_summary.csv
```

Each figure gets its own directory so data, scripts, and figure-specific notes remain isolated and easy to track.
