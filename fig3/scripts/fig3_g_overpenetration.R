library(dplyr)
library(ggplot2)
library(ggbeeswarm)

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
}
script_path <- script_dir()
setwd(script_path)

# ---------------------------------------------------------------------------
# Fig. 3g reconstructed over-penetration panel.
# The original per-pistil central-cell percentages are not available in the
# main count table, so this panel is generated from reconstructed percentages.
# ---------------------------------------------------------------------------

df <- read.csv(
  file.path(script_path, "..", "data", "fig3_overpenetration_reconstructed.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "ec1")),
    central_cell_pct = as.numeric(central_cell_pct)
  )

summary_df <- df %>%
  group_by(genotype) %>%
  summarise(
    q1 = quantile(central_cell_pct, 0.25),
    med = median(central_cell_pct),
    q3 = quantile(central_cell_pct, 0.75),
    .groups = "drop"
  )

p_value <- t.test(central_cell_pct ~ genotype, data = df, var.equal = FALSE)$p.value
pval_df <- tibble(
  trait = "In central cell",
  p_value = p_value,
  label = sprintf("P = %.2e", p_value)
)

p <- ggplot(df, aes(x = genotype, y = central_cell_pct, color = genotype)) +
  geom_beeswarm(size = 2.2, alpha = 0.85, priority = "density") +
  geom_linerange(
    data = summary_df,
    aes(x = genotype, ymin = q1, ymax = q3),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype, ymin = med, ymax = med),
    width = 0.22,
    linewidth = 0.9,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_text(
    data = pval_df,
    aes(x = 1.5, y = max(df$central_cell_pct) + 2, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_color_manual(values = c("WT" = "grey40", "ec1" = "#D95F02")) +
  scale_y_continuous(name = "Over-penetration (%)", limits = c(0, 40)) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11)
  )

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

ggsave("../plots/fig3_g_overpenetration_summary.pdf", p, width = 3.5, height = 3)
write.csv(summary_df, "../results/fig3_g_overpenetration_summary.csv", row.names = FALSE)
write.csv(pval_df, "../results/fig3_g_overpenetration_pvalues.csv", row.names = FALSE)
