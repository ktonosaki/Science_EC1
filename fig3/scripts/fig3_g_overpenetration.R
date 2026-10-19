library(dplyr)
library(tidyr)
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
# Fig. 3g over-penetration analysis based on per-pistil count data.
# P-values are taken from binomial logistic regression likelihood-ratio tests.
# ---------------------------------------------------------------------------

df <- read.csv(
  file.path(script_path, "..", "data", "fig3_g_overpenetration_counts.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "ec1")),
    total = as.integer(total),
    abnormal = as.integer(abnormal)
  )

df_long <- df %>%
  pivot_longer(
    cols = c(abnormal),
    names_to = "pheno",
    values_to = "abnormal"
  ) %>%
  mutate(
    rate = 100 * abnormal / total,
    pheno = factor(pheno, levels = "abnormal", labels = "central cell")
  )

summary_df <- df_long %>%
  group_by(pheno, genotype) %>%
  summarise(
    q1 = quantile(rate, 0.25),
    med = median(rate),
    q3 = quantile(rate, 0.75),
    .groups = "drop"
  )

logit_results <- df_long %>%
  group_by(pheno) %>%
  group_modify(~ {
    fit_null <- glm(
      cbind(abnormal, total - abnormal) ~ 1,
      family = binomial(link = "logit"),
      data = .x
    )
    fit_full <- glm(
      cbind(abnormal, total - abnormal) ~ genotype,
      family = binomial(link = "logit"),
      data = .x
    )
    lrt <- anova(fit_null, fit_full, test = "LRT")

    tibble(
      chi_square = lrt$Deviance[2],
      df = lrt$Df[2],
      p_value = lrt$`Pr(>Chi)`[2]
    )
  }) %>%
  ungroup() %>%
  mutate(label = sprintf("P = %.2e", p_value))

anno_df <- tibble(
  pheno = levels(df_long$pheno),
  x1 = 1,
  x2 = 2,
  y = max(df_long$rate) + 5,
  label = logit_results$label
)

p <- ggplot(df_long, aes(x = genotype, y = rate, color = genotype)) +
  geom_quasirandom(groupOnX = TRUE, width = 0.2, size = 2, alpha = 0.85) +
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
  geom_errorbar(
    data = summary_df,
    aes(x = genotype, ymin = q1, ymax = q1),
    width = 0.10,
    linewidth = 0.6,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype, ymin = q3, ymax = q3),
    width = 0.10,
    linewidth = 0.6,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x1, xend = x1, y = y, yend = y - 1),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x2, xend = x2, y = y, yend = y - 1),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_text(
    data = anno_df,
    aes(x = 1.5, y = y + 1.5, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ pheno, nrow = 1, scales = "fixed") +
  scale_color_manual(values = c("WT" = "grey40", "ec1" = "#D95F02")) +
  scale_y_continuous(name = "Abnormal ovules (%)", limits = c(0, 30)) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "none"
  )

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

ggsave("../plots/fig3_g_overpenetration_summary.pdf", p, width = 3, height = 3)
write.csv(summary_df, "../results/fig3_g_overpenetration_summary.csv", row.names = FALSE)
write.csv(logit_results, "../results/fig3_g_overpenetration_pvalues.csv", row.names = FALSE)
