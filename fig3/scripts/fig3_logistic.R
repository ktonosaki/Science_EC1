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
# Fig. 3f abnormal ovule phenotype analysis.
# Keep the raw counts in CSV form and load them here so the analysis script does
# not duplicate the tabulated data in code.
# ---------------------------------------------------------------------------

df <- read.csv(file.path(script_path, "..", "data", "fig3_abnormal_ovule.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "ec1")),
    total = as.integer(total),
    normal = as.integer(normal),
    abnormal_1 = as.integer(abnormal_1),
    abnormal_2 = as.integer(abnormal_2)
  )

# Convert the count table to a long data frame for plotting and summary stats.
df_long <- df %>%
  pivot_longer(
    cols = c(abnormal_1, abnormal_2),
    names_to = "pheno",
    values_to = "abnormal"
  ) %>%
  mutate(
    rate = 100 * abnormal / total,
    genotype = factor(genotype, levels = c("WT", "ec1")),
    pheno = factor(
      pheno,
      levels = c("abnormal_2", "abnormal_1"),
      labels = c("Middle of rSY", "Proximal to EC")
    )
  )

# Summary statistics used for the confidence-interval-style overlay.
sum_df <- df_long %>%
  group_by(pheno, genotype) %>%
  summarise(
    q1  = quantile(rate, 0.25),
    med = median(rate),
    q3  = quantile(rate, 0.75),
    .groups = "drop"
  )

# Keep the figure annotations aligned with the paper's statistical validation:
# the main annotation labels are based on the binomial logistic regression model.
# Fisher exact tests are retained only as a complementary sanity check.
mat1 <- matrix(c(
  sum(df$abnormal_1[df$genotype == "ec1"]),
  sum(df$total[df$genotype == "ec1"] - df$abnormal_1[df$genotype == "ec1"]),
  sum(df$abnormal_1[df$genotype == "WT"]),
  sum(df$total[df$genotype == "WT"] - df$abnormal_1[df$genotype == "WT"])
), nrow = 2, byrow = TRUE)

p1 <- fisher.test(mat1)$p.value

mat2 <- matrix(c(
  sum(df$abnormal_2[df$genotype == "ec1"]),
  sum(df$total[df$genotype == "ec1"] - df$abnormal_2[df$genotype == "ec1"]),
  sum(df$abnormal_2[df$genotype == "WT"]),
  sum(df$total[df$genotype == "WT"] - df$abnormal_2[df$genotype == "WT"])
), nrow = 2, byrow = TRUE)

p2 <- fisher.test(mat2)$p.value

fisher_results <- tibble(
 pheno = levels(df_long$pheno),
 fisher_p_value = c(p1, p2)
)

# Show the exact p-value in the publication-style annotation instead of asterisks.
p_to_label <- function(p) sprintf("P = %.2e", p)

# Logistic regression is the validation used for the figure annotation in the paper.
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
  mutate(label = sapply(p_value, p_to_label))

# Use the model-based p-values for the annotation labels shown in the figure.
anno_df <- data.frame(
  pheno = levels(df_long$pheno),
  x1 = 1,
  x2 = 2,
  y = c(
    max(df_long$rate[df_long$pheno == levels(df_long$pheno)[1]]) + 5,
    max(df_long$rate[df_long$pheno == levels(df_long$pheno)[2]]) + 5
  ),
  label = logit_results$label
)

# Create the publication-style boxplot with beeswarm points and significance bars.
p <- ggplot(df_long, aes(x = genotype, y = rate, color = genotype)) +
  geom_beeswarm(
    size = 2.2,
    alpha = 0.85,
    priority = "density"
  ) +
  geom_linerange(
    data = sum_df,
    aes(x = genotype, ymin = q1, ymax = q3),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_errorbar(
    data = sum_df,
    aes(x = genotype, ymin = med, ymax = med),
    width = 0.22,
    linewidth = 0.9,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbar(
    data = sum_df,
    aes(x = genotype, ymin = q1, ymax = q1),
    width = 0.10,
    linewidth = 0.6,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbar(
    data = sum_df,
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
  scale_y_continuous(name = "Abnormal ovules (%)", limits = c(0, 40)) +
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
ggsave("../plots/fig3_abnormal_ovule_summary.pdf", p, width = 3.5, height = 3)

# ---------------------------------------------------------------------------
# Binomial logistic regression with a likelihood-ratio test.
# These are the validation results used to match the paper's annotation logic.
# ---------------------------------------------------------------------------
logit_results <- logit_results %>%
  arrange(factor(pheno, levels = levels(df_long$pheno)))

print(logit_results)

if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

pval_df <- logit_results %>%
  left_join(fisher_results, by = "pheno") %>%
  select(pheno, p_value, fisher_p_value, label)

write.csv(sum_df, "../results/fig3_abnormal_ovule_summary.csv", row.names = FALSE)
write.csv(pval_df, "../results/fig3_abnormal_ovule_pvalues.csv", row.names = FALSE)
