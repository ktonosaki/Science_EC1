library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

# ---------------------------------------------------------------------------
# Data reconstruction for Fig. 3 abnormal ovule phenotype analysis.
# This script reproduces the original selection of wild-type and ec1 counts
# and generates the boxplot summary used in the publication.
# ---------------------------------------------------------------------------

ec1 <- data.frame(
  genotype = "ec1",
  total = c(72,74,74,66,66,67,63,67,61,68,61,68),
  normal = c(5,11,5,5,11,15,6,4,6,6,10,9),
  abnormal_1 = c(8,4,4,9,11,11,8,9,13,9,13,8),
  abnormal_2 = c(16,16,15,12,18,16,16,15,14,20,14,16)
)

WT <- data.frame(
  genotype = "WT",
  total = c(47,64,64,69,63,67,70,59,73,56,67),
  normal = c(27,27,33,21,36,42,42,32,23,31,38),
  abnormal_1 = c(0,0,0,0,3,0,0,0,0,1,1),
  abnormal_2 = c(0,0,0,0,0,0,0,0,0,0,0)
)

df <- bind_rows(ec1, WT)

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
      labels = c("Near egg cell", "In synergid")
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

print(p)

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
ggsave("../plots/abnormal_ovule_synergid.pdf", p, width = 3.5, height = 3)

# ---------------------------------------------------------------------------
# Binomial logistic regression with a likelihood-ratio test.
# These are the validation results used to match the paper's annotation logic.
# ---------------------------------------------------------------------------
logit_results <- logit_results %>%
  arrange(factor(pheno, levels = levels(df_long$pheno)))

print(logit_results)
