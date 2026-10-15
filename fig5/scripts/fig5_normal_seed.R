library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

# Read the raw Fig.5 seed counts
# The original script was written directly as a pasted table; we preserve the values in CSV format.
df <- read.csv("../data/fig5_normal_seed.csv", stringsAsFactors = FALSE)

# Compute total and treatment rate
# The original script groups 'normal + aborted' as the successful developing fraction,
# while undeveloped is treated as the non-developing fraction.
df <- df %>%
  mutate(
    total = normal + aborted + undeveloped,
    developing = normal + aborted,
    rate = 100 * developing / total
  )

# Standardize genotype labels used in the original analysis
# WT, ec1, NM, M and also individual NM_1/2/3 variants are available in the raw data.
df <- df %>%
  mutate(
    genotype2 = case_when(
      genotype == "WT" ~ "WT",
      genotype == "ec1_null" ~ "ec1",
      grepl("SUP35NM", genotype) ~ "NM",
      genotype == "ec1_null+SUP35M" ~ "M",
      TRUE ~ genotype
    )
  ) %>%
  mutate(genotype2 = factor(genotype2, levels = c("WT", "ec1", "NM", "M")))

# Summary stats for boxplot-style whiskers
summary_df <- df %>%
  group_by(genotype2) %>%
  summarise(
    q1 = quantile(rate, 0.25),
    med = median(rate),
    q3 = quantile(rate, 0.75),
    .groups = "drop"
  )

# Pairwise comparisons using Fisher exact tests against the original script logic.
# Compare the fraction of developing seeds (normal + aborted) vs undeveloped seeds.
# This mirrors the original pattern of tests for genotype effects.
contrast_pairs <- list(
  c("ec1", "NM"),
  c("ec1", "M"),
  c("NM", "M")
)

pvals <- lapply(contrast_pairs, function(pair) {
  g1 <- pair[1]
  g2 <- pair[2]
  sub <- dplyr::filter(df, genotype2 %in% c(g1, g2))
  mat <- matrix(c(
    sum(sub$developing[sub$genotype2 == g1]),
    sum(sub$total[sub$genotype2 == g1] - sub$developing[sub$genotype2 == g1]),
    sum(sub$developing[sub$genotype2 == g2]),
    sum(sub$total[sub$genotype2 == g2] - sub$developing[sub$genotype2 == g2])
  ), nrow = 2, byrow = TRUE)
  tibble(
    contrast = paste(g1, g2, sep = " - "),
    p_value = fisher.test(mat)$p.value
  )
}) %>%
  bind_rows()

p_to_symbol <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

pvals <- pvals %>%
  mutate(label = sapply(p_value, p_to_symbol))

# Annotation positions matching the original publication-style plots.
# Because the plot contains three genotype groups only, positions must stay within
# x = 1, 2, 3 to avoid warnings from out-of-range segments.
max_rate <- max(df$rate)
anno_df <- tibble(
  x1 = c(1, 1, 2),
  x2 = c(2, 3, 3),
  y = c(max_rate + 3, max_rate + 8, max_rate + 13),
  label = c(
    pvals$label[pvals$contrast == "ec1 - NM"],
    pvals$label[pvals$contrast == "ec1 - M"],
    pvals$label[pvals$contrast == "NM - M"]
  )
)

p <- ggplot(df, aes(x = genotype2, y = rate, color = genotype2)) +
  geom_quasirandom(width = 0.2, size = 2, alpha = 0.85, groupOnX = TRUE) +
  geom_linerange(
    data = summary_df,
    aes(x = genotype2, ymin = q1, ymax = q3),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype2, ymin = med, ymax = med),
    width = 0.24,
    linewidth = 0.9,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype2, ymin = q1, ymax = q1),
    width = 0.10,
    linewidth = 0.6,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype2, ymin = q3, ymax = q3),
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
    aes(x = x1, xend = x1, y = y, yend = y - 2),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x2, xend = x2, y = y, yend = y - 2),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_text(
    data = anno_df,
    aes(x = 1.5, y = y + 2, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_color_manual(values = c("WT" = "grey40", "ec1" = "#D95F02", "NM" = "#7570B3", "M" = "#1B9E77")) +
  scale_y_continuous(name = "Normal seeds (%)", limits = c(0, 100)) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11)
  )

print(p)

ggsave("../plots/fig5_normal_seed_summary.pdf", p, width = 5.5, height = 3.5)

# Output a concise result table
write.csv(summary_df, "../results/fig5_normal_seed_summary.csv", row.names = FALSE)
write.csv(pvals, "../results/fig5_normal_seed_pvalues.csv", row.names = FALSE)
