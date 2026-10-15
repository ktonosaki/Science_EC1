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

# Extended Data Fig. 9: independent Sup35NM replacement lines
# The raw values in this archive are reconstructed from the original analysis script.

df <- read.csv(file.path(script_path, "..", "data", "extdata_fig9_seed_by_line.csv"), stringsAsFactors = FALSE)

df <- df %>%
  mutate(
    total = normal + aborted + undeveloped,
    developing = normal + aborted,
    rate = 100 * developing / total,
    genotype2 = case_when(
      genotype == "ec1_null" ~ "ec1",
      genotype == "ec1_null+SUP35NM1" ~ "NM_1",
      genotype == "ec1_null+SUP35NM2" ~ "NM_2",
      genotype == "ec1_null+SUP35NM3" ~ "NM_3",
      TRUE ~ genotype
    )
  ) %>%
  filter(genotype2 %in% c("ec1", "NM_1", "NM_2", "NM_3")) %>%
  mutate(genotype2 = factor(genotype2, levels = c("ec1", "NM_1", "NM_2", "NM_3")))

summary_df <- df %>%
  group_by(genotype2) %>%
  summarise(
    q1 = quantile(rate, 0.25),
    med = median(rate),
    q3 = quantile(rate, 0.75),
    .groups = "drop"
  )

contrasts <- list(
  c("ec1", "NM_1"),
  c("ec1", "NM_2"),
  c("ec1", "NM_3")
)

pvals <- lapply(contrasts, function(pair) {
  g1 <- pair[1]
  g2 <- pair[2]
  sub <- dplyr::filter(df, genotype2 %in% c(g1, g2)) %>%
    mutate(genotype2 = factor(genotype2, levels = c(g1, g2)))

  fit_null <- glm(
    cbind(developing, total - developing) ~ 1,
    family = binomial(link = "logit"),
    data = sub
  )
  fit_full <- glm(
    cbind(developing, total - developing) ~ genotype2,
    family = binomial(link = "logit"),
    data = sub
  )
  lrt <- anova(fit_null, fit_full, test = "LRT")

  tibble(
    contrast = paste(g1, g2, sep = " - "),
    p_value = lrt$`Pr(>Chi)`[2]
  )
}) %>% bind_rows()

fishers <- lapply(contrasts, function(pair) {
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
    fisher_p_value = fisher.test(mat)$p.value
  )
}) %>% bind_rows()

p_to_label <- function(p) sprintf("P = %.2e", p)

pvals <- pvals %>%
  left_join(fishers, by = "contrast") %>%
  mutate(label = sapply(p_value, p_to_label))

max_rate <- max(df$rate)
anno_df <- tibble(
  x1 = c(1, 1, 1),
  x2 = c(2, 3, 4),
  y = c(
    min(75, max_rate * 0.82),
    min(85, max_rate * 0.90),
    min(95, max_rate * 0.97)
  ),
  label = c(
    pvals$label[pvals$contrast == "ec1 - NM_1"],
    pvals$label[pvals$contrast == "ec1 - NM_2"],
    pvals$label[pvals$contrast == "ec1 - NM_3"]
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
    aes(x = 2.5, y = y + 2, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_color_manual(values = c("ec1" = "#D95F02", "NM_1" = "#7570B3", "NM_2" = "#E7298A", "NM_3" = "#66A61E")) +
  scale_y_continuous(name = "Normal seed (%)", limits = c(0, 100)) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11)
  )

print(p)

ggsave("../plots/extdata_fig9_seed_by_line_summary.pdf", p, width = 5.5, height = 3.5)

# The annotation values are drawn from the logistic regression validation, while the
# Fisher exact values are included for cross-checking with the original analysis.
write.csv(summary_df, "../results/extdata_fig9_seed_by_line_summary.csv", row.names = FALSE)
write.csv(pvals, "../results/extdata_fig9_seed_by_line_pvalues.csv", row.names = FALSE)
