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

df <- read.csv(file.path(script_path, "..", "data", "fig5_normal_seed.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    total = normal + aborted + undeveloped,
    developing = normal + aborted,
    rate = 100 * developing / total,
    genotype2 = case_when(
      genotype == "WT" ~ "WT",
      genotype == "ec1_null" ~ "ec1",
      grepl("SUP35NM", genotype) ~ "NM",
      genotype == "ec1_null+SUP35M" ~ "M",
      TRUE ~ genotype
    ),
    genotype2 = factor(genotype2, levels = c("WT", "ec1", "NM", "M"))
  )

summary_df <- df %>%
  group_by(genotype2) %>%
  summarise(
    q1 = quantile(rate, 0.25),
    med = median(rate),
    q3 = quantile(rate, 0.75),
    .groups = "drop"
  )

contrast_pairs <- list(
  c("ec1", "NM"),
  c("ec1", "M"),
  c("NM", "M")
)

logit_results <- lapply(contrast_pairs, function(pair) {
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
    chi_square = lrt$Deviance[2],
    df = lrt$Df[2],
    p_value_raw = lrt$`Pr(>Chi)`[2]
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_value = p.adjust(p_value_raw, method = "BH"),
    label = sprintf("P = %.2e", p_value)
  )

max_rate <- max(df$rate)
anno_df <- tibble(
  x1 = c(2, 2, 3),
  x2 = c(3, 4, 4),
  y = c(
    min(90, max_rate * 0.85),
    min(94, max_rate * 0.90),
    min(98, max_rate * 0.97)
  ),
  label = c(
    logit_results$label[logit_results$contrast == "ec1 - NM"],
    logit_results$label[logit_results$contrast == "ec1 - M"],
    logit_results$label[logit_results$contrast == "NM - M"]
  )
)

p <- ggplot(df, aes(x = genotype2, y = rate, color = genotype2)) +
  geom_quasirandom(groupOnX = TRUE, width = 0.2, size = 2, alpha = 0.85) +
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
    aes(x = (x1 + x2) / 2, y = y + 2, label = label),
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

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

ggsave("../plots/fig5_normal_seed_summary.pdf", p, width = 5.5, height = 3.5)
write.csv(summary_df, "../results/fig5_normal_seed_summary.csv", row.names = FALSE)
write.csv(logit_results, "../results/fig5_normal_seed_pvalues.csv", row.names = FALSE)
