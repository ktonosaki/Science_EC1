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

df <- read.csv(
  file.path(script_path, "..", "data", "fig5_abnormal_ovule_counts.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "ec1", "NM_1", "NM_2", "NM_3", "M_1")),
    normal = as.integer(normal),
    abnormal = as.integer(abnormal),
    total = normal + abnormal,
    rate = 100 * abnormal / total,
    genotype2 = case_when(
      genotype %in% c("NM_1", "NM_2", "NM_3") ~ "NM",
      TRUE ~ as.character(genotype)
    ),
    genotype2 = factor(genotype2, levels = c("WT", "ec1", "NM", "M_1"))
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
  c("ec1", "M_1"),
  c("NM", "M_1")
)

pvals <- lapply(contrast_pairs, function(pair) {
  g1 <- pair[1]
  g2 <- pair[2]
  sub <- dplyr::filter(df, genotype2 %in% c(g1, g2)) %>%
    mutate(genotype2 = factor(genotype2, levels = c(g1, g2)))

  fit_null <- glm(
    cbind(abnormal, normal) ~ 1,
    family = binomial(link = "logit"),
    data = sub
  )
  fit_full <- glm(
    cbind(abnormal, normal) ~ genotype2,
    family = binomial(link = "logit"),
    data = sub
  )
  lrt <- anova(fit_null, fit_full, test = "LRT")

  tibble(
    contrast = paste(g1, g2, sep = " - "),
    chi_square = lrt$Deviance[2],
    df = lrt$Df[2],
    p_value = lrt$`Pr(>Chi)`[2]
  )
}) %>%
  bind_rows() %>%
  mutate(label = sprintf("P = %.2e", p_value))

max_rate <- max(df$rate)
anno_df <- tibble(
  x1 = c(2, 2, 3),
  x2 = c(3, 4, 4),
  y = c(
    min(26, max_rate * 0.80),
    min(28, max_rate * 0.88),
    min(30, max_rate * 0.96)
  ),
  label = c(
    pvals$label[pvals$contrast == "ec1 - NM"],
    pvals$label[pvals$contrast == "ec1 - M_1"],
    pvals$label[pvals$contrast == "NM - M_1"]
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
    aes(x = x1, xend = x1, y = y, yend = y - 0.5),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = anno_df,
    aes(x = x2, xend = x2, y = y, yend = y - 0.5),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_text(
    data = anno_df,
    aes(x = (x1 + x2) / 2, y = y + 0.8, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_color_manual(values = c("WT" = "grey40", "ec1" = "#D95F02", "NM" = "#7570B3", "M_1" = "#1B9E77")) +
  scale_y_continuous(name = "Abnormal ovules (%)", limits = c(0, 32)) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11)
  )

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

ggsave("../plots/fig5_abnormal_ovule_summary.pdf", p, width = 5.5, height = 3.5)
write.csv(summary_df, "../results/fig5_abnormal_ovule_summary.csv", row.names = FALSE)
write.csv(pvals, "../results/fig5_abnormal_ovule_pvalues.csv", row.names = FALSE)
