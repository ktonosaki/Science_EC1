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
  file.path(script_path, "..", "data", "extdata_fig9_abnormal_by_line_counts.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    genotype = factor(genotype, levels = c("ec1", "NM_1", "NM_2", "NM_3")),
    normal = as.integer(normal),
    abnormal = as.integer(abnormal),
    total = normal + abnormal,
    rate = 100 * abnormal / total
  )

summary_df <- df %>%
  group_by(genotype) %>%
  summarise(
    q1 = quantile(rate, 0.25),
    med = median(rate),
    q3 = quantile(rate, 0.75),
    .groups = "drop"
  )

contrast_pairs <- list(
  c("ec1", "NM_1"),
  c("ec1", "NM_2"),
  c("ec1", "NM_3")
)

pvals <- lapply(contrast_pairs, function(pair) {
  g1 <- pair[1]
  g2 <- pair[2]
  sub <- dplyr::filter(df, genotype %in% c(g1, g2)) %>%
    mutate(genotype = factor(genotype, levels = c(g1, g2)))

  fit_null <- glm(
    cbind(abnormal, normal) ~ 1,
    family = binomial(link = "logit"),
    data = sub
  )
  fit_full <- glm(
    cbind(abnormal, normal) ~ genotype,
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
  x1 = c(1, 1, 1),
  x2 = c(2, 3, 4),
  y = c(
    min(26, max_rate * 0.78),
    min(28, max_rate * 0.86),
    min(30, max_rate * 0.94)
  ),
  label = c(
    pvals$label[pvals$contrast == "ec1 - NM_1"],
    pvals$label[pvals$contrast == "ec1 - NM_2"],
    pvals$label[pvals$contrast == "ec1 - NM_3"]
  )
)

p <- ggplot(df, aes(x = genotype, y = rate, color = genotype)) +
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
    width = 0.24,
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
  scale_color_manual(values = c("ec1" = "#D95F02", "NM_1" = "#7570B3", "NM_2" = "#E7298A", "NM_3" = "#66A61E")) +
  scale_y_continuous(name = "Abnormal ovules (%)", limits = c(0, 32)) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11)
  )

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

ggsave("../plots/extdata_fig9_abnormal_by_line_summary.pdf", p, width = 5.5, height = 3.5)
write.csv(summary_df, "../results/extdata_fig9_abnormal_by_line_summary.csv", row.names = FALSE)
write.csv(pvals, "../results/extdata_fig9_abnormal_by_line_pvalues.csv", row.names = FALSE)
