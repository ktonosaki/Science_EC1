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
  file.path(script_path, "..", "data", "fig5e_abnormal_ovule_counts.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "ec1", "NM_1", "NM_2", "NM_3", "M_1")),
    total = normal + abnormal,
    rate = 100 * abnormal / total,
    genotype2 = ifelse(genotype %in% c("NM_1", "NM_2", "NM_3"), "NM", as.character(genotype)),
    genotype2 = factor(genotype2, levels = c("WT", "ec1", "NM", "M_1"))
  )

sum_df <- df %>%
  group_by(genotype2) %>%
  summarise(
    q1 = quantile(rate, 0.25),
    med = median(rate),
    q3 = quantile(rate, 0.75),
    .groups = "drop"
  )

df_cmp <- df %>%
  filter(genotype2 %in% c("ec1", "NM", "M_1")) %>%
  mutate(genotype2 = droplevels(genotype2))

fit <- glm(
  cbind(normal, abnormal) ~ genotype2,
  family = binomial(link = "logit"),
  data = df_cmp
)

coef_vec <- coef(fit)
vcov_mat <- vcov(fit)
contrast_mat <- rbind(
  `ec1 - NM` = c(0, -1, 0),
  `ec1 - M_1` = c(0, 0, -1),
  `NM - M_1` = c(0, 1, -1)
)

estimate <- as.vector(contrast_mat %*% coef_vec)
std_error <- sqrt(diag(contrast_mat %*% vcov_mat %*% t(contrast_mat)))
z_value <- estimate / std_error
p_value_raw <- 2 * pnorm(abs(z_value), lower.tail = FALSE)

logit_results <- tibble(
  contrast = rownames(contrast_mat),
  estimate = estimate,
  std_error = std_error,
  z_value = z_value,
  p_value_raw = p_value_raw
) %>%
  mutate(
    p_value = p.adjust(p_value_raw, method = "BH"),
    label = sprintf("P = %.2e", p_value)
  )

ymax <- max(df$rate)
anno_df <- tibble(
  x1 = c(2, 2, 3),
  x2 = c(4, 3, 4),
  y = c(ymax + 11, ymax + 3, ymax + 7),
  label = c(
    logit_results$label[logit_results$contrast == "ec1 - M_1"],
    logit_results$label[logit_results$contrast == "ec1 - NM"],
    logit_results$label[logit_results$contrast == "NM - M_1"]
  )
)

p <- ggplot(df, aes(x = genotype2, y = rate, color = genotype2)) +
  geom_quasirandom(groupOnX = TRUE, width = 0.2, size = 1, alpha = 0.85) +
  geom_linerange(
    data = sum_df,
    aes(x = genotype2, ymin = q1, ymax = q3),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_errorbar(
    data = sum_df,
    aes(x = genotype2, ymin = med, ymax = med),
    width = 0.24,
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = sum_df,
    aes(x = genotype2, ymin = q1, ymax = q1),
    width = 0.10,
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = sum_df,
    aes(x = genotype2, ymin = q3, ymax = q3),
    width = 0.10,
    linewidth = 0.6,
    inherit.aes = FALSE
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
  scale_y_continuous(
    name = "Abnormal ovules (%)",
    limits = c(0, ymax + 15),
    expand = expansion(mult = c(0, 0.02))
  ) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 5)
  )

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

ggsave("../plots/fig5e_abnormal_ovule.pdf", p, width = 5, height = 3)
write.csv(sum_df, "../results/fig5e_abnormal_ovule_summary.csv", row.names = FALSE)
write.csv(logit_results, "../results/fig5e_abnormal_ovule_pvalues.csv", row.names = FALSE)
