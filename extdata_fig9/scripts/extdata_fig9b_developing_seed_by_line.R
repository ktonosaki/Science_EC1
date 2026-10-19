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

df <- read.csv(file.path(script_path, "..", "data", "extdata_fig9b_developing_seed_by_line.csv"), stringsAsFactors = FALSE) %>%
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

contrast_pairs <- list(
  c("ec1", "NM_1"),
  c("ec1", "NM_2"),
  c("ec1", "NM_3")
)

fit <- glm(
  cbind(developing, total - developing) ~ genotype2,
  family = binomial(link = "logit"),
  data = df
)

coef_vec <- coef(fit)
vcov_mat <- vcov(fit)
contrast_mat <- rbind(
  `ec1 - NM_1` = c(0, -1, 0, 0),
  `ec1 - NM_2` = c(0, 0, -1, 0),
  `ec1 - NM_3` = c(0, 0, 0, -1),
  `NM_1 - NM_2` = c(0, 1, -1, 0),
  `NM_1 - NM_3` = c(0, 1, 0, -1),
  `NM_2 - NM_3` = c(0, 0, 1, -1)
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
    label = ifelse(p_value >= 0.01, sprintf("P = %.3f", p_value), sprintf("P = %.2e", p_value))
  )

max_rate <- max(df$rate)
anno_df <- tibble(
  x1 = c(1, 1, 1),
  x2 = c(2, 3, 4),
  y = c(79, 89, 102),
  label = c(
    logit_results$label[logit_results$contrast == "ec1 - NM_1"],
    logit_results$label[logit_results$contrast == "ec1 - NM_2"],
    logit_results$label[logit_results$contrast == "ec1 - NM_3"]
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
  scale_color_manual(values = c("ec1" = "#F8766D", "NM_1" = "#7CAE00", "NM_2" = "#00BFC4", "NM_3" = "#C77CFF")) +
  scale_x_discrete(labels = c(
    "ec1" = "ec1\nnull\nn = 11",
    "NM_1" = "line 1\nn = 13",
    "NM_2" = "line 2\nn = 15",
    "NM_3" = "line 3\nn = 21"
  )) +
  scale_y_continuous(name = "Developing seeds (%)", limits = c(0, 105)) +
  xlab("Sup35NM repl") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11),
    axis.title.x = element_text(face = "italic", size = 11),
    axis.text.x = element_text(size = 9)
  )

if (!dir.exists("../plots")) dir.create("../plots", recursive = TRUE)
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)

ggsave("../plots/extdata_fig9b_developing_seed_by_line.pdf", p, width = 5.5, height = 3.5)
write.csv(summary_df, "../results/extdata_fig9b_developing_seed_by_line_summary.csv", row.names = FALSE)
write.csv(logit_results, "../results/extdata_fig9b_developing_seed_by_line_pvalues.csv", row.names = FALSE)
