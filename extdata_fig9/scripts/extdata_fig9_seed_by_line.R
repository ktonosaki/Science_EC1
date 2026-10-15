library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# Extended Data Fig. 9: independent Sup35NM replacement lines
# The raw values in this archive are reconstructed from the original analysis script.

df <- read.csv("../data/extdata_fig9_seed_by_line.csv", stringsAsFactors = FALSE)

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
}) %>% bind_rows()

p_to_symbol <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

pvals <- pvals %>% mutate(label = sapply(p_value, p_to_symbol))

max_rate <- max(df$rate)
anno_df <- tibble(
  x1 = c(1, 1, 1),
  x2 = c(2, 3, 4),
  y = c(max_rate + 3, max_rate + 8, max_rate + 13),
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

write.csv(summary_df, "../results/extdata_fig9_seed_by_line_summary.csv", row.names = FALSE)
write.csv(pvals, "../results/extdata_fig9_seed_by_line_pvalues.csv", row.names = FALSE)
