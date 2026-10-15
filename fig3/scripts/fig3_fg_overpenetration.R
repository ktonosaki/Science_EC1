library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

# Reconstructed data from the published figure panels f and g.
# Values are approximated from the plotted points in the supplied image and
# are intended for reproducible plotting and downstream analysis.

df <- read.csv("../data/fig3_overpenetration.csv", stringsAsFactors = FALSE)

df_long <- df %>%
  pivot_longer(
    cols = c(middle_rSY_pct, proximal_to_EC_pct, central_cell_pct),
    names_to = "trait",
    values_to = "pct"
  ) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "ec1")),
    trait = factor(
      trait,
      levels = c("middle_rSY_pct", "proximal_to_EC_pct", "central_cell_pct"),
      labels = c("Middle of rSY", "Proximal to EC", "In central cell")
    )
  )

# Summary stats for the panel-style representation
summary_df <- df_long %>%
  group_by(trait, genotype) %>%
  summarise(
    q1 = quantile(pct, 0.25),
    med = median(pct),
    q3 = quantile(pct, 0.75),
    .groups = "drop"
  )

# Pairwise Welch t-tests for each trait
pvals <- df_long %>%
  group_by(trait) %>%
  summarise(
    p_value = t.test(pct ~ genotype, var.equal = FALSE)$p.value,
    .groups = "drop"
  )

pval_text <- pvals %>%
  mutate(
    label = case_when(
      p_value < 1e-4 ~ "P < 1e-4",
      TRUE ~ sprintf("P = %.2e", p_value)
    )
  )

# Plot: panel f + g in a combined format
p <- ggplot(df_long, aes(x = genotype, y = pct, color = genotype)) +
  geom_beeswarm(size = 2.2, alpha = 0.85, priority = "density") +
  geom_linerange(
    data = summary_df,
    aes(x = genotype, ymin = q1, ymax = q3),
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype, ymin = med, ymax = med),
    width = 0.22,
    linewidth = 0.9,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype, ymin = q1, ymax = q1),
    width = 0.10,
    linewidth = 0.7,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = summary_df,
    aes(x = genotype, ymin = q3, ymax = q3),
    width = 0.10,
    linewidth = 0.7,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ trait, nrow = 1) +
  scale_color_manual(values = c("WT" = "grey40", "ec1" = "#D95F02")) +
  scale_y_continuous(name = "Percentage (%)", limits = c(0, 40)) +
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
ggsave("../plots/fig3_fg_overpenetration.pdf", p, width = 7, height = 3.2)

central_df <- df_long %>% filter(trait == "In central cell")
central_summary <- summary_df %>% filter(trait == "In central cell")

p2 <- ggplot(central_df, aes(x = genotype, y = pct, color = genotype)) +
  geom_beeswarm(size = 2.3, alpha = 0.85, priority = "density") +
  geom_linerange(
    data = central_summary,
    aes(x = genotype, ymin = q1, ymax = q3),
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = central_summary,
    aes(x = genotype, ymin = med, ymax = med),
    width = 0.22,
    linewidth = 0.9,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("WT" = "grey40", "ec1" = "#D95F02")) +
  geom_text(
    data = pval_text %>% filter(trait == "In central cell"),
    mapping = aes(x = 1.5, y = 35, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_y_continuous(name = "Over-penetration (%)", limits = c(0, 40)) +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 11)
  )

print(p2)
ggsave("../plots/fig3_fg_central_cell_only.pdf", p2, width = 3.5, height = 3)

# Save summary stats and p-values as a report
write.csv(pvals, "../results/fig3_fg_pvalues.csv", row.names = FALSE)
write.csv(summary_df, "../results/fig3_fg_summary_stats.csv", row.names = FALSE)
