source("config.R")

df <- fread('../data/gtex_v8/spliceai_window_size.csv')

chosen_models <- c(
  "SpliceAI (50bp)",
  "SpliceAI (100bp)", 
  "SpliceAI (150bp)",
  "SpliceAI (200bp)",
  "SpliceAI (250bp)"
)

df <- process_models_pr_curve(df, chosen_models)
df$model <- factor(df$model, levels = chosen_models)

SI_fig_1 <- (
  ggplot(df, aes(x=recall, y=precision, color=model))
  + geom_step(direction = "vh")
  + scale_colour_brewer(palette = "Set1")
  + theme_cowplot(font_size = fontsize)
  # + background_grid(major = "xy", minor = "xy")
  + scale_x_continuous(limits = c(-0.01, 0.6), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2))
  + scale_y_continuous(limits = c(-0.01, 0.2), breaks = seq(0, 1, 0.1), minor_breaks = seq(0, 1, 0.02))
  + labs(x='Recall', y='Precision', color='')
)
SI_fig_1

ggsave(plot=SI_fig_1, filename="out/Supplementary_Figures/Supplementary_Figure_1/Supplementary_Figure_1.pdf", height=unit(7, "cm"), width=unit(12, "cm"), dpi=450)
ggsave(plot=SI_fig_1, filename="out/Supplementary_Figures/Supplementary_Figure_1/Supplementary_Figure_1.png", height=unit(7, "cm"), width=unit(12, "cm"), dpi=450)
