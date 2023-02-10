source("config.R")

chosen_models <- c(
  'AbSplice-DNA',
  'SpliceAI + SpliceMap + Ψ_ref',
  'SpliceAI + SpliceMap',
  'SpliceAI',
  'MMSplice + SpliceMap + Ψ_ref',
  'MMSplice + SpliceMap',
  'MMSplice'
)
color_list = absplice_model_colors[chosen_models]

# SI Fig. 7a
df <- fread('../data/LOEUF/matched_recall.csv')

df$model <- factor(df$model, levels=chosen_models)

# df$pval_rounded <- round(df$pval_spliceai, 2)
# df$pval_rounded_label <- ifelse(((df$pval_rounded <= 0.05) && df$significance_label_spliceai != ''), paste0('P=', df$pval_rounded), 'N.S.')

# df$pval_spliceai_label <- factor(df$pval_spliceai_label)

df$pval_spliceai_label <- ifelse(df$pval_spliceai_label <= 0.05, paste0('P=', df$pval_spliceai_label), 'N.S.')

g_LOEUF_matched_recall <- (
  ggplot(df)
  + aes_string(
   x='LOEUF_decile', 
   y='depletion_high_impact_variants',
   fill='model', 
   # label="significance_label_spliceai"
   label="pval_spliceai_label"
  )
  + geom_text(
   # nudge_x = -0.2,
   nudge_y = 0.3,
   size=4, 
   color="black")
  + geom_bar(stat='identity', position='dodge')
  + labs(x='LOEUF decile', y='Depletion \nof high impact SNVs')
  + scale_fill_manual(values=color_list)
  + theme_cowplot(font_size=fontsize)
  # + background_grid()
  + scale_x_continuous(
   breaks = seq(0, 10, 1),
  )
  + theme(
   legend.title = element_blank(),
   legend.position= c(0.7, 0.8)
  )
)
g_LOEUF_matched_recall
ggsave(plot=g_LOEUF_matched_recall, filename="out/Supplementary_Figures/Supplementary_Figure_7/Supplementary_Figure_7a.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)

# SI Fig. 7b
df <- fread('../data/LOEUF/matched_top_n.csv')
df$model <- factor(df$model, levels=chosen_models)

df$pval_rounded_label <- ifelse(df$pval_rounded <= 0.05, paste0('P=', df$pval_rounded), 'N.S.')

g_LOEUF_matched_top_n <- (
  ggplot(df)
  + aes_string(
    x='LOEUF_decile', 
    y='depletion_high_impact_variants',
    fill='model', 
    # label="significance_label_spliceai"
    label="pval_spliceai_label"
  )
  + geom_text(
    # nudge_x = -0.2,
    nudge_y = 0.3,
    size=4, 
    color="black")
  + geom_bar(stat='identity', position='dodge')
  + labs(x='LOEUF decile', y='Depletion \nof high impact SNVs')
  + scale_fill_manual(values=color_list)
  + theme_cowplot(font_size=fontsize)
  # + background_grid()
  + scale_x_continuous(
    breaks = seq(0, 10, 1),
  )
  + theme(
    legend.title = element_blank(),
    legend.position= c(0.7, 0.8)
  )
)
g_LOEUF_matched_top_n
ggsave(plot=g_LOEUF_matched_top_n, filename="out/Supplementary_Figures/Supplementary_Figure_7/Supplementary_Figure_7b.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)
g_LOEUF_matched_top_n <- g_LOEUF_matched_top_n + theme(legend.position = "none")

SI_fig_7 <- ggarrange(
  nrow=2, widths = c(1, 1), labels = c('a', 'b'), 
  # common.legend = TRUE, legend="right",
  g_LOEUF_matched_recall,
  g_LOEUF_matched_top_n
)
SI_fig_7
ggsave(plot=SI_fig_7, filename="out/Supplementary_Figures/Supplementary_Figure_7/Supplementary_Figure_7_with_pval.svg", height=unit(12, "cm"), width=unit(12, "cm"), dpi=450)
