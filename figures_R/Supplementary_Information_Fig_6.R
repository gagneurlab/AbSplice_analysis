source("config.R")

chosen_models <- c(
  'AbSplice-DNA',
  'MMSplice',
  'MMSplice + SpliceMap',
  'MMSplice + SpliceMap + Ψ_ref',
  'SpliceAI',
  'SpliceAI + SpliceMap',
  'SpliceAI + SpliceMap + Ψ_ref'
)
color_list = absplice_model_colors[chosen_models]

# SI Fig. 6a
df <- fread('../data/als/pr_curve.csv')

g_als_pr_curve_faceted <- (
  ggplot(df, aes(x=recall, y=precision, color=model))
  + geom_step(direction = "vh")
  + theme_cowplot(font_size = fontsize)
  # + background_grid(major = "xy", minor = "xy")
  + scale_x_continuous(limits = c(-0.01, 0.5), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2))
  + scale_y_continuous(limits = c(-0.01, 0.3), breaks = seq(0, 1, 0.1), minor_breaks = seq(0, 1, 0.01))
  + scale_color_manual(values=color_list, guide=guide_legend(reverse = FALSE))
  + theme(
    legend.position = c(0.2, 0.7),
  )
  + labs(
    x='Recall',
    y='Precision',
    color=''
  )
  + facet_wrap('SpliceMap')
  + theme(
    strip.background = element_rect(colour="white", fill="white")
  )
  + scale_linetype_manual(guide = guide_legend(reverse = FALSE))
)
g_als_pr_curve_faceted
ggsave(plot=g_als_pr_curve_faceted, filename="out/Supplementary_Figures/Supplementary_Figure_6/Supplementary_Figure_6a.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)


# SI Fig. 6b
df_points <- fread('../data/als/tp_among_preds_points.csv')
df_line <- fread('../data/als/tp_among_preds_line.csv')


chosen_cutoffs <- c(
  "high",
  "medium",
  "low"
)
df_points <- df_points[model %in% chosen_models]
df_points$model <- factor(df_points$model, levels=chosen_models)
df_points$cutoff <- factor(df_points$cutoff, levels=chosen_cutoffs)

df_line <- df_line[model %in% chosen_models]
df_line$model <- factor(df_line$model, levels=chosen_models)

g_tp_total_als <- (
  ggplot()
  + geom_line(data=df_line, aes(x=rank_mean, y=tp_sum, color=model, alpha=0.5), size=1)
  + geom_point(data=df_points, aes(x=rank_mean, y=tp_sum, color=model, shape=cutoff), size=3)
  + scale_color_manual(values=color_list)
  + ylim(0,120)
  + xlim(0,6000)
  + theme_cowplot(font_size = 14)
  + labs(x='Total predictions',
         y='True positive \npredictions',
         color='',
         shape='cutoff'
  )
  + theme(
    legend.direction = "vertical",
    legend.position=c(0.65,0.3),
  )
  + guides(
    shape = guide_legend(
      nrow=3
    ),
    color = guide_legend(
      title = element_blank(),
      reverse = TRUE, 
      nrow=3
    )
  )
  + guides(
    color = guide_legend(order = 1, reverse = TRUE),
    shape = guide_legend(order = 2)
  )
  + guides(
    alpha = "none"
  )
  
)
g_tp_total_als
ggsave(plot=g_tp_total_als, filename="out/Supplementary_Figures/Supplementary_Figure_6/Supplementary_Figure_6b.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)


# SI Fig. 6c
df <- fread('../data/als/known_als_genes.csv')

chosen_models <- c(
  "MMSplice",                                     
  "MMSplice + SpliceMap",
  "MMSplice + SpliceMap + Ψ_ref",
  "SpliceAI",
  "SpliceAI + SpliceMap",
  "SpliceAI + SpliceMap + Ψ_ref",
  "AbSplice-DNA \n(SpliceMap \nGTEx brain)",
  "AbSplice-DNA \n(SpliceMap ALS \nmotor neuron)",
  "random"
)

color_list = absplice_model_colors[chosen_models]

df$model <- factor(df$model, levels=chosen_models)
df$model_name <- paste0(df$model, '\nn=', df$above_cutoff_all)

df_mapping <- df %>% select(model, model_name) %>% unique()
label_dict <- setNames(df_mapping$model_name, df_mapping$model)

# df$model_name <- factor(df$model_name, levels=df$model_name[order(chosen_models)], ordered=TRUE)
df$pval_rounded <- round(df$pval, 2)
df$pval_rounded_label <- ifelse(df$pval_rounded <= 0.05, paste0('P=', df$pval_rounded), 'N.S.')
# df$pval_rounded <- factor(df$pval_rounded)

g_als_enrichment <- (
  ggplot(df[df$model != 'random'])
  + aes_string(
    x='model_name', 
    y='enrichment', 
    fill='model', 
    # label="significance_label"
    label="pval_rounded_label"
    )
  + geom_text(
    nudge_y = +0.08,
    # nudge_x = -1,
    size=4, color="black"
  )
  + geom_bar(
    stat = "identity", 
    width=0.5,
    position = "dodge2"
  )
  + scale_fill_manual(values=color_list)
  + theme_cowplot(font_size = fontsize)
  + scale_x_discrete(labels = function(x) {
    x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
  })
  # + background_grid()
  + labs(
    x=element_blank(),
    y="Enrichment \nin ALS genes",
    title=element_blank()
  )
  + guides(
    fill = "none"
  )
  + geom_hline(yintercept=1, linetype="dashed")
  # + scale_x_discrete(labels = label_dict)
)
g_als_enrichment
ggsave(plot=g_als_enrichment, filename="out/Supplementary_Figures/Supplementary_Figure_6/Supplementary_Figure_6c.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)


# SI Fig. 6d
df <- fread('../data/als/proteomics.csv')

chosen_cutoffs <- c(
  "high"
)

df <- df[cutoff_category %in% chosen_cutoffs] 

library(binom)
binom.confint(x = c(18), n = 58, tol = 1e-8, method="prop.test")
df_test <- binom.confint(
  x = df$num_protein_low, 
  n = df$num, 
  tol = 1e-8, 
  # method="prop.test"
  method="exact"
)

df_plot <- cbind(df, df_test)
df_plot$model_name <- paste0(df_plot$model, '\nn=', df_plot$num)


g_als_proteomics_true_preds <- (
  ggplot(df_plot, aes(
    x=model_name,
    y=proportion_true
  ))
  + scale_fill_manual(values=color_list)
  + scale_color_manual(values=color_list)
  + geom_bar(aes(
    fill=model
  ),
  stat = "identity", 
  width=0.5,
  position = "dodge2",
  size=0.3
  )
  + geom_errorbar(
    aes(
      ymin=lower,
      ymax=upper
    ),
    position = position_dodge(width = 0.9),
    width = 0.2
  )
  + scale_x_discrete(labels = function(x) {
    x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
  })
  + theme_cowplot(font_size = fontsize)
  + background_grid()
  + labs(
    x=element_blank(),
    y="Proportion \nhigh impact predictions \nvalidated via proteomics",
    title=element_blank(),
    color='SpliceMap',
    fill=element_blank()
  )
  + guides(
    fill = "none"
  )
  # + coord_flip()
)
g_als_proteomics_true_preds
ggsave(plot=g_als_proteomics_true_preds, filename="out/Supplementary_Figures/Supplementary_Figure_6/Supplementary_Figure_6d.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)


# SI Fig. 6
SI_fig_6 <- ggarrange(
  nrow = 4, heights=c(1,1,1,1),
  ggarrange(labels = c('a'),
            g_als_pr_curve_faceted
  ),
  ggarrange(labels = c('b'),
            g_tp_total_als
  ),
  ggarrange(labels = c('c'),
            g_als_enrichment
  ),
  ggarrange(labels = c('d'),
            g_als_proteomics_true_preds
  )
)
SI_fig_6

ggsave(plot=SI_fig_6, filename="out/Supplementary_Figures/Supplementary_Figure_6/Supplementary_Figure_6.svg", height=unit(16, "cm"), width=unit(14, "cm"), dpi=450)


