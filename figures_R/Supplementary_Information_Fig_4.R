source("config.R")

chosen_models <- c(
  'MMSplice',
  'MMSplice + SpliceMap',
  'MMSplice + SpliceMap + Ψ_ref',
  'SpliceAI',
  'SpliceAI + SpliceMap',
  'SpliceAI + SpliceMap + Ψ_ref',
  'AbSplice-DNA'
)
color_list = absplice_model_colors[chosen_models]
color_list

# SI Fig. 4a
df <- fread('../data/mito/pr_jackknife.csv')
df$model <- factor(df$model, levels=chosen_models)

g_mito_jackknife_revision <- (
  ggplot(df, aes(x=model, y=mean))
  + scale_fill_manual(values=color_list, guide=guide_legend(reverse = TRUE))
  + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE))
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
      ymin=mean - sem,
      ymax=mean + sem
    ),
    position = position_dodge(width = 0.9),
    width = 0.2
  )
  + scale_x_discrete(labels = function(x) {
    x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
  })
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    x=element_blank(),
    y="auPRC",
    title=element_blank(),
    color='SpliceMap',
    fill=element_blank()
  )
  + guides(
    fill = "none"
  )
)
g_mito_jackknife_revision
ggsave(plot=g_mito_jackknife_revision, filename="out/Supplementary_Figures/Supplementary_Figure_4/Supplementary_Figure_4a.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)


# SI Fig. 4b
df_points <- fread('../data/mito/tp_among_preds_points.csv')
df_line <- fread('../data/mito/tp_among_preds_line.csv')

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

g_tp_total_mito <- (
  ggplot()
  + geom_line(data=df_line, aes(x=rank_mean, y=tp_sum, color=model, alpha=0.5), size=1)
 + geom_point(data=df_points, aes(x=rank_mean, y=tp_sum, color=model, shape=cutoff), size=3)
 + scale_color_manual(values=color_list)
 + ylim(0,25)
 + xlim(0,800)
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
g_tp_total_mito
ggsave(plot=g_tp_total_mito, filename="out/Supplementary_Figures/Supplementary_Figure_4/Supplementary_Figure_4b.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)


# SI Fig. 4
g_mito_all_models <- ggarrange(
  ggarrange(nrow=2, widths = c(1, 1), labels = c('a', 'b'),
            g_mito_jackknife_revision,
            g_tp_total_mito
  )
)
g_mito_all_models

ggsave(plot=g_mito_all_models, filename="out/Supplementary_Figures/Supplementary_Figure_4/Supplementary_Figure_4.svg", height=unit(12, "cm"), width=unit(12, "cm"), dpi=450)

