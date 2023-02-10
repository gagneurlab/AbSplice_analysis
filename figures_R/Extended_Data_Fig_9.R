source("config.R")

# Extended Data Figure 9a
chosen_models <- c(
  "SQUIRLS",
  "CADD-Splice",
  'MTSplice',
  'MMSplice',
  'SpliceAI'
)
color_list = absplice_model_colors[chosen_models]
color_list

df <- fread('../data/gtex_v8/pr_curve_all_tissues.csv')
df <- process_models_pr_curve(df, chosen_models)

# df <- sample_frac(df, 0.01)

Ext_Data_Fig_9a <- plot_pr_curve(
  df,
  ylim=0.2,
  color_list
)
Ext_Data_Fig_9a
ggsave(plot=Ext_Data_Fig_9a, filename="out/Extended_Data_Figures/Extended_Data_Figure_9/Extended_Data_Figure_9a.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Extended Data Figure 9b
df <- fread('../data/gtex_v8/absplice_different_features.csv')

color_list <- c(
  'AbSplice-DNA \n(+ CADD-Splice)' = 'grey',
  'AbSplice-DNA' = 'red',
  'AbSplice-DNA \n(+ SQUIRLS)' = 'grey',
  'AbSplice-DNA \n(+ MTSplice)' = 'grey'
)

plot_box_plot <- function(
  df, color_list, title="Across tissues",
  comparisons = list(
    c('AbSplice-DNA \n(+ CADD-Splice)',
      'AbSplice-DNA')
  ),
  coord_flip=FALSE,
  x_lim=max(df$`Average Precision Score`) * 1.3
) {
  g <- (
    ggplot(df, aes(x=`Average Precision Score`, y=reorder(model, mean), fill=model))
    + scale_fill_manual(values=color_list)
    + geom_boxplot()
    + xlim(0, x_lim)
    + theme_cowplot(font_size = fontsize)
    # + background_grid()
    + theme(legend.title = element_blank())
    + theme(legend.position = "none")
    + labs(
      x="auPRC",
      y=element_blank(),
      title=title
    )
    + {if (coord_flip == TRUE)coord_flip()}
    + stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      paired = TRUE,
      method.args = list(alternative = "two.sided")
    )
  )
  return(g)
}

Ext_Data_Fig_9b <- plot_box_plot(
  df,
  color_list,
  comparisons=list(
    c(
      'AbSplice-DNA \n(+ CADD-Splice)',
      'AbSplice-DNA'
    ),
    c(
      'AbSplice-DNA \n(+ MTSplice)',
      'AbSplice-DNA'
    ),
    c(
      'AbSplice-DNA \n(+ SQUIRLS)',
      'AbSplice-DNA'
    )
  ),
  coord_flip=FALSE,
  x_lim = 0.4
)
Ext_Data_Fig_9b
ggsave(plot=Ext_Data_Fig_9b, filename="out/Extended_Data_Figures/Extended_Data_Figure_9/Extended_Data_Figure_9b.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Extended Data Figure 9c
chosen_models <- c(
  "AbSplice-DNA (logisitic regression)",
  "AbSplice-DNA (random forest)",
  "AbSplice-DNA (GAM)"
)
color_list = absplice_model_colors[chosen_models]
color_list

df <- fread('../data/gtex_v8/absplice_different_classifiers/pr_curve.csv')
df <- process_models_pr_curve(df, chosen_models)
# df <- sample_frac(df, 0.01)

Ext_Data_Fig_9c <- plot_pr_curve(
  df,
  ylim=0.5,
  color_list
)
Ext_Data_Fig_9c
ggsave(plot=Ext_Data_Fig_9c, filename="out/Extended_Data_Figures/Extended_Data_Figure_9/Extended_Data_Figure_9c.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)

# Extended Data Figure 9d
df <- fread('../data/gtex_v8/absplice_different_classifiers/boxplot.csv')

plot_box_plot <- function(
  df, color_list, title="Across tissues",
  comparisons = list(
    c('AbSplice-DNA (GAM)', 'AbSplice-DNA (random forest)')
  ),
  coord_flip=FALSE,
  x_lim=max(df$`Average Precision Score`) * 1.3
) {
  g <- (
    ggplot(df, aes(x=`Average Precision Score`, y=reorder(model, mean), fill=model))
    + geom_boxplot()
    + xlim(0, x_lim)
    + scale_y_discrete(labels = function(x) {
      x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
    })
    + theme_cowplot(font_size = fontsize)
    + background_grid()
    + theme(legend.title = element_blank())
    + theme(legend.position = "none")
    + labs(
      x="auPRC",
      y=element_blank(),
      title=title
    )
    + {if (coord_flip == TRUE)coord_flip()}
    + stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      paired = TRUE,
      method.args = list(alternative = "greater")
    )
  )
  return(g)
}

Ext_Data_Fig_9d <- plot_box_plot(
  df,
  comparisons=list(
    c('AbSplice-DNA (GAM)'),
    c('AbSplice-DNA (random forest)')
  )
)
Ext_Data_Fig_9d
ggsave(plot=Ext_Data_Fig_9d, filename="out/Extended_Data_Figures/Extended_Data_Figure_9/Extended_Data_Figure_9d.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)

# Extended Data Figure 9
Ext_Data_Fig_9 <- ggarrange(
  nrow = 2, heights=c(1,1),
  ggarrange(
    ncol = 2, labels = c('a', 'b'),
    Ext_Data_Fig_9a,
    Ext_Data_Fig_9b
  ),
  ggarrange(ncol=2, widths = c(1, 1), labels = c('c', 'd'),
            Ext_Data_Fig_9c,
            Ext_Data_Fig_9d
  )
)
ggsave(plot=Ext_Data_Fig_9, filename="out/Extended_Data_Figures/Extended_Data_Figure_9/Extended_Data_Figure_9.png", height=unit(10, "cm"), width=unit(12, "cm"), dpi=450)
ggsave(plot=Ext_Data_Fig_9, filename="out/Extended_Data_Figures/Extended_Data_Figure_9/Extended_Data_Figure_9.pdf", height=unit(10, "cm"), width=unit(12, "cm"), dpi=450)

