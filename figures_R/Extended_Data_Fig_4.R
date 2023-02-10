source("config.R")

chosen_models <- c(
  'MMSplice',
  'MMSplice + SpliceMap',
  'MMSplice + SpliceMap + Ψ_ref',
  'SpliceAI',
  'SpliceAI + SpliceMap',
  'SpliceAI + SpliceMap + Ψ_ref',
  "AbSplice-DNA"
)

color_list = absplice_model_colors[chosen_models]
color_list

df <- fread('../data/gtex_v8/outlier_filters.csv')
df <- process_models_pr_curve(df, chosen_models)

# Extended Data Fig. 4a
df_plot <- df[df$outlier_filter %in% c('default cutoffs')]

Ext_Data_Fig_4a <- plot_pr_curve(
  df_plot,
  ylim=0.4,
  color_list,
  title='Filter 1: \ndefault cutoffs'
)
ggsave(plot=Ext_Data_Fig_4a, filename="out/Extended_Data_Figures/Extended_Data_Figure_4/Extended_Data_Figure_4a.svg", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)
Ext_Data_Fig_4a <- Ext_Data_Fig_4a + theme(legend.position = "none")

# Extended Data Fig. 4b
df_plot <- df[df$outlier_filter %in% c('default cutoffs \n& multiple tissues')]

Ext_Data_Fig_4b <- plot_pr_curve(
  df_plot,
  ylim=0.4,
  color_list,
  title='Filter 2: \ndefault cutoffs \n& multiple tissues'
)
ggsave(plot=Ext_Data_Fig_4b, filename="out/Extended_Data_Figures/Extended_Data_Figure_4/Extended_Data_Figure_4b.svg", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)
Ext_Data_Fig_4b <- Ext_Data_Fig_4b + theme(legend.position = "none")

# Extended Data Fig. 4c
df_plot <- df[df$outlier_filter %in% c('default cutoffs \n& rare variants')]

Ext_Data_Fig_4c <- plot_pr_curve(
  df_plot,
  ylim=0.4,
  color_list,
  title='Filter 3: \ndefault cutoffs \n& rare variants'
)
ggsave(plot=Ext_Data_Fig_4c, filename="out/Extended_Data_Figures/Extended_Data_Figure_4/Extended_Data_Figure_4c.svg", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)
Ext_Data_Fig_4c <- Ext_Data_Fig_4c + theme(legend.position = "none")

# Extended Data Fig. 4d
df <- fread('../data/gtex_v8/replication.csv')
df <- process_models_pr_curve(df, chosen_models)

plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh")
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE))
    + theme_cowplot(font_size = fontsize)
    # + background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + theme(
      legend.position = c(0.15, 0.6),
    )
    + labs(x='Recall', y='Precision', color='', title=title)
    + facet_wrap('replication')
    + theme(
      legend.key = element_blank(),
      strip.background = element_rect(colour="white", fill="white")
    )
  )
  return(g)
}

Ext_Data_Fig_4d <- plot_pr_curve(
  df,
  ylim=0.4,
  color_list,
  title=''
)

ggsave(plot=Ext_Data_Fig_4d, filename="out/Extended_Data_Figures/Extended_Data_Figure_4/Extended_Data_Figure_4d.svg", height=unit(4, "cm"), width=unit(12, "cm"), dpi=450)


# Extended Data Fig. 4
Ext_Data_Fig_4 <- ggarrange(nrow = 2,
                      ggarrange(
                        ncol = 3,
                        ggarrange(labels = c('a'),
                                  Ext_Data_Fig_4a
                        ),
                        ggarrange(labels = c('b'),
                                  Ext_Data_Fig_4b
                        ),
                        ggarrange(labels = c('c'),
                                  Ext_Data_Fig_4c
                        )
                      ),
                      ggarrange(labels = c('d'),
                                Ext_Data_Fig_4d
                      )
)

ggsave(plot=Ext_Data_Fig_4, filename="out/Extended_Data_Figures/Extended_Data_Figure_4/Extended_Data_Figure_4.svg", height=unit(8, "cm"), width=unit(12, "cm"), dpi=450)
