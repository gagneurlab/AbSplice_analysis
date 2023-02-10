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

# Extended Data Figure 1a
plot_box_plot <- function(
  df, color_list, title="Across tissues",
  comparisons = list(
    c("AbSplice-DNA", "MMSplice + SpliceMap + Ψ_ref")
  ),
  y_lim=max(df$`Average Precision Score`) * 1.3
) {
  g <- (
    ggplot(df, aes(y=`Average Precision Score`, x=model, fill=model))
    + scale_fill_manual(values=color_list)
    + geom_boxplot()
    + facet_wrap('outlier_origin', nrow=3)
    + ylim(0, y_lim)
    + scale_x_discrete(labels = function(x) {
      x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
    })
    
    + theme_cowplot(font_size = fontsize)
    # + background_grid()
    + theme(
      legend.position = "none",
      legend.title = element_blank()
    )
    + labs(
      y="auPRC",
      x=element_blank(),
      title=title
    )
    + stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      paired = TRUE,
      method.args = list(alternative = "greater")
    )
    + coord_flip()
    + theme(
      legend.key = element_blank(),
      legend.position = "none",
      strip.background = element_rect(colour="white", fill="white")
    )
  )
  return(g)
}

chosen_outliers <- c(
  "LeafcutterMD",
  "SPOT",
  "FRASER"
)

df <- fread('../data/gtex_v8/leafcutter_spot_fraser.csv')
df <- df[df$model %in% chosen_models]
df <- df[df$outlier_origin %in% chosen_outliers]

df <- process_models_box_plot(df, chosen_models)

Ext_Data_Fig_1a <- plot_box_plot(
  df,
  color_list, 
  comparisons=list(
    c("AbSplice-DNA", "SpliceAI + SpliceMap + Ψ_ref")
  ),
  y_lim = 0.07,
  title=''
)
Ext_Data_Fig_1a
ggsave(plot=Ext_Data_Fig_1a, filename="out/Extended_Data_Figures/Extended_Data_Figure_1/Extended_Data_Figure_1a.pdf", height=unit(15, "cm"), width=unit(12, "cm"), dpi=450)

# Extended Data Figure 1b
df <- fread('../data/gtex_v8/delta_psi_cutoff.csv')
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
    + labs(x='Recall', y='Precision', color='', title=title)
    + facet_wrap('delta_psi_cutoff_name')
    # + facet_wrap('var_category', nrow=1, scales='free_y', labeller=labeller(var_category=mapping_list))
    + theme(
      legend.key = element_blank(),
      legend.position = "right",
      strip.background = element_rect(colour="white", fill="white")
    )
  )
  return(g)
}

# df <- sample_frac(df, 0.01) 
df$delta_psi_cutoff_name <- paste0('ΔΨ = ', df$delta_psi_cutoff)

Ext_Data_Fig_1b <- plot_pr_curve(
  df,
  ylim=0.5,
  color_list,
  title=''
)
Ext_Data_Fig_1b
ggsave(plot=Ext_Data_Fig_1b, filename="out/Extended_Data_Figures/Extended_Data_Figure_1/Extended_Data_Figure_1b.svg", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)

# Extended Data Figure 1
Ext_Data_Fig_1 <- ggarrange(
  nrow = 2, labels = c('a', 'b'), heights = c(3,1),
  Ext_Data_Fig_1a,
  Ext_Data_Fig_1b
)
Ext_Data_Fig_1


ggsave(plot=Ext_Data_Fig_1, filename="out/Extended_Data_Figures/Extended_Data_Figure_1/Extended_Data_Figure_1.svg", height=unit(15, "cm"), width=unit(12, "cm"), dpi=450)
