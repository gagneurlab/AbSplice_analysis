source("config.R")

chosen_models <- c(
  "SpliceAI",
  "SpliceAI + SpliceMap",
  "SpliceAI + SpliceMap + Ψ_ref",
  "MMSplice",
  "MMSplice + SpliceMap",
  "MMSplice + SpliceMap + Ψ_ref",
  "AbSplice-DNA"
)
color_list = absplice_model_colors[chosen_models]

# Fig. 3a
sashimi_ref_psi <- readPNG('../data/sashimi/ref_psi.png')
g_sashimi_ref_psi <- as.raster(sashimi_ref_psi)
g_sashimi_ref_psi <- rasterGrob(g_sashimi_ref_psi, interpolate = FALSE)

ggsave(plot=g_sashimi_ref_psi, filename="out/main_figures/Figure_3/Figure_3a.pdf", height=unit(7, "cm"), width=unit(6, "cm"), dpi=450)

# Fig. 3b
df <- fread(gtex_v8_pr_curve)
df <- process_models_pr_curve(df, chosen_models)

df_thresholds <- fread('../data/gtex_v8/pr_curve_thresholds.csv')

chosen_models_points <- c(
  "SpliceAI",
  "MMSplice",
  "AbSplice-DNA"
)

chosen_cutoffs <- c(
  "high",
  "medium",
  "low"
)

df_points <- df_thresholds[model %in% chosen_models_points]
df_points$model <- factor(df_points$model, levels=chosen_models_points)
df_points$cutoff <- factor(df_points$cutoff, levels=chosen_cutoffs)

plot_pr_curve_thresholds <- function(df, df_points, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh")
    + geom_point(data=df_points, aes(x=recall, y=precision, color=model, shape=cutoff), size=3)
    + theme_cowplot(font_size = fontsize)
    # + background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + guides(
      color = guide_legend(order = 1, reverse = TRUE),
      shape = guide_legend(order = 2)
    )
    + scale_color_manual(values=color_list)
    + theme(
      legend.position = c(0.95, 1.1),
      legend.justification = c("right", "top")
    )
    + labs(
      x='Recall', 
      y='Precision', 
      color='', 
      shape='cutoff',
      title=title
    )
  )
  return(g)
}

g_gtex_v8_DNA_pr_curve_fig3 <- plot_pr_curve_thresholds(
  df,
  df_points,
  ylim=0.4,
  color_list
)
ggsave(plot=g_gtex_v8_DNA_pr_curve_fig3, filename="out/main_figures/Figure_3/Figure_3b.svg", height=unit(4, "cm"), width=unit(9, "cm"), dpi=450)

# Fig. 3c
df <- fread(gtex_v8_boxplot)
df <- process_models_box_plot(df, chosen_models)

g_gtex_v8_DNA_box_plot_fig3 <- plot_box_plot(
  df,
  color_list, 
  comparisons=list(
    c("AbSplice-DNA", "MMSplice + SpliceMap + Ψ_ref"),
    c("MMSplice + SpliceMap + Ψ_ref", "MMSplice + SpliceMap"),
    c("SpliceAI + SpliceMap + Ψ_ref", "SpliceAI + SpliceMap")
  ),
  coord_flip=TRUE
)
g_gtex_v8_DNA_box_plot_fig3
ggsave(plot=g_gtex_v8_DNA_box_plot_fig3, filename="out/main_figures/Figure_3/Figure_3c.pdf", height=unit(4, "cm"), width=unit(9, "cm"), dpi=450)


# Fig. 3d
chosen_models <- c(
  "SQUIRLS",
  "CADD-Splice",
  "MMSplice",
  "SpliceAI",
  "AbSplice-DNA"
)
color_list = absplice_model_colors[chosen_models]

var_category_order <- c(
  "Splice acceptor",
  "Splice donor",
  "Splice region",
  "Exon",
  "Intron",
  "All"
)

plot_pr_curve_variant <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh")
    + scale_color_manual(
      values=color_list, 
      guide=guide_legend(
        reverse = TRUE
        )
      )
    + theme_cowplot(font_size = fontsize)
    # + background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(x='Recall', y='Precision', color='', title=title)
    + facet_wrap('var_category', nrow=1, scales='free_y', labeller=labeller(var_category=mapping_list))
    # + theme(legend.text = element_text(
    #   margin = margin(r = 30, unit = "pt")))
    # + panel_border()
    + theme(
      legend.key = element_blank(),
      legend.position = "none",
      strip.background = element_rect(colour="white", fill="white")
    )
  )
  return(g)
}

df <- fread('../data/gtex_v8/var_categories.csv')
df_stats <- fread('../data/gtex_v8/var_categories_stats.csv')
df <- merge(df, df_stats[, c('var_category', 'num_outliers')], by='var_category', all.x =TRUE)

df$var_category_name <- paste0(df$var_category, '\nn=', df$num_outliers)

df <- process_models_pr_curve(df, chosen_models)

# define the order of the facets
df$var_category <- factor(df$var_category, levels=var_category_order)
df_mapping <- df %>% select(var_category, var_category_name) %>% unique()
mapping_list <- setNames(df_mapping$var_category_name, df_mapping$var_category)

g_gtex_v8_DNA_pr_curve_var_categories_all <- plot_pr_curve_variant(
  df,
  ylim=0.5,
  color_list, 
  title=''
)
ggsave(plot=g_gtex_v8_DNA_pr_curve_var_categories_all, filename="out/main_figures/Figure_3/Figure_3d.pdf", height=unit(4, "cm"), width=unit(15, "cm"), dpi=450)


# Fig. 3e
plot_pr_curve_outlier <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='All tissues') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh")
    + scale_color_manual(
      values=color_list
      )
    + theme_cowplot(font_size = fontsize)
    # + background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(
      x='Recall', 
      y='Precision', 
      title=title
      )
    + facet_wrap('outlier_category', nrow=1, scales='free_y', labeller=labeller(outlier_category=mapping_list))
    + theme(
      legend.key = element_blank(),
      strip.background = element_rect(colour="white", fill="white")
    )
    + theme(legend.position = "none")
  )
  return(g)
}

outlier_category_order <- c(
  'Exon elongation',
  'Exon truncation',
  'Exon skipping',
  'Any splicing efficiency outlier',
  'Any alternative donor or acceptor choice',
  'All'
)
outlier_category_rename <- c(
  "Exon elongation" = "Exon elongation",
  "Exon truncation" = "Exon truncation",
  "Exon skipping" = "Exon skipping",  
  "Any splicing efficiency outlier" = "Any splicing \nefficiency outlier",
  "Any alternative donor or acceptor choice" = "Any alternative donor \nor acceptor choice",
  "All" = "All"
)

df <- fread('../data/gtex_v8/outlier_categories.csv')
df_stats <- fread('../data/gtex_v8/outlier_categories_stats.csv')
df <- merge(df, df_stats[, c('outlier_category', 'num_outliers')], by='outlier_category', all.x =TRUE)

df$outlier_category_rename <- revalue(df$outlier_category, outlier_category_rename, warn_missing = FALSE)
df$outlier_category_name <- paste0(df$outlier_category_rename, '\nn=', df$num_outliers)

df <- process_models_pr_curve(df, chosen_models)

# define the order of the facets
df$outlier_category <- factor(df$outlier_category, levels=outlier_category_order)
df_mapping <- df %>% select(outlier_category, outlier_category_name) %>% unique()
mapping_list <- setNames(df_mapping$outlier_category_name, df_mapping$outlier_category)

g_gtex_v8_DNA_pr_curve_outlier_categories_all <- plot_pr_curve_outlier(
  df,
  ylim=0.5,
  color_list,
  title=''
)
ggsave(plot=g_gtex_v8_DNA_pr_curve_outlier_categories_all, filename="out/main_figures/Figure_3/Figure_3e.pdf", height=unit(4, "cm"), width=unit(15, "cm"), dpi=450)

# Fig 3
fig_3 <- ggarrange(
    nrow = 2, heights = c(1,1),
    ggarrange(ncol = 2, widths = c(0.8, 1),
              ggarrange(g_sashimi_ref_psi, labels = c('a')),
              ggarrange(nrow = 2, labels = c('b', 'c'), heights = c(1,1),
                         g_gtex_v8_DNA_pr_curve_fig3,
                         g_gtex_v8_DNA_box_plot_fig3
                       )
             ),
    ggarrange(
      ncol=1, nrow=2, labels=c('d', 'e'), common.legend = TRUE, legend="bottom",
      g_gtex_v8_DNA_pr_curve_var_categories_all, 
      g_gtex_v8_DNA_pr_curve_outlier_categories_all
    )
)

ggsave(plot=fig_3, filename="out/main_figures/Figure_3/Figure_3.png", height=unit(16, "cm"), width=unit(15, "cm"), dpi=dpi)
ggsave(plot=fig_3, filename="out/main_figures/Figure_3/Figure_3.svg", height=unit(16, "cm"), width=unit(15, "cm"), dpi=dpi)



