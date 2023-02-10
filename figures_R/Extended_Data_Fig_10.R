source("config.R")

# Extended Data Figure 10a
df <- fread('../data/gtex_v8/rna/all_cat_compared/cats_single.csv')

chosen_models <- c(
  'AbSplice-DNA',
  'CAT p-value (blood)',
  'CAT p-value (fibroblasts)',
  'CAT p-value (lymphocytes)',
  'AbSplice-RNA (blood)',
  'AbSplice-RNA (fibroblasts)',
  'AbSplice-RNA (lymphocytes)'
)

absplice_model_names_cat <- c(
  'AbSplice-DNA' = 'AbSplice-DNA',
  'CAT p-value (blood)' = 'CAT p-value',
  'CAT p-value (fibroblasts)' = 'CAT p-value',
  'CAT p-value (lymphocytes)' = 'CAT p-value',
  'AbSplice-RNA (blood)' = 'AbSplice-RNA',
  'AbSplice-RNA (fibroblasts)' = 'AbSplice-RNA',
  'AbSplice-RNA (lymphocytes)' = 'AbSplice-RNA'
)

cat_map <- c(
  "Blood" = "Blood (N=12975)",
  "Fibroblasts" = "Fibroblasts (N=8688)",
  "Lymphocytes" = "Lymphocytes (N=3268)"
)

df$model <- revalue(df$model, absplice_model_names_cat, warn_missing = FALSE)
df$cat <- revalue(df$cat, cat_map, warn_missing = FALSE)

df <- process_models_pr_curve(df, c("AbSplice-DNA", "CAT p-value",  "AbSplice-RNA"))

color_list = absplice_model_colors[absplice_model_names_cat[chosen_models]]
color_list

ylim=1
breaks_x=0.2
minor_breaks_x=0.2
breaks_y=0.2
minor_breaks_y=0.1
title=''

g_compare_cats <- (
  ggplot(df, aes(x=recall, y=precision, color=model))
  + geom_step(direction = "vh")
  + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE))
  + theme_cowplot(font_size = fontsize)
  + facet_wrap('cat')
  # + background_grid(major = "xy", minor = "xy")
  + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
  + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
  + theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )
  + labs(x='Recall', y='Precision', color='', title=title)
  + theme(
    strip.background = element_rect(colour="white", fill="white")
  )
)
# g_compare_cats
ggsave(plot=g_compare_cats, filename="out/Extended_Data_Figures/Extended_Data_Figure_10/Extended_Data_Figure_10a.pdf", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)

# Extended Data Figure 10b
df <- fread('../data/gtex_v8/rna/all_cat_compared/blood_fibro.csv')

chosen_models <- c(
  # 'AbSplice-DNA',
  'CAT p-value (blood)',
  'CAT p-value (fibroblasts)',
  'AbSplice-RNA (blood)',
  'AbSplice-RNA (fibroblasts)',
  'AbSplice-RNA (all CATs)'
)

absplice_model_colors_d = c(
  # 'AbSplice-DNA' = "#d62728",
  'CAT p-value (blood)' ='lightblue',
  'CAT p-value (fibroblasts)' = "#ffd33f",
  'AbSplice-RNA (blood)' = 'royalblue',
  'AbSplice-RNA (fibroblasts)' = "#b28c4d",
  'AbSplice-RNA (all CATs)' = 'brown'
)

color_list = absplice_model_colors_d[chosen_models]
color_list

df <- process_models_pr_curve(df, chosen_models)

g_blood_fibro <- plot_pr_curve(
  df,
  ylim=1,
  color_list,
  breaks_y=0.2, minor_breaks_y=0.1,
  title='Blood & Fibroblasts (N=7490)'
)
g_blood_fibro

g_blood_fibro <- g_blood_fibro + theme(plot.title = element_text(face = "plain"))
ggsave(plot=g_blood_fibro, filename="out/Extended_Data_Figures/Extended_Data_Figure_10/Extended_Data_Figure_10b.pdf", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)


# Extended Data Figure 10c
df <-fread('../data/gtex_v8/rna/cat_expression_levels.csv')

chosen_models <- c(
  'MMSplice',
  'SpliceAI',
  'AbSplice-DNA',
  'CAT p-value (fibroblasts)',
  'AbSplice-RNA (fibroblasts)'
)

color_list = absplice_model_colors[chosen_models]
df <- process_models_pr_curve(df, chosen_models)

plot_pr_curve <- function(df, color_list, ylim=0.2, breaks_x=0.2, minor_breaks_x=0.2, breaks_y=0.1, minor_breaks_y=0.02, title='') {
  g <- (
    ggplot(df, aes(x=recall, y=precision, color=model))
    + geom_step(direction = "vh")
    + scale_color_manual(values=color_list, guide=guide_legend(reverse = TRUE))
    + theme_cowplot(font_size = fontsize)
    # + background_grid(major = "xy", minor = "xy")
    + scale_x_continuous(limits = c(-0.01, 1.05), breaks = seq(0, 1, breaks_x), minor_breaks = seq(0, 1, minor_breaks_x))
    + scale_y_continuous(limits = c(-0.01, ylim), breaks = seq(0, 1, breaks_y), minor_breaks = seq(0, 1, minor_breaks_y))
    + labs(
      x='Recall', 
      y='Precision', 
      color='',
      title=title
    )
    + facet_wrap('above_cutoff')
    + theme(
      legend.key = element_blank(),
      # legend.position = "none",
      strip.background = element_rect(colour="white", fill="white")
    )
  )
  return(g)
}

Ext_Data_Fig_10c <- plot_pr_curve(
  df,
  ylim=1,
  color_list,
  breaks_y=0.2, minor_breaks_y=0.1,
)
ggsave(plot=Ext_Data_Fig_10c, filename="out/Extended_Data_Figures/Extended_Data_Figure_10/Extended_Data_Figure_10c.pdf", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)

# Extended Data Figure 10
Ext_Data_Fig_10 <-ggarrange(
   ncol = 1, nrow = 3, labels = c('a', 'b', 'c'),
   g_compare_cats,
   g_blood_fibro,
   Ext_Data_Fig_10c
)
ggsave(plot=Ext_Data_Fig_10, filename="out/Extended_Data_Figures/Extended_Data_Figure_10/Extended_Data_Figure_10.pdf", height=unit(15, "cm"), width=unit(12, "cm"), dpi=450)
ggsave(plot=Ext_Data_Fig_10, filename="out/Extended_Data_Figures/Extended_Data_Figure_10/Extended_Data_Figure_10.png", height=unit(15, "cm"), width=unit(12, "cm"), dpi=450)
