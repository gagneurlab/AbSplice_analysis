source("config.R")

# Fig. 5a
sashimi_cat <- readPNG('../data/sashimi/clinically_accessible_tissue.png')
g_sashimi_cat <- as.raster(sashimi_cat)
g_sashimi_cat <- rasterGrob(g_sashimi_cat, interpolate = FALSE)
ggsave(plot=g_sashimi_cat, filename="out/main_figures/Figure_5/Figure_5a.pdf", height=unit(7, "cm"), width=unit(6, "cm"), dpi=450)

# Fig. 5b
df <- fread('../data/gtex_v8/rna/shared_sites_cat.csv')

df <- filter(df, tissue_target!='Brain') # show all brain
df <- filter(df, tissue_target!='allBrainTissues') # show all brain
# df <- df %>% filter(!grepl('Brain_', tissue_target)) # show only Brain union

df <- df %>% filter(!grepl('Whole_Blood', tissue_target))
df <- df %>% filter(!grepl('Cells_Cultured_fibroblasts', tissue_target))
df <- df %>% filter(!grepl('Cells_EBV_transformed_lymphocytes', tissue_target))


df$tissue_target <- revalue(df$tissue_target, tissue_mapping, warn_missing = FALSE)
df$tissue_cat <- revalue(df$tissue_cat, tissue_mapping, warn_missing = FALSE)

df$tissue_target <- factor(df$tissue_target, levels=str_sort(unique(df$tissue_target), decreasing = TRUE))

fontsize_5b = 14

g_shared_sites_with_cat <- (ggplot(df)
                            + aes_string(x = "tissue_cat", "tissue_target")
                            + geom_tile(aes(fill = percentage))
                            + scale_fill_distiller(
                              palette = 8, 
                              direction=1
                            )
                            + theme_cowplot(font_size=fontsize_5b)
                            + theme(axis.text.x = element_text(angle = 30, hjust = 1, size=fontsize_5b),
                                    axis.text.y = element_text(angle = 0, hjust = 1, size=fontsize_5b),
                                    legend.position = 'right', 
                                    legend.title=element_text("Proportion"),
                                    legend.text=element_text(size=fontsize_5b)
                            )
                            + labs(
                              x=element_blank(),
                              y='Clinically non-accessible tissues'
                            )
                            + labs(fill = "Proportion of \nsplice sites \nin CAT")
)
g_shared_sites_with_cat
ggsave(plot=g_shared_sites_with_cat, filename="out/main_figures/Figure_5/Figure_5b.pdf", height=unit(12, "cm"), width=unit(9, "cm"), dpi=450)


# Fig. 5c
df <- fread('../data/gtex_v8/rna/pr_curve_all_tissues.csv')

chosen_models <- c(
  'MMSplice', 
  'SpliceAI',
  'AbSplice-DNA', 
  'CAT p-value (fibroblasts)',
  'AbSplice-RNA (fibroblasts)' 
)
color_list = absplice_model_colors[chosen_models]
color_list

font_size_cat <- 15

df <- process_models_pr_curve(df, chosen_models)

g_gtex_v8_RNA_pr_curve <- plot_pr_curve(
  df,
  ylim=1,
  color_list,
  breaks_y=0.2, minor_breaks_y=0.1,
)
g_gtex_v8_RNA_pr_curve
ggsave(plot=g_gtex_v8_RNA_pr_curve, filename="out/main_figures/Figure_5/Figure_5c.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 5d
df <- fread('../data/gtex_v8/rna/boxplot_across_tissues.csv')
# df <- df %>% filter(grepl('Brain_', tissue)) # only show Brain tissues

chosen_models <- c(
  'MMSplice',
  'SpliceAI',
  'AbSplice-DNA',
  # 'CAT p-value (blood)',
  'CAT p-value (fibroblasts)',
  # 'AbSplice-RNA (blood)',
  'AbSplice-RNA (fibroblasts)'
  
  # 'AbSplice-RNA (all CATs)'
)


color_list = absplice_model_colors[chosen_models]
df <- process_models_box_plot(df, chosen_models)

g_gtex_v8_RNA_box_plot <- plot_box_plot(
  df,
  color_list, 
  comparisons=list(
    c("AbSplice-RNA (fibroblasts)", "AbSplice-DNA"),
    c("CAT p-value (fibroblasts)", "AbSplice-DNA")
    
    # c("AbSplice-RNA (fibroblasts)", "AbSplice-RNA (blood)")
  ),
  x_lim=0.5,
  title="All tissues",
  coord_flip=FALSE
  # jitter=TRUE
)
g_gtex_v8_RNA_box_plot
ggsave(plot=g_gtex_v8_RNA_box_plot, filename="out/main_figures/Figure_5/Figure_5d.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Figure 5
fig_5 <- ggarrange(
  ncol = 1, nrow = 2, heights = c(5,3),
  ggarrange(ncol = 2, labels = c('a', 'b'), widths = c(1, 1),
            g_sashimi_cat,
            g_shared_sites_with_cat
  ),
  ggarrange(ncol = 2, labels = c('c', 'd'), widths = c(1,1),
            g_gtex_v8_RNA_pr_curve,
            g_gtex_v8_RNA_box_plot
  )
)
fig_5

ggsave(plot=fig_5, filename="out/main_figures/Figure_5/Figure_5.pdf", height=unit(14, "cm"), width=unit(13, "cm"), dpi=dpi)
ggsave(plot=fig_5, filename="out/main_figures/Figure_5/Figure_5.svg", height=unit(14, "cm"), width=unit(13, "cm"), dpi=dpi)
ggsave(plot=fig_5, filename="out/main_figures/Figure_5/Figure_5.png", height=unit(14, "cm"), width=unit(13, "cm"), dpi=dpi)
