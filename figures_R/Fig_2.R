source("config.R")

chosen_models <- c(
  "MMSplice",
  "MMSplice + SpliceMap",
  "SpliceAI",
  "SpliceAI + SpliceMap"
)
color_list = absplice_model_colors[chosen_models]

# Fig. 2a-c
sashimi_gene_unexpressed <- readPNG('../data/sashimi/unexpressed_gene.png')
sashimi_exon_unused <- readPNG('../data/sashimi/unexpressed_exon.png')
sashimi_weak_site <- readPNG('../data/sashimi/weak_site.png')

fig_2a_sashimi_gene_unexpressed <- as.raster(sashimi_gene_unexpressed)
fig_2a_sashimi_gene_unexpressed <- rasterGrob(fig_2a_sashimi_gene_unexpressed, interpolate = FALSE)

fig_2b_sashimi_exon_unused <- as.raster(sashimi_exon_unused)
fig_2b_sashimi_exon_unused <- rasterGrob(fig_2b_sashimi_exon_unused, interpolate = FALSE)

fig_2c_sashimi_weak_site <- as.raster(sashimi_weak_site)
fig_2c_sashimi_weak_site <- rasterGrob(fig_2c_sashimi_weak_site, interpolate = FALSE)

ggsave(plot=fig_2a_sashimi_gene_unexpressed, filename="out/main_figures/Figure_2/Figure_2a.pdf", height=unit(7, "cm"), width=unit(6, "cm"), dpi=450)
ggsave(plot=fig_2b_sashimi_exon_unused, filename="out/main_figures/Figure_2/Figure_2b.pdf", height=unit(7, "cm"), width=unit(6, "cm"), dpi=450)
ggsave(plot=fig_2c_sashimi_weak_site, filename="out/main_figures/Figure_2/Figure_2c.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 2d
df <- fread('../data/gtex_v8/venn_gencode.csv')
splice_site_gencode <- df[df$annotation_source == 'gencode', splice_site]
splice_site_splicemap <- df[df$annotation_source == 'splicemap', splice_site]
venn_diag <- list(
  'GENCODE' = splice_site_gencode,
  'SpliceMap\n(all GTEx tissues)' = splice_site_splicemap
)

fig_2d_gtex_gencode_venn = as.ggplot(plot(
  euler(
    venn_diag
  ), 
  key = TRUE, 
  counts = TRUE, 
  quantities = list(
    type = c("counts", "percent"), 
    font=fontsize,
    round=2
  ), 
  fills =list(
    fill=c("#e7b8b9", "#bbdcc2", "#e7ddcd")
  ), 
  edges=list(lty = 1, alpha=0.2),
  factor_names = TRUE, 
  labels=list(
    font=fontsize
  )
  
))
fig_2d_gtex_gencode_venn

ggsave(plot=fig_2d_gtex_gencode_venn, filename="out/main_figures/Figure_2/Figure_2d.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 2e
df <- fread(gtex_v8_pr_curve)
df <- process_models_pr_curve(df, chosen_models)

fig_2e_gtex_v8_pr_curve <- plot_pr_curve(
  df,
  ylim=0.4,
  color_list
)

ggsave(plot=fig_2e_gtex_v8_pr_curve, filename="out/main_figures/Figure_2/Figure_2e.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)

# Fig. 2f
df <- fread(gtex_v8_boxplot)
df <- process_models_box_plot(df, chosen_models)

fig_2f_gtex_v8_boxplot <- plot_box_plot(
  df,
  color_list, 
  comparisons=list(
    c("MMSplice + SpliceMap", "MMSplice"),
    c("SpliceAI + SpliceMap", "SpliceAI")
  ),
  coord_flip=TRUE
)

ggsave(plot=fig_2f_gtex_v8_boxplot, filename="out/main_figures/Figure_2/Figure_2f.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)

# Arrange Fig. 2
fig_2 <- ggarrange(
  nrow = 3, widths = c(2,2), heights = c(5,3,3),
  ggarrange(ncol = 2, labels = c('a', 'b'), 
            fig_2a_sashimi_gene_unexpressed,
            fig_2b_sashimi_exon_unused
  ),
  ggarrange(ncol = 2, labels = c('c', 'd'),
            fig_2c_sashimi_weak_site,
            as.ggplot(plot(fig_2d_gtex_gencode_venn))
  ),
  ggarrange(ncol = 2, labels = c('e', 'f'),
            fig_2e_gtex_v8_pr_curve,
            fig_2f_gtex_v8_boxplot
  )
)
fig_2

ggsave(plot=fig_2, filename="out/main_figures/Figure_2/Fig_2.svg", height=unit(15, "cm"), width=unit(13, "cm"), dpi=450)
ggsave(plot=fig_2, filename="out/main_figures/Figure_2/Fig_2.pdf", height=unit(15, "cm"), width=unit(13, "cm"), dpi=450)
ggsave(plot=fig_2, filename="out/main_figures/Figure_2/Fig_2.png", height=unit(15, "cm"), width=unit(13, "cm"), dpi=450)
