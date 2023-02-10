
# SI Fig. 5a
df <- fread('../data/als/venn_all_tissues.csv')

df$tissue_gtex <- revalue(df$tissue_gtex, tissue_mapping, warn_missing = FALSE)
tissue_list <- list(df[order(-Jaccard),]$tissue_gtex)
df$tissue_gtex <- factor(df$tissue_gtex, levels=tissue_list[[1]])

g_als_splicemap_splice_sites_all_gtex_tissues <- (
  ggplot(df, aes(x=tissue_gtex, y=Jaccard, color=brain_tissue))
  + geom_point()
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    x='GTEx tissue',
    y="Jaccard of splice sites \n[ALS motor neurons vs GTEx tissue]",
    title=element_blank(),
    color='Brain tissue'
  )
  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  + theme(
    legend.position = c(1, 1),
    legend.justification = c("right", "top")
  )
)
g_als_splicemap_splice_sites_all_gtex_tissues
ggsave(plot=g_als_splicemap_splice_sites_all_gtex_tissues, filename="out/Supplementary_Figures/Supplementary_Figure_5/Supplementary_Figure_5a.pdf", height=unit(8, "cm"), width=unit(12, "cm"), dpi=450)

# SI Fig. 5b
df <- fread('../data/als/venn.csv')

splice_site_als <- df[df$annotation_source == 'als', splice_site]
splice_site_splicemap <- df[df$annotation_source == 'gtex', splice_site]
venn_diag <- list(
  'SpliceMap \nALS \n(motor neurons)' = splice_site_als,
  'SpliceMap \nGTEx \n(brain)' = splice_site_splicemap
)

g_als_venn_splice_sites = as.ggplot(plot(
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
g_als_venn_splice_sites
ggsave(plot=g_als_venn_splice_sites, filename="out/Supplementary_Figures/Supplementary_Figure_5/Supplementary_Figure_5b.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)

# SI Fig. 5c
df <- fread('../data/als/ref_psi_corr.csv')

g_als_splicemap_ref_psi_correlation <- (
  ggplot(
    df,
    aes(x=ref_psi_GTEx, y=ref_psi_ALS)
  )
  + geom_hex(aes(fill = stat(log(count))))
  + scale_fill_gradientn(colors = brewer.pal(10,"Blues"))
  + labs(x='Ψ ref \n(GTEx Brain)', y='Ψ ref \n(ALS motor neurons)', 
  )
  + theme_cowplot()
)
g_als_splicemap_ref_psi_correlation
ggsave(plot=g_als_splicemap_ref_psi_correlation, filename="out/Supplementary_Figures/Supplementary_Figure_5/Supplementary_Figure_5c.svg", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)

# SI Fig. 5
SI_fig_5 <- ggarrange(
  nrow = 2, heights=c(2.5,1),
  ggarrange(labels = c('a'),
            g_als_splicemap_splice_sites_all_gtex_tissues
  ),
  ggarrange(ncol=2, labels = c('b', 'c'),
            g_als_venn_splice_sites,
            g_als_splicemap_ref_psi_correlation
  )
)
SI_fig_5
ggsave(plot=SI_fig_5, filename="out/Supplementary_Figures/Supplementary_Figure_5/Supplementary_Figure_5.svg", height=unit(12, "cm"), width=unit(12, "cm"), dpi=450)



