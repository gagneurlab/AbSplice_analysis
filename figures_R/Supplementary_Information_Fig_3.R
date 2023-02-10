df <- fread('../data/mito/venn_all_tissues.csv')

df$tissue_gtex <- revalue(df$tissue_gtex, tissue_mapping, warn_missing = FALSE)
tissue_list <- list(df[order(-Jaccard),]$tissue_gtex)
df$tissue_gtex <- factor(df$tissue_gtex, levels=tissue_list[[1]])

SI_Fig_3 <- (
  ggplot(df, aes(x=tissue_gtex, y=Jaccard))
  + geom_point()
  # + geom_point(data=g1, colour="red")
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    x='GTEx tissue',
    y="Jaccard of splice sites [Mito-disease fibroblasts vs GTEx tissue]",
    title=element_blank(),
  )
  + theme(axis.text.x = element_text(
    angle = 90, 
    hjust = 1
    ))
)
SI_Fig_3

ggsave(plot=SI_Fig_3, filename="out/Supplementary_Figures/Supplementary_Figure_3/Supplementary_Figure_3.pdf", height=unit(10, "cm"), width=unit(12, "cm"), dpi=450)
ggsave(plot=SI_Fig_3, filename="out/Supplementary_Figures/Supplementary_Figure_3/Supplementary_Figure_3.png", height=unit(10, "cm"), width=unit(12, "cm"), dpi=450)
