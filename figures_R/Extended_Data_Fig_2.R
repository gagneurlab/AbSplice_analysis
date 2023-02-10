source("config.R")

# Extended Data Figure 2a
df <- fread('../data/gtex_v8/enrichment_outlier_rare_variant_distance.csv')
df$num_tissue_cutoff <- as.character(df$num_tissue_cutoff)

Ext_Data_Fig_2a <- (
  ggplot(df, aes(x=abs_Distance, y=local_enrichment, color=num_tissue_cutoff))
  + geom_line()
  + xlim(0, 1000)
  + theme_cowplot(font_size = fontsize)
  # + background_grid(major = "xy", minor = "xy")
  + labs(
    color='Number of tissues',
    x='Distance between variant and outlier [bp]',
    y='Enrichment of \nreplicated splicing outliers'
  )
  + theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )
  + geom_hline(
    yintercept=1,
    linetype="dashed"
  ) 
)
Ext_Data_Fig_2a
ggsave(plot=Ext_Data_Fig_2a, filename="out/Extended_Data_Figures/Extended_Data_Figure_2/Extended_Data_Figure_2a.pdf", height=unit(6, "cm"), width=unit(8, "cm"), dpi=450)

# Extended Data Figure 2b
df <- fread('../data/gtex_v8/fraser_singleton_reduction.csv')

Ext_Data_Fig_2b <- (
  ggplot(df, aes(
    x=outlier_with_variant, 
    y=proportion, 
    fill=outlier_with_variant
  ))
  + geom_boxplot()
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    y="Proportion of \nreplicated/non-replicated \noutliers",
    x=element_blank(),
    fill='Outlier with \nrare variant'
  )
  # + coord_flip()
  # + theme(
  #   axis.text.y=element_blank(),
  #   axis.ticks.y=element_blank()
  #   )
  # + theme(
  #   legend.position = c(0.95, 0.95),
  #   legend.justification = c("right", "top")
  # )
)
Ext_Data_Fig_2b
ggsave(plot=Ext_Data_Fig_2b, filename="out/Extended_Data_Figures/Extended_Data_Figure_2/Extended_Data_Figure_2b.pdf", height=unit(5, "cm"), width=unit(5, "cm"), dpi=450)

# Extended Data Figure 2c
df <- fread('../data/gtex_v8/fraser_singleton_proportion.csv')
df$tissue <- revalue(df$tissue, tissue_mapping, warn_missing = FALSE)

df_ <- df[df$outlier_with_variant == 'FALSE']
tissue_list <- list(df_[order(-percentage_of_singleton),]$tissue)
df$tissue <- factor(df$tissue, levels=tissue_list[[1]])

Ext_Data_Fig_2c <- (
  ggplot(df, aes(y=tissue, x=percentage_of_singleton))
  + geom_bar(aes(
    fill=outlier_with_variant
  ),
  stat = "identity", 
  position = "dodge2",
  size=0.5
  )
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    y=element_blank(),
    x="Percentage of singletons",
    title=element_blank(),
    fill='Outlier with \nrare variant'
  )
  + coord_flip()
  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # + theme(
  #   legend.position = c(0.95, 0.95),
  #   legend.justification = c("right", "top")
  # )
)
Ext_Data_Fig_2c
ggsave(plot=Ext_Data_Fig_2c, filename="out/Extended_Data_Figures/Extended_Data_Figure_2/Extended_Data_Figure_2c.pdf", height=unit(6, "cm"), width=unit(12, "cm"), dpi=450)

# Extended Data Figure 2
Ext_Data_Fig_2 <- ggarrange(nrow=2, heights=c(2,5),
                            ggarrange(
                              ncol = 2, widths = c(2, 1), labels = c('a', 'b'), align='h',
                              Ext_Data_Fig_2a,
                              Ext_Data_Fig_2b
                            ),
  ggarrange(Ext_Data_Fig_2c)
)
Ext_Data_Fig_2

ggsave(plot=Ext_Data_Fig_2, filename="out/Extended_Data_Figures/Extended_Data_Figure_2/Extended_Data_Figure_2.svg", height=unit(12, "cm"), width=unit(12, "cm"), dpi=450)
