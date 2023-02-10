source("config.R")

# Fig. 4a
df <- fread('../data/mito/venn.csv')

splice_site_prokisch <- df[df$annotation_source == 'prokisch', splice_site]
splice_site_splicemap <- df[df$annotation_source == 'gtex', splice_site]
venn_diag <- list(
  'SpliceMap \nMito-disease \n(fibroblasts)' = splice_site_prokisch,
  'SpliceMap \nGTEx \n(fibroblasts)' = splice_site_splicemap
)

g_mito_venn_splice_sites = as.ggplot(plot(
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
ggsave(plot=g_mito_venn_splice_sites, filename="out/main_figures/Figure_4/Figure_4a.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 4b
df <- fread('../data/mito/ref_psi_corr.csv')

g_mito_splicemap_ref_psi_correlation <- (
  ggplot(
    df,
    aes(x=ref_psi_GTEx, y=ref_psi_mito)
    )
  + geom_hex(aes(fill = stat(log(count))))
  + scale_fill_gradientn(colors = brewer.pal(10,"Blues"))
  + labs(x='Ψ ref \n(GTEx fibroblasts)', y='Ψ ref \n(Mito-disease fibroblasts)', 
  )
  + theme_cowplot()
)
g_mito_splicemap_ref_psi_correlation
ggsave(plot=g_mito_splicemap_ref_psi_correlation, filename="out/main_figures/Figure_4/Figure_4b.svg", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)

# Fig. 4c
df <- fread('../data/mito/pr_jackknife.csv')

chosen_models <- c(
  "MMSplice",
  "SpliceAI",
  "AbSplice-DNA"
)
model_rename <- c(
  "MMSplice" = "MMSplice",
  "SpliceAI" = "SpliceAI",
  "AbSplice-DNA" = "AbSplice-DNA \n(SpliceMap \n(GTEx fibroblasts)"
)

color_list = absplice_model_colors[chosen_models]

df <- df[df$model %in% chosen_models]
df$model <- factor(df$model, levels=chosen_models)


df$model_names <- revalue(df$model, model_rename, warn_missing = FALSE)
df$model_names <- paste0(df$model_names, '\nn=', df$number_of_samples)
label_dict <- setNames(df$model_names, df$model)

g_mito_jackknife <- (
  ggplot(df, aes(x=model, y=mean))
  + scale_fill_manual(values=color_list)
  + scale_color_manual(values=color_list)
  + geom_bar(aes(
      fill=model
    ),
    stat = "identity", 
    width=0.5,
    position = "dodge2",
    size=0.3
  )
  + geom_errorbar(
    aes(
      ymin=mean - sem,
      ymax=mean + sem
    ),
    position = position_dodge(width = 0.9),
    width = 0.2
  )
  + scale_x_discrete(labels = function(x) {
    x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
  })
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    x=element_blank(),
    y="auPRC",
    title=element_blank(),
    color='SpliceMap',
    fill=element_blank()
  )
  + guides(
    fill = "none"
  )
  + scale_x_discrete(labels = label_dict)
)
g_mito_jackknife
ggsave(plot=g_mito_jackknife, filename="out/main_figures/Figure_4/Figure_4c.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 4d
df <- fread('../data/als/known_als_genes.csv')

chosen_models <- c(
  'MMSplice',
  'SpliceAI',
  'AbSplice-DNA \n(SpliceMap \nGTEx brain)',
  'AbSplice-DNA \n(SpliceMap ALS \nmotor neuron)',
  'random'
)
color_list = absplice_model_colors[chosen_models]
df <- df[df$model %in% chosen_models]
df$model <- factor(df$model, levels=chosen_models)

df$model_name <- paste0(df$model, '\nn=', df$above_cutoff_all)
df_mapping <- df %>% select(model, model_name) %>% unique()
label_dict <- setNames(df_mapping$model_name, df_mapping$model)

df$significance_label_pval <- ifelse(df$pval <= 0.05, paste0('P=', df$significance_label_pval), 'N.S.')

g_als_enrichment <- (
  ggplot(df[df$model != 'random'])
  + aes_string(
    x='model', 
    y='enrichment', 
    fill='model', 
    # label="significance_label"
    label="significance_label_pval"
    )
  + geom_text(
    nudge_y = +0.1,
    # nudge_x = -1,
    size=4, color="black"
  )
  + geom_bar(
    stat = "identity", 
    width=0.5,
    position = "dodge2"
  )
  + scale_fill_manual(values=color_list)
  + theme_cowplot(font_size = fontsize)
  + scale_x_discrete(labels = function(x) {
    x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
  })
  # + background_grid()
  + labs(
    x=element_blank(),
    y="Enrichment \nin ALS genes",
    title=element_blank()
  )
  + guides(
    fill = "none"
  )
  + geom_hline(yintercept=1, linetype="dashed")
  + scale_x_discrete(labels = label_dict)
)
g_als_enrichment
ggsave(plot=g_als_enrichment, filename="out/main_figures/Figure_4/Figure_4d.pdf", height=unit(4, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 4e
df <- fread('../data/als/proteomics.csv')

chosen_models <- c(
  'All genes',
  'MMSplice',
  'SpliceAI',
  'AbSplice-DNA \n(GTEx brain)',
  'AbSplice-DNA \n(motor neuron)'
)
color_list = absplice_model_colors[chosen_models]
df <- df[df$model %in% chosen_models]
df$model <- factor(df$model, levels=chosen_models)

chosen_cutoffs <- c(
  "high"
)

df <- df[cutoff_category %in% chosen_cutoffs] 

library(binom)
# binom.confint(x = c(18), n = 58, tol = 1e-8, method="prop.test")
df_test <- binom.confint(
  x = df$num_protein_low, 
  n = df$num, 
  tol = 1e-8, 
  # method="prop.test"
  method="exact"
)

df_plot <- cbind(df, df_test)

df_plot$model_name <- paste0(df_plot$model, '\nn=', df$num)
df_mapping <- df_plot %>% select(model, model_name) %>% unique()
label_dict <- setNames(df_mapping$model_name, df_mapping$model)

g_als_proteomics_true_preds <- (
  ggplot(df_plot, aes(
    x=model,
    y=proportion_true
  ))
  + scale_fill_manual(values=color_list)
  + scale_color_manual(values=color_list)
  + geom_bar(aes(
    fill=model
  ),
  stat = "identity", 
  width=0.5,
  position = "dodge2",
  size=0.3
  )
  + geom_errorbar(
    aes(
      ymin=lower,
      ymax=upper
    ),
    position = position_dodge(width = 0.9),
    width = 0.2
  )
  + scale_x_discrete(labels = function(x) {
    x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
  })
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    x=element_blank(),
    y="Proportion high impact predictions \nvalidated via proteomics",
    title=element_blank(),
    color='SpliceMap',
    fill=element_blank()
  )
  + guides(
    fill = "none"
  )
  + scale_x_discrete(labels = label_dict)
  + coord_flip()
)
g_als_proteomics_true_preds
ggsave(plot=g_als_proteomics_true_preds, filename="out/main_figures/Figure_4/Figure_4e.pdf", height=unit(6, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 4f
df <- fread('../data/LOEUF/matched_recall.csv')
chosen_models <- c(
  'AbSplice-DNA', 
  'SpliceAI',
  'MMSplice'
)
color_list = absplice_model_colors[chosen_models]
df <- df[df$model %in% chosen_models]
df$model <- factor(df$model, levels=chosen_models)

df$significance_label_pval <- ifelse(df$pval_spliceai_label <= 0.05, paste0('P=', df$pval_spliceai_label), 'N.S.')

g_LOEUF_norm <- (ggplot(df)
                 + aes_string(
                   x='LOEUF_decile', 
                   y='depletion_high_impact_variants',
                   fill='model', 
                   # label="significance_label_spliceai"
                   label="significance_label_pval"
                 )
                 + geom_text(
                   # nudge_x = -0.2,
                   nudge_y = 0.2,
                   size=3, 
                   color="black")
                 + geom_bar(stat='identity', position='dodge')
                 + labs(x='LOEUF decile', y='Depletion \nof high impact SNVs')
                 + scale_fill_manual(values=color_list)
                 + theme_cowplot(font_size=fontsize)
                 # + background_grid()
                 + scale_x_continuous(
                   breaks = seq(0, 10, 1),
                 )
                 + theme(
                   legend.title = element_blank(),
                   legend.position= c(0.7, 0.8)
                 )
)
g_LOEUF_norm
ggsave(plot=g_LOEUF_norm, filename="out/main_figures/Figure_4/Figure_4f.pdf", height=unit(6, "cm"), width=unit(6, "cm"), dpi=450)


# Fig. 4
fig_4 <- ggarrange(
  nrow = 3, align = 'v',
  ggarrange(ncol=2, widths = c(1, 1), labels = c('a', 'b'),
            g_mito_venn_splice_sites,
            g_mito_splicemap_ref_psi_correlation
  ),
  ggarrange(ncol=2, widths = c(1, 1), labels = c('c', 'd'),
            g_mito_jackknife,
            g_als_enrichment
  ),
  ggarrange(ncol=2, labels = c('e', 'f'),
            g_als_proteomics_true_preds,
            g_LOEUF_norm
  )
)
fig_4

ggsave(plot=fig_4, filename="out/main_figures/Figure_4/Figure_4.svg", height=unit(12, "cm"), width=unit(12, "cm"), dpi=dpi)
