# -*- coding: utf-8 -*-
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(cowplot)
  library(ggthemes)
  library(ggpubr)
  library(ggforce)
  library(ggrepel)
  library(gplots)
  library(png)
  library(ggvenn)
  library(eulerr)
  library(ggnewscale)
  require(ggbeeswarm)
  library(plyr)
  library(dplyr)
  library(RColorBrewer)
  library(heatmaply)
  library(ggplotify)
  library(UpSetR)
  library(Rmisc)
  library(ggupset)
  library(tidyverse, warn.conflicts = FALSE)
  library(stringr)
})

fontsize <- 14
linesize <- 1
dpi=450

gtex_v8_pr_curve = '../data/gtex_v8/pr_curve_all_tissues.csv'
gtex_v8_boxplot = '../data/gtex_v8/boxplot_across_tissues.csv'

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
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top")
    )
    + labs(x='Recall', y='Precision', color='', title=title)
  )
  return(g)
}


plot_box_plot <- function(
  df, color_list, title="Across tissues",
  comparisons = list(
    c("AbSplice-DNA", "MMSplice + SpliceMap + Ψ_ref")
  ),
  coord_flip=FALSE,
  jitter=FALSE,
  x_lim=max(df$`Average Precision Score`) * 1.3
  ) {
  g <- (
    ggplot(df, aes(x=`Average Precision Score`, y=model, fill=model))
    + scale_fill_manual(values=color_list)
    + geom_boxplot()
    + xlim(0, x_lim)
    + scale_y_discrete(labels = function(x) {
      x %>% str_replace_all(" \\+", "\n+") %>% str_replace_all(" \\(", "\n(")
    })
    + theme_cowplot(font_size = fontsize)
    # + background_grid()
    + theme(legend.title = element_blank())
    + theme(legend.position = "none")
    + labs(
      x="auPRC",
      # x="area under the precision-recall curve",
      y=element_blank(),
      title=title
    )
    + {if (coord_flip == TRUE)coord_flip()}
    + stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      paired = TRUE,
      method.args = list(alternative = "greater")
    )
    + {if (jitter == TRUE)geom_jitter(color="black", size=0.4, alpha=0.9)}
  )
  return(g)
}

set_start_precision_zero <- function(df, chosen_models) {
  for (model_temp in chosen_models) {
    df_model = df[df[,df$model==model_temp]]
    df$precision[df$model==model_temp][nrow(df_model)] <- 0
  }
  return(df)
}

process_models_pr_curve <- function(df, chosen_models) {
    df <- df[df$model %in% chosen_models]
    df$model <- factor(df$model, levels = chosen_models,ordered = TRUE)
    df <- set_start_precision_zero(df, chosen_models)
    return(df)
}

process_models_box_plot <- function(df, chosen_models) {
    df <- df[df$model %in% chosen_models]
    df$model <- factor(df$model, levels = chosen_models,ordered = TRUE)
    return(df)
}

absplice_model_colors = c(
  "CADD-Splice" = '#20b22a',
  "SQUIRLS" = "#d96c0d",
  "MTSplice" = "#28666c",
  "SpliceAI" = "#adacac",
  
  "SpliceAI + SpliceMap" = "#797677",
  "SpliceAI + SpliceMap \n(SpliceMap \nGTEx brain)" = "#797677",
  "SpliceAI + SpliceMap \n(SpliceMap \nALS motor neuron)" = "#797677",
  
  "SpliceAI + SpliceMap + Ψ_ref" = "black",
  "SpliceAI + SpliceMap + Ψ_ref \n(SpliceMap \nGTEx brain)" = "black",     
  "SpliceAI + SpliceMap + Ψ_ref \n(SpliceMap \nALS motor neuron)" = "black",

  "MMSplice" = "paleturquoise3",# "#00BEFF", ##A6DEFF",
  
  "MMSplice + SpliceMap" = "royalblue",#"#00BEFF",
  "MMSplice + SpliceMap \n(SpliceMap \nGTEx brain)" = "royalblue",#"#00BEFF",
  "MMSplice + SpliceMap \n(SpliceMap \nALS motor neuron)" = "royalblue",#"#00BEFF",
  
  "MMSplice + SpliceMap + Ψ_ref" = "royalblue4",#"#006199",
  "MMSplice + SpliceMap + Ψ_ref \n(SpliceMap \nGTEx brain)"  = "royalblue4",#"#006199",
  "MMSplice + SpliceMap + Ψ_ref \n(SpliceMap \nALS motor neuron)"  = "royalblue4",#"#006199",
  
  "AbSplice-DNA" = "#d62728",
  "AbSplice-DNA \n(SpliceMap \nGTEx brain)" = "#d62728",
  "AbSplice-DNA \n(GTEx brain)" = "#d62728",
  "AbSplice-DNA \n(SpliceMap \nALS motor neuron)" = "#d62728",
  "AbSplice-DNA \n(SpliceMap ALS \nmotor neuron)" = "#d62728",
  "AbSplice-DNA \n(motor neuron)" = "#d62728",
  "AbSplice-DNA \n(SpliceMap \nGTEx fibroblasts)" = "#d62728",
  
  "AbSplice-DNA (GAM)" = "#d62728",
  "AbSplice-DNA (random forest)" = "#619CFF",
  "AbSplice-DNA (logisitic regression)" = "#00BA38",
  
  "random" = "black",
  "All genes" = "black",
  
  'CAT p-value (fibroblasts)' = "#ffd33f",
  "CAT p-value (lymphocytes)" = "#ffd33f",
  "CAT p-value (blood)" = "#ffd33f",
  'CAT p-value' = "#ffd33f",
  
  "AbSplice-RNA (fibroblasts)" = "#b28c4d",
  "AbSplice-RNA (lymphocytes)" = "#b28c4d",
  "AbSplice-RNA (blood)" = "#b28c4d",
  "AbSplice-RNA" = "#b28c4d",
  "AbSplice-RNA (all CATs)" = "brown"
)


tissue_mapping = c(
  'Adipose_Subcutaneous' = "Adipose Subcutaneous",
  'Adipose_Visceral_Omentum' = "Adipose Visceral Omentum",
  'Adrenal_Gland' = "Adrenal Gland",   
  'Artery_Aorta' ="Artery Aorta",
  'Artery_Coronary' = "Artery Coronary", 
  'Artery_Tibial' = "Artery Tibial",
  'Brain_Amygdala' = 'Brain Amygdala',
  'Brain_Anterior_cingulate_cortex_BA24' = 'Brain Ant. cing. cortex BA24',
  'Brain_Caudate_basal_ganglia' = 'Brain Caudate basal ganglia',
  'Brain_Cerebellar_Hemisphere' = 'Brain Cerebellar Hemisphere',
  'Brain_Cerebellum' = 'Brain Cerebellum',
  'Brain_Cortex' = 'Brain Cortex',
  'Brain_Frontal_Cortex_BA9' = 'Brain Frontal Cortex BA9',
  'Brain_Hippocampus' = 'Brain Hippocampus',
  'Brain_Hypothalamus' = 'Brain Hypothalamus',
  'Brain_Nucleus_accumbens_basal_ganglia' = 'Brain Nuc. basal ganglia',
  'Brain_Putamen_basal_ganglia' = 'Brain Putamen basal ganglia',
  'Brain_Spinal_cord_cervical_c_1' = 'Brain Spinal cord cervical c1',
  'Brain_Substantia_nigra' = 'Brain Substantia nigra',
  'Breast_Mammary_Tissue' = "Breast",
  'Colon_Sigmoid' = 'Colon Sigmoid',
  'Colon_Transverse' = 'Colon Transverse',
  'Esophagus_Gastroesophageal_Junction' = 'Esophagus GJ',
  'Esophagus_Mucosa' = 'Esophagus Mucosa',
  'Esophagus_Muscularis' = 'Esophagus Muscularis',
  'Heart_Atrial_Appendage' = 'Heart Atrial Appendage', 
  'Heart_Left_Ventricle' = 'Heart Left Ventricle',
  'Kidney_Cortex' = 'Kidney Cortex',
  'Liver' = 'Liver',
  'Lung' = 'Lung',
  'Minor_Salivary_Gland' = 'Minor Salivary Gland',
  'Muscle_Skeletal' = 'Muscle Skeletal',
  'Nerve_Tibial' = 'Nerve Tibial',
  'Ovary' = 'Ovary',
  'Pancreas' = 'Pancreas',
  'Pituitary' = 'Pituitary',
  'Prostate' = 'Prostate',
  'Skin_Not_Sun_Exposed_Suprapubic' = 'Skin No Sun',
  'Skin_Sun_Exposed_Lower_leg' = 'Skin Sun',
  'Small_Intestine_Terminal_Ileum' = 'Small Intestine',
  'Spleen' = 'Spleen',
  'Stomach' = 'Stomach',
  'Testis' = 'Testis',
  'Thyroid' = 'Thyroid',
  'Uterus' = 'Uterus',
  'Vagina' = 'Vagina',
  'allBrainTissues' = 'Brain',
  'Cells_Cultured_fibroblasts' = 'Fibroblasts',
  'Cells_EBV_transformed_lymphocytes' = 'Lymphocytes',
  'Whole_Blood' = 'Blood'
)

