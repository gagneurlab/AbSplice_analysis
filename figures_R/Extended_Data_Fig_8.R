binwidth = 0.1
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
logit <- function(x) log10( (x) / (1-x) )

dt_all <- fread('../data/gtex_v8/absplice_calibration.csv')
dt_all$target = dt_all$outlier
dt_all[, logits_AbSplice_DNA := logit(AbSplice_DNA)]

dt <- dt_all[, c(
  "AbSplice_DNA", 
  "logits_AbSplice_DNA", 
  "target" 
)]

dt[, bins := cut(logits_AbSplice_DNA, breaks=seq(min(logits_AbSplice_DNA)-1, max(logits_AbSplice_DNA), binwidth))]
dt[, mean_bins := mean(logits_AbSplice_DNA), by=bins]

props <- dt[, `:=`(prop=mean(target==1),
                   mean_bins= mean(logits_AbSplice_DNA),
                   size=.N)
            , by=bins]
dt[, odds := sum(target==1)/sum(target==0), by=bins]


# props <- sample_frac(props, 0.01)

# ============================================ log odds ===================================
p_log_odds <- (
  ggplot(props, aes(mean_bins, odds,
                    color=(!(odds==0 | is.infinite(odds)))))
  + geom_point()
  + geom_abline(slope=1)
  + ylab("Odds")
  + scale_y_log10(breaks=breaks, minor_breaks =minor_breaks)
  # + scale_y_log10()
  
  # + scale_color_manual(values=c("#999999","#000000"))
  + scale_color_manual(values=c("#000000"))
  
  + labs(
    x="mean bins of logit(AbSplice-DNA)",
  )
  + theme_cowplot(font_size = fontsize)
  # + background_grid(major = "xy", minor = "xy")
  + theme(legend.position = "none")
  + annotation_logticks(side="l")
)
p_log_odds
ggsave(plot=p_log_odds, filename="out/Extended_Data_Figures/Extended_Data_Figure_8/Extended_Data_Figure_8c.pdf", height=unit(4, "cm"), width=unit(12, "cm"), dpi=450)

# ============================================ histograms ===================================
hist_0 <- (
  ggplot(dt[target==0], aes(logits_AbSplice_DNA))
  + geom_histogram(binwidth = binwidth)
  + geom_vline(aes(xintercept = median(logits_AbSplice_DNA)),color="red", linetype="dashed", size=0.5)
  + labs(
    y="Number of non-outliers",
    x='logit(AbSplice-DNA)'
  )
  + theme_cowplot(font_size = fontsize)
  # + background_grid(major = "xy", minor = "xy")
)
ggsave(plot=hist_0, filename="out/Extended_Data_Figures/Extended_Data_Figure_8/Extended_Data_Figure_8a.pdf", height=unit(4, "cm"), width=unit(12, "cm"), dpi=450)


hist_1 <- (
  ggplot(dt[target==1], aes(logits_AbSplice_DNA))
  + geom_histogram(binwidth = binwidth)
  + geom_vline(aes(xintercept = median(logits_AbSplice_DNA)),color="red", linetype="dashed", size=0.5)
  + labs(
    y="Number of outliers",
    x='logit(AbSplice-DNA)'
  )
  + theme_cowplot(font_size = fontsize)
  # + background_grid(major = "xy", minor = "xy")
)
ggsave(plot=hist_1, filename="out/Extended_Data_Figures/Extended_Data_Figure_8/Extended_Data_Figure_8b.pdf", height=unit(4, "cm"), width=unit(12, "cm"), dpi=450)


# Extended Data Fig. 8
plot_calibration_AbSplice_DNA <- ggarrange(
  ncol = 1, nrow = 3, labels = c('a', 'b', 'c'),
  hist_0,
  hist_1,
  p_log_odds
)
ggsave(plot=plot_calibration_AbSplice_DNA, filename="out/Extended_Data_Figures/Extended_Data_Figure_8/Extended_Data_Figure_8.pdf", height=unit(12, "cm"), width=unit(12, "cm"), dpi=450)
ggsave(plot=plot_calibration_AbSplice_DNA, filename="out/Extended_Data_Figures/Extended_Data_Figure_8/Extended_Data_Figure_8.png", height=unit(12, "cm"), width=unit(12, "cm"), dpi=450)

