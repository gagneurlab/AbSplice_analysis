source("config.R")

# Extended Data Fig. 7a
linspace <- function(x1, x2, n=100) {
  stopifnot(is.numeric(x1), is.numeric(x2), length(x1)==1, length(x2)==1)
  n <- floor(n)
  if (n <= 1) x2
  else seq(x1, x2, length.out=n)
}

sigmoid = function(x) {
  1 / (1 + exp(-x))
}

delta_logit_psi = 1.5
y1 = 0.7
y2 = 0.1
k = 1
a = -(log(1/y1 - 1))
b = -(log(1/y2 - 1))
x1 = a/k
x2 = b/k
x1_var = x1 - delta_logit_psi
x2_var = x2 - delta_logit_psi
y1_var = 1/(1 + exp(-k*x1_var))
y2_var = 1/(1 + exp(-k*x2_var))

x <- linspace(-10, 10, 1000)
y = 1/(1 + exp(-k*x))
df <- data.frame('x' = x, 'y' =y)

point_size = 3

segments_dt = data.table(
  x = c(x1, x1_var, x2, x2_var),
  y = c(y1, y1, y2, y2),
  xend = c(x1_var, x1_var, x2_var, x2_var),
  yend = c(y1, y1_var, y2, y2_var),
  "mutation effect" = c("Δlogit(Ψ)", "ΔΨ", "Δlogit(Ψ)", "ΔΨ")
)

reference="reference"
alternative="alternative"
points_dt = data.table(
  x=c(x1, x2, x1_var, x2_var),
  y=c(y1, y2, y1_var, y2_var),
  `sample type`=c(reference, reference, alternative, alternative),
  tissue=c("lung", "brain", "lung", "brain")
)
points_dt
+ geom_point(data=data.frame(x1,y1), aes(x1,y1), colour='#ff0000', size=point_size)
+ geom_point(data=data.frame(x2,y2), aes(x2,y2), colour='#2204ff', size=point_size)
+ geom_point(data=data.frame(x1_var,y1_var), aes(x1_var, y1_var), colour='#8a1308', size=point_size)
+ geom_point(data=data.frame(x2_var,y2_var), aes(x2_var, y2_var), colour='#15058b', size=point_size)

show_point_shapes(16)

points_dt$`sample type` <- factor(points_dt$`sample type`, levels = c("reference", "alternative"))

Ext_Data_Fig_7a <- (
  ggplot(data=df, aes(x, y)) 
  + geom_line()
  + geom_point(data=points_dt, aes_string(
    x="x", 
    y="y",
    color="`sample type`"
  ), size=4)
  + theme_cowplot(font_size = fontsize)
  # + background_grid()
  + labs(
    x = 'Δ logit(Ψ)', 
    y = 'Ψ', 
    color=''
  )
  + scale_color_manual(values=c('black', 'red'))
  + geom_segment(
    data=segments_dt, 
    aes_string(
      x="x", 
      y="y", 
      xend="xend", 
      yend="yend"
    ), 
    lineend="round", 
    arrow = arrow(type="closed", length=unit(0.2, "cm"))
  )
  + theme(
    legend.position = c(0.1, 0.9)
  )
  
)
ggsave(plot=Ext_Data_Fig_7a, filename="out/Extended_Data_Figures/Extended_Data_Figure_7/Extended_Data_Figure_7a.svg", height=unit(4, "cm"), width=unit(4, "cm"), dpi=450)

# Extended Data Fig. 7b
df <- fread('../data/gtex_v8/ref_psi_distribution.csv')

Ext_Data_Fig_7b <- (
  ggplot(df, aes(x=ref_psi))
  + geom_histogram(bins = 30)
  + theme_cowplot(font_size = fontsize)
  # + background_grid(major = "xy", minor = "xy")
  + facet_wrap('psi')
  + theme(
    strip.background = element_rect(colour="white", fill="white")
  )
  + labs(
    x='Ψ_ref',
    y='Introns'
  )
  + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
)
ggsave(plot=Ext_Data_Fig_7b, filename="out/Extended_Data_Figures/Extended_Data_Figure_7/Extended_Data_Figure_7b.svg", height=unit(4, "cm"), width=unit(4, "cm"), dpi=450)

# Extended Data Fig. 7c
df <- fread('../data/gtex_v8/ref_psi_distribution_cumulative.csv')

Ext_Data_Fig_7c <- (
  ggplot(df, aes(x=ref_psi_diff, color=psi))
  + scale_color_manual(values=c('royalblue4', '#d62728'))
  + stat_ecdf(geom = "step")
  + theme_cowplot(font_size = fontsize)
  + labs(
    x='max(Ψ_ref) - min(Ψ_ref)',
    y='CDF',
    color=''
  )
  # + background_grid(major = "xy", minor = "xy")
  + scale_y_continuous(limits = c(0.8, 1.01))
  + theme(
    legend.position = c(0.5, 0.5)
  )
)
ggsave(plot=Ext_Data_Fig_7c, filename="out/Extended_Data_Figures/Extended_Data_Figure_7/Extended_Data_Figure_7c.svg", height=unit(4, "cm"), width=unit(4, "cm"), dpi=450)

# Extended Data Fig. 7d
img_ref_psi_heatmap <- readPNG('../data/gtex_v8/ref_psi_heatmap.png')
Ext_Data_Fig_7d <- as.raster(img_ref_psi_heatmap)
Ext_Data_Fig_7d <- rasterGrob(Ext_Data_Fig_7d, interpolate = FALSE)
ggsave(plot=Ext_Data_Fig_7d, filename="out/Extended_Data_Figures/Extended_Data_Figure_7/Extended_Data_Figure_7d.pdf", height=unit(8, "cm"), width=unit(12, "cm"), dpi=450)

# Extended Data Fig. 7
Ext_Data_Fig_7 <- ggarrange(
  nrow = 2, heights=c(1,2),
  ggarrange(ncol = 3, labels = c('a', 'b', 'c'),
            Ext_Data_Fig_7a,
            Ext_Data_Fig_7b,
            Ext_Data_Fig_7c
  ),
  ggarrange(ncol = 1, labels = c('d'),
            Ext_Data_Fig_7d
  )
)
ggsave(plot=Ext_Data_Fig_7, filename="out/Extended_Data_Figures/Extended_Data_Figure_7/Extended_Data_Figure_7.svg", height=unit(12, "cm"), width=unit(12, "cm"), dpi=450)
