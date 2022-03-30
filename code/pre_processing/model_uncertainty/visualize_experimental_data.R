rm(list=ls(all=TRUE))

library(rmatio)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(dplyr)
library(reshape2)
library(tibble)
library(latex2exp)
# library(scatterPlotMatrix)
library("GGally")


# Set working directory to this file, works only with RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# =================================================================
# Input and options
# =================================================================
data_dir = file.path('..', '..', '..', 'data')
fig_dir = 'results/'

data_file_name = 'beam_p_adj.mat'

dpi = 400

# =================================================================
# Pre-process
# =================================================================
# read data
data_path = file.path(data_dir, data_file_name)
exp_data = read.mat(data_path)

V_Rexp      = exp_data$Vu

# in [MPa], according to par. 7.2.3 (for our purpose, we use the fc;test value; otherwise: fck = fcmean - 8)
fc          = exp_data$fc         

# in [mm], width of the beam
b           = exp_data$b        

# in [mm], effective height
d           = exp_data$d          

rho         = exp_data$rho/100

# note that this done because dg is fixed in the calibration!!
dg          = 16*rep(1, length(rho))

# in [mm2], area of tensile reinforcement in considered section
Asl         = rho*b*d       

exp_data_df = tibble(V_Rexp=V_Rexp, fc=fc, b=b, d=d, rho=rho)

column_labels = c(
  expression(italic("V")*""["exp"] * " [kN]"),
  expression(italic("f")*""["c"] * " [MPa]"),
  expression(italic("b")*""["w"]* " [mm]"),
  expression(italic("d")* " [mm]"),
  expression(italic(rho)*""["l"]* "")
)

column_labels2 = c()

for (column_label in column_labels) {
  dp = deparse(column_label)
  dp = gsub(pattern="^expression\\(", replacement="", x=dp)
  dp = gsub(pattern="\\)$", replacement="", x=dp)
  column_labels2 = append(column_labels2, dp)
}

# =================================================================
# Visualize
# =================================================================

# ................................................................
# Histograms with ggplot2
# ................................................................
df = melt(exp_data_df)
levels(df$variable) <- column_labels

g = ggplot(df, aes(x=value))
g = g + geom_histogram()
g = g + facet_grid(. ~ variable, scales="free", labeller = label_parsed)

plot(g)

# ................................................................
# Scatter plot with scatterPlotMatrix
# ................................................................
# Pairwise plots with histograms - interactive
# scatterPlotMatrix(exp_data_df, corrPlotType = "Text")


# ................................................................
# Scatter plot with ggpairs
# ................................................................
gp = ggpairs(
  exp_data_df,
  diag=list(continuous="barDiag"),
  upper=list(continuous = wrap(
    "points", alpha = 0.3, color="#565656", shape=16)),
  lower = list(continuous = wrap(ggally_cor, stars = F)),
  columnLabels = column_labels2,
  labeller = "label_parsed"
)

fname = file.path(fig_dir, "scatter_plot_matrix_exp_data.png")
ggsave(fname, plot=gp, dpi=dpi, units=c("mm"), width=200, height=200)

