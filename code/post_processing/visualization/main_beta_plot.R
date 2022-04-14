rm(list=ls(all=TRUE))

library(readr)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(dplyr)

# Set working directory to this file, works only with RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# =================================================================
# INPUT & OPTIONS
# =================================================================

fig_dir                 = 'results/'

ID                      = '2018-Jul-19_15.08.36'

beta_t                  = 4.7

# =================================================================
# PRE-PROCESS
# =================================================================

fpath                   = paste('input/beta_calibration_run_', ID, '.csv', sep = '')
df_beta                 = read_csv(fpath)

df_beta$load_comb       = as.factor(df_beta$load_comb)
df_beta$f_cck           = as.factor(df_beta$f_cck)
df_beta$rho             = as.factor(df_beta$rho)
df_beta$d               = as.factor(df_beta$d)
df_beta$weight_logic    = df_beta$weight > 0
# df_beta$chi2            = as.factor(df_beta$chi2)


# =================================================================
# VISUALIZE
# =================================================================

# -----------------------------------------------------------------
# BETA - CHI1
# -----------------------------------------------------------------

f_cck_ii                = levels(df_beta$f_cck)[2]
d_ii                    = levels(df_beta$d)[2]

df                      = df_beta %>% filter(f_cck == f_cck_ii & d == d_ii)
# df = df_beta

# for overlay plot - indicate weight range (not the best solution)
dfw                     = df
idx                     = df$weight_logic
dfw                     = dfw[!idx,]

g = ggplot(df, aes(x = chi1, y = beta, color = as.factor(chi2)))
g = g + geom_point(mapping = aes(size = weight), shape = 16, alpha = 0.7)
g = g + geom_path()
g = g + geom_hline(yintercept = beta_t)
g = g + facet_grid(rho ~ load_comb, labeller = label_both)
g = g + geom_point(data = dfw, color = 'white', fill = 'white', size = 1, stroke = 0)
g = g + scale_color_discrete(name = 'chi2')
g = g + scale_size_continuous(range = c(1.5,4))

g = g + ggtitle(paste('f_cck = ', f_cck_ii, '; d = ', d_ii, sep = ''))

g = g + theme_bw()
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g = g + xlim(c(0,1))


print(g)

# SAVE
fig_name  = paste('chi1_beta_', ID, '.png')
fig_path = file.path(fig_dir, fig_name)
ggsave(fig_path,
       width = 25, height = 15, units = 'cm',
       dpi = 600)

# -----------------------------------------------------------------
# BETA - CHI2
# -----------------------------------------------------------------


# for overlay plot - indicate weight range (not the best solution)
dfw                     = df
idx                     = df$weight_logic
dfw                     = dfw[!idx,]

g = ggplot(df, aes(x = chi2, y = beta, color = as.factor(chi1)))
g = g + geom_point(mapping = aes(size = weight), shape = 16, alpha = 0.7)
g = g + geom_path()
g = g + geom_hline(yintercept = beta_t)
g = g + facet_grid(rho ~ load_comb, labeller = label_both)
g = g + geom_point(data = dfw, color = 'white', fill = 'white', size = 1, stroke = 0)
g = g + scale_color_discrete(name = 'chi1')
g = g + scale_size_continuous(range = c(1.5,4))

g = g + ggtitle(paste('f_cck = ', f_cck_ii, '; d = ', d_ii, sep = ''))

g = g + theme_bw()
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g = g + xlim(c(0,1))


print(g)

# SAVE
fig_name  = paste('chi2_beta_', ID, '.png')
fig_path = file.path(fig_dir, fig_name)
ggsave(fig_path,
       width = 25, height = 15, units = 'cm',
       dpi = 600)