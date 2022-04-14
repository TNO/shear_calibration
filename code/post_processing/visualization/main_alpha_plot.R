rm(list=ls(all=TRUE))

library(readr)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(dplyr)
library(tibble)
library(latex2exp)

# Set working directory to this file, works only with RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# =================================================================
# INPUT & OPTIONS
# =================================================================

fig_dir                 = 'results/'

ID                      = '2018-Jul-19_15.08.36'

# =================================================================
# PRE-PROCESS
# =================================================================

fpath                   = paste('input/alpha2_calibration_run_', ID, '.csv', sep = '')
df_alpha                = read_csv(fpath)

# df_alpha = as.tibble(df_alpha)

#df_alpha$alpha_labels    = as.factor(df_alpha$alpha_labels)
df_alpha$load_comb_r     = as.factor(df_alpha$load_comb_r)
df_alpha$f_cck_r         = as.factor(df_alpha$f_cck_r)
df_alpha$rho_r           = as.factor(df_alpha$rho_r)
df_alpha$d_r             = as.factor(df_alpha$d_r)
df_alpha$weight_logic_r  = df_alpha$weight_r > 0
# df_alpha$chi1_r          = as.factor(df_alpha$chi1_r)
# df_alpha$chi2_r          = as.factor(df_alpha$chi2_r)


# =================================================================
# VISUALIZE
# =================================================================

# -----------------------------------------------------------------
# ALPHA'S
# -----------------------------------------------------------------

# f_cck_ii                = levels(df_alpha$f_cck_r)[4]
# d_ii                    = levels(df_alpha$d_r)[3]
# rho_ii                  = levels(df_alpha$rho_r)[3]
# chi1_ii                 = levels(df_alpha$chi1_r)[6]
# chi2_ii                 = levels(df_alpha$chi2_r)[1]
# load_comb_ii            = levels(df_alpha$load_comb_r)[3]

f_cck_ii                = levels(df_alpha$f_cck_r)[1]
d_ii                    = levels(df_alpha$d_r)[1]
rho_ii                  = levels(df_alpha$rho_r)[1]
chi1_ii                 = 0.2
chi2_ii                 = 0.6
load_comb_ii            = levels(df_alpha$load_comb_r)[1] 


df_chi1                 = df_alpha %>% filter(f_cck_r == f_cck_ii & d_r == d_ii & rho_r == rho_ii & chi2_r == chi2_ii & load_comb_r == load_comb_ii)
df_chi2                 = df_alpha %>% filter(f_cck_r == f_cck_ii & d_r == d_ii & rho_r == rho_ii & chi1_r == chi1_ii & load_comb_r == load_comb_ii)


df_chi1_new = df_chi1[,c('chi1_r', 'alpha2_r', 'alpha_labels')]
df_chi2_new = df_chi2[,c('chi2_r', 'alpha2_r', 'alpha_labels')]



g = ggplot(df_chi1_new, aes(x=chi1_r, y=alpha2_r, fill=alpha_labels))
g = g + geom_area(alpha=0.6 , size=1, colour="black")
g = g + ggtitle(paste('f_cck = ', f_cck_ii, '; d = ', d_ii, '; rho = ', rho_ii, '; chi2 = ', chi2_ii, '; load_comb = ', load_comb_ii, sep = ''))
g = g + xlab(TeX('$\\chi_1$'))
g = g + ylab(TeX('$\\alpha^2$'))
g = g + labs(fill="RV's") 
g = g + scale_fill_brewer(palette="Paired")
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g = g + scale_x_continuous(breaks = seq(0.1,0.9, by=0.05), expand = c(0,0))
g = g + scale_y_continuous(limits = seq(0,1.0), expand = c(0,0))

print(g)


# SAVE
fig_name  = paste('chi1_alpha2_', ID, '.png')
fig_path = file.path(fig_dir, fig_name)
ggsave(fig_path,
       width = 25, height = 15, units = 'cm',
       dpi = 600)



g = ggplot(df_chi2_new, aes(x=chi2_r, y=alpha2_r, fill=alpha_labels))
g = g + geom_area(alpha=0.6 , size=1, colour="black")
g = g + ggtitle(paste('f_cck = ', f_cck_ii, '; d = ', d_ii, '; rho = ', rho_ii, '; chi1 = ', chi1_ii, '; load_comb = ', load_comb_ii, sep = ''))
g = g + xlab(TeX('$\\chi_2$'))
g = g + ylab(TeX('$\\alpha^2$'))
g = g + labs(fill="RV's") 
g = g + scale_fill_brewer(palette="Paired")
g = g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g = g + scale_x_continuous(breaks = seq(0.1,0.9, by=0.1), expand = c(0,0))
g = g + scale_y_continuous(limits = seq(0,1.0), expand = c(0,0))

print(g)


# SAVE
fig_name  = paste('chi2_alpha2_', ID, '.png')
fig_path = file.path(fig_dir, fig_name)
ggsave(fig_path,
       width = 25, height = 15, units = 'cm',
       dpi = 600)
