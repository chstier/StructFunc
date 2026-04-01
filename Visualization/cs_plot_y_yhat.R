### This script plots the original, averaged powerspectrum / AC function 
### against the predicted
### Christina Stier, 2025

## R version 4.2.2 (2022-10-31)
## RStudio 2023.3.0.386 for macOS

# Load required packages

library(tidyverse)
library(qgraph)
library(corrplot)
library(igraph)
require(igraph)
library(compute.es)
library(effects)
library(multcomp)
library(pastecs)
library(psych)
library(Hmisc)
library(car)
library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(lme4)
library(optimx)
library(plyr)
library(doBy)
library(boot)
library(lmPerm)
library(R.matlab)
library(plyr)
library(abind)
library(reshape2)
library(tidyr)
library(ggseg)
library(ggplot2)
library(multcomp)
library(devtools)
library(reshape2)
library(tidyr)
library(plyr)
library(gcookbook) 
library(dplyr)
library(tidyr)
library(extrafont)

rm(list = ls())

# Read files
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/data/'

# select
res = 'absolute_10ncomp'
#res = 'zscore_10ncomp'
resname = 'Absolute (scaled) power'
resname = 'Power z-scored'

# Get y and yhat
results = readMat(paste(root, 'orig_pred_powspctrm_allselected_', res, '.mat', sep=""))
Y = as.data.frame(results$Y)
Y1 = as.data.frame(results$Y1)

# Read parcel names
regions = read_csv('/Users/path/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv')
reg = as.factor(regions$`ROI Name`)

# Now plot for 200 parcels
for (p in 1:200){

  pp = as.data.frame(t(rbind(Y[p,], Y1[p,])))
  names(pp)[1] = 'original'
  names(pp)[2] = 'predicted'
  
  pp$roi = reg[p]
  pp$frequency = 1:60
  
  # Pivot to long format
  df_long <- pp %>%
    gather(key = "spectrum", value = "value", original, predicted)
  
  # Plot
  g = ggplot(df_long, aes(x = frequency, y = value, color = spectrum)) +
    geom_line(aes(colour=spectrum)) + 
    #geom_ribbon(aes(ymin=low_freq, ymax=up_freq, fill = spectrum), alpha=0.15, colour = NA) + 
    geom_path(size = 1.2) +
    # scale_x_log10() + 
    # scale_y_log10() + 
    xlab("Frequency (Hz)") + 
    ylab(resname) +
    scale_color_manual(values=c("#293352", "#D16103"), 
                       labels=c("Original", "Predicted"),
                       breaks=c("original", "predicted"), name = "Spectrum") + 
    scale_fill_manual(values=c("#293352", "#D16103"), breaks=c("original", "predicted"), guide="none") +
    theme_pubr(base_size = 20, legend = c(0.8, 0.8)) 
  
  ggsave(file=paste("Spctrm_", res, reg[p],'.png', sep=''),plot=g, dpi = 300, limitsize = TRUE, width = 6, height = 6)
  # "#293352", "#C3D7A4"
}



################################### Same for autocorrelation
# Read files
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/data/'

# select
#res = 'abs_10ncomp'
res = '10ncomp_noleak'
#resname = 'Absolute (scaled) power'
resname = 'Autocorrelation'

# Get y and yhat
results = readMat(paste(root, 'orig_pred_ACall_allselected_', res, '.mat', sep=""))
Y = as.data.frame(results$Y)
Y1 = as.data.frame(results$Y1)

# Read parcel names
regions = read_csv('/Users/path/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv')
reg = as.factor(regions$`ROI Name`)

# Now plot for 200 parcels
for (p in 1:200){
  
  pp = as.data.frame(t(rbind(Y[p,], Y1[p,])))
  names(pp)[1] = 'original'
  names(pp)[2] = 'predicted'
  
  pp$roi = reg[p]
  pp$frequency = 1:40
  
  # Pivot to long format
  df_long <- pp %>%
    gather(key = "autocorrelation", value = "value", original, predicted)
  
  # Plot
  g = ggplot(df_long, aes(x = frequency, y = value, color = autocorrelation)) +
    geom_line(aes(colour=autocorrelation)) + 
    #geom_ribbon(aes(ymin=low_freq, ymax=up_freq, fill = spectrum), alpha=0.15, colour = NA) + 
    geom_path(size = 1.2) +
    # scale_x_log10() + 
    # scale_y_log10() + 
    xlab("Time lags") + 
    ylab(resname) +
    scale_color_manual(values=c("#293352", "#D16103"), 
                       labels=c("original", "predicted"),
                       breaks=c("original", "predicted"), name = "autocorrelation") + 
    scale_fill_manual(values=c("#293352", "#D16103"), breaks=c("original", "predicted"), guide="none") +
    theme_pubr(base_size = 20, legend = c(0.8, 0.8)) 
  
  ggsave(file=paste("AC", res, reg[p],'.png', sep=''),plot=g, dpi = 300, limitsize = TRUE, width = 6, height = 6)
  # "#293352", "#C3D7A4"
}




