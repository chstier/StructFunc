### This script plots ......
### Christina Stier, 2025

## R version 4.2.2 (2022-10-31)
## RStudio 2023.3.0.386 for macOS

rm(list = ls())

## load packages
install.packages("corrplot")
install.packages("igraph")
install.packages("qgraph")
install.packages("car")
install.packages("compute.es")
install.packages("effects")
install.packages("compute.es")
install.packages("ggplot2")
install.packages("multcomp")
install.packages("pastecs")
install.packages("psych")
install.packages("Hmisc")
install.packages("Rcmdr")
install.packages("splines")
install.packages("gridExtra")
install.packages("grid")	
install.packages("ggpubr")
install.packages("cowplot")
install.packages("optimx")
install.packages("plyr")
install.packages("doBy")
install.packages("boot")
install.packages("lmPerm")
install.packages('R.matlab')
install.packages('abind')
install.packages("devtools")
install.packages("plotrix")
install.packages("gcookbook")
install.packages("hrbrthemes")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggseg")
install.packages("readxl")

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
library(tidyverse)    
library(tidymodels)   
library(readxl)
library(hrbrthemes)


##### set path for power results
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/'
root1 = '~/sciebo/Aging_multiv/'
res1 = 'ageeffects_powerspctrmCorr_weights'
#res2 = 'maps_allselected_Powadjusted_zscore_10ncomp_noleak'

##### get age effects on PLS weights
p = readMat(paste(root,'results/PLS_allmaps/Individual/', res1, '.mat', sep=""))

r = as.data.frame(array(unlist(p$r), dim = c(55, 60)))
w = p$weights
r2 = p$r2.all

##### get labels/classes in original order
file = paste(root1, 'maps_classes_selected.csv', sep="")
labels_orig1 = read.table(file, header = TRUE, sep=";")
labels_orig1 = as.data.frame(labels_orig1)

file = paste(root1, 'maps_classes_selected_plotting.csv', sep="")
labels_orig = read.table(file, header = TRUE, sep=";")
labels_orig = as.data.frame(labels_orig)

# # reorder plotting label according to the results
labels_plot <- labels_orig[match(labels_orig1$annotation, labels_orig$annotation), ]

################################################## Plot BETA WEIGHTS
# check in which maps show the largest age effects
listr = as.numeric(apply(abs(r), 1, mean))

weights_sorted <- r[order(listr, decreasing = TRUE), ]
labels_sorted = labels_plot[order(listr, decreasing = TRUE), ]

# Add map labels from the vip-sorted labels
weights_sorted$map <- labels_sorted$tags_plotting

weights_long <- weights_sorted %>%
  pivot_longer(cols = -map, names_to = "frequency", values_to = "weight")

# Make sure 'map' is a factor with desired order
weights_long$map <- factor(weights_long$map, levels = rev(unique(weights_sorted$map)))
weights_long$frequency <- factor(weights_long$frequency, levels = unique(weights_long$frequency))
weights_long$frequency <- as.numeric(gsub("V", "", weights_long$frequency))


# Calculate the maximum absolute age effect
max_abs <- max(abs(weights_long$weight))

# Adapt scaling to AC results to provide head-to-head comparison
ggplot(weights_long, aes(x = frequency, y = map, fill = weight)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B", limits = c(-0.66, 0.66)) +
  scale_x_continuous(breaks = seq(5, 55, by = 5)) +  # only label every 5th frequency
  labs(x = "Frequency (Hz)", y = "Map", fill = "Age effects (r) on weights") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave('~/sciebo/path/ageeffects_weights_powCorr_scaled2.png', width = 10, height = 11, dpi = 300)

###### Load age
load("~/sciebo/Age/R_camcan/datasets/coh_img_data_absmean_global_allcohorts.Rda")

head(all_imcoh)
str(all_imcoh)

# change format of a few variables
age = as.numeric(all_imcoh$age)

# get map of interest
idx = which(labels_orig$tags == 'COX-1')
l = 'COX-1'
data = as.data.frame(p$weights[,idx,15])
data = cbind(data, age)
names(data)[1] = 'map'

freqname = '15Hz'
name = paste("Individual beta weights", l, freqname, sep=" ") 

g = ggplot(data,aes(x = age, y = map)) + geom_point(aes(col="#393939"), alpha = 1) + #stat_smooth(method = "lm", formula = y ~ x + I(x^2), se= FALSE, size = 1.5, color = "gray") + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1, color = "#525252") + labs(y=name, x="Age")
g = g + scale_color_viridis_d(option = "inferno", guide = FALSE) 
g = g + scale_x_continuous(breaks=c(20,30,40,50,60,70,80))
g <- g + theme(text = element_text(size = 22), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))	
g


freq = 1:60
idx = which(labels_orig$tags == 'ev-expansion ')
l = 'ev-expansion'


for (i in 1:length(freq)){
  freqname = freq[i]
  #data = data_long[data_long$freq == freqname,]
  data = as.data.frame(p$weights[,idx,freq[i]])
  data = cbind(data, age)
  names(data)[1] = 'map'
  
  #name = "Indidivual weights (PLSR, Power)"
  name = paste("Individual beta weights", l, freqname, sep=" ") 
  g = ggplot(data,aes(x = age, y = map)) + geom_point(aes(col="#393939"), alpha = 1) + #stat_smooth(method = "lm", formula = y ~ x + I(x^2), se= FALSE, size = 1.5, color = "gray") + 
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1, color = "#525252") + labs(y=name, x="Age")
  g = g + scale_color_viridis_d(option = "inferno", guide = FALSE) 
  g = g + scale_x_continuous(breaks=c(20,30,40,50,60,70,80))
  g <- g + theme(text = element_text(size = 22), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))	
    
  ggsave(paste("ageeffects_map_", l,"_" ,freq[i], ".png", sep=""), plot=g, dpi = 300, limitsize = TRUE, width = 8, height = 8)
  #plot_list[[i]] = g
}


dat = as.data.frame(cbind(age,r2))
names(dat)[2] = 'r2'
ggplot(dat, aes(x = age, y = r2)) + 
  geom_point(color = "#393939", alpha = 0.6, size = 2.5) + 
  stat_smooth(method = "lm", formula = y ~ x, 
              se = FALSE, size = 1.2, color = "#525252", linetype = "dashed") +
  labs(y = "Individual deviations: predicted R² (Power)", x = "Age") +
  scale_x_continuous(breaks = seq(20, 80, by = 10)) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 24),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5)
  )

ggsave('~/sciebo/path/overall_individualr2_power.png', width = 10, height = 8, dpi = 300)

mean_list = list()
for (i in 1:length(weights_sorted[,1])){
  absm = mean(as.numeric(abs(weights_sorted[i,1:60])))
  mean_list[[i]] = absm
}


##########################################################################
###### Same for AC ####################################################### 
rm(list = ls())
##### set path
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/'
root1 = '~/sciebo/Aging_multiv/'
res1 = 'ageeffects_ACall_weights'
#res2 = 'maps_allselected_Powadjusted_zscore_10ncomp_noleak'

##### get age effects on PLS weights
p = readMat(paste(root, 'results/PLS_AC/ACall_individual/', res1, '.mat', sep=""))
r = as.data.frame(array(unlist(p$r), dim = c(55, 40)))
w = p$weights
r2 = p$r2.all

##### get labels/classes in original order
file = paste(root1, 'maps_classes_selected.csv', sep="")
labels_orig1 = read.table(file, header = TRUE, sep=";")
labels_orig1 = as.data.frame(labels_orig1)

file = paste(root1, 'maps_classes_selected_plotting.csv', sep="")
labels_orig = read.table(file, header = TRUE, sep=";")
labels_orig = as.data.frame(labels_orig)

# # reorder plotting label according to the results
labels_plot <- labels_orig[match(labels_orig1$annotation, labels_orig$annotation), ]

##### read feature labels
v = (1:40)
ms = (v/300)*1000  
ms = round(ms, digits = 1)
time = as.factor(ms)

################################################## Plot BETA WEIGHTS
# check in which maps show the largest age effects
listr = as.numeric(apply(abs(r), 1, mean))

weights_sorted <- r[order(listr, decreasing = TRUE), ]
labels_sorted = labels_plot[order(listr, decreasing = TRUE), ]

# # Add map labels from the vip-sorted labels
# weights_sorted$map <- labels_sorted$tags
# Add map labels from the vip-sorted labels
weights_sorted$map <- labels_sorted$tags_plotting

#weights_sorted$time <- time

weights_long <- weights_sorted %>%
  pivot_longer(cols = -map, names_to = "frequency", values_to = "weight")

# Make sure 'map' is a factor with desired order
weights_long$map <- factor(weights_long$map, levels = rev(unique(weights_sorted$map)))
weights_long$frequency <- factor(weights_long$frequency, levels = unique(weights_long$frequency))
weights_long$frequency <- as.numeric(gsub("V", "", weights_long$frequency))
weights_long$time <- rep(time,55)

# Calculate the maximum absolute age effect
max_abs <- max(abs(weights_long$weight))

ggplot(weights_long, aes(x = time, y = map, fill = weight)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B", limits = c(-0.66, 0.66)) +
  scale_x_discrete(breaks = c(3.3, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130)) +
  labs(x = "Autocorrelation time lags (ms)", y = "Map", fill = "Age effects (r) on weights") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave('~/sciebo/path/ageeffects_weights_ACall_scaled2.png', width = 10, height = 11, dpi = 300)

mean_list = list()
for (i in 1:length(weights_sorted[,1])){
  absm = mean(as.numeric(abs(weights_sorted[i,1:40])))
  mean_list[[i]] = absm
}


###### Load age
load("~/sciebo/Age/path/coh_img_data_absmean_global_allcohorts.Rda")

head(all_imcoh)
str(all_imcoh)

# change format of a few variables
age = as.numeric(all_imcoh$age)

# get map of interest
idx = which(labels_orig$tags == 'COX-1')
l = 'COX-1'
data = as.data.frame(p$weights[,idx,15])
data = cbind(data, age)
names(data)[1] = 'map'

freqname = '15Hz'
name = paste("Individual beta weights", l, freqname, sep=" ") 

g = ggplot(data,aes(x = age, y = map)) + geom_point(aes(col="#393939"), alpha = 1) + #stat_smooth(method = "lm", formula = y ~ x + I(x^2), se= FALSE, size = 1.5, color = "gray") + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1, color = "#525252") + labs(y=name, x="Age")
g = g + scale_color_viridis_d(option = "inferno", guide = FALSE) 
g = g + scale_x_continuous(breaks=c(20,30,40,50,60,70,80))
g <- g + theme(text = element_text(size = 22), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))	
g


freq = 1:40
idx = which(labels_orig$tags == 'ev-expansion ')
l = 'ev-expansion'


for (i in 1:length(freq)){
  freqname = freq[i]
  #data = data_long[data_long$freq == freqname,]
  data = as.data.frame(p$weights[,idx,freq[i]])
  data = cbind(data, age)
  names(data)[1] = 'map'
  
  
  name = paste("Individual beta weights", l, freqname, sep=" ") 
  #name = "Indidivual weights (PLSR, AC)"
  g = ggplot(data,aes(x = age, y = map)) + geom_point(aes(col="#393939"), alpha = 1) + #stat_smooth(method = "lm", formula = y ~ x + I(x^2), se= FALSE, size = 1.5, color = "gray") + 
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1, color = "#525252") + labs(y=name, x="Age")
  g = g + scale_color_viridis_d(option = "inferno", guide = FALSE) 
  g = g + scale_x_continuous(breaks=c(20,30,40,50,60,70,80))
  g <- g + theme(text = element_text(size = 22), legend.title = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50"), plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))	
  
  ggsave(paste("ageeffects_map_", l,"_" ,freq[i], ".png", sep=""), plot=g, dpi = 300, limitsize = TRUE, width = 8, height = 8)
  #plot_list[[i]] = g
}

dat = as.data.frame(cbind(age,r2))
names(dat)[2] = 'r2'
ggplot(dat, aes(x = age, y = r2)) + 
  geom_point(color = "#393939", alpha = 0.6, size = 2.5) + 
  stat_smooth(method = "lm", formula = y ~ x, 
              se = FALSE, size = 1.2, color = "#525252", linetype = "dashed") +
  labs(y = "Individual-level predicted R² (AC)", x = "Age") +
  scale_x_continuous(breaks = seq(20, 80, by = 10)) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 24),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5)
  )

ggsave('~/sciebo/path/overall_individualr2_AC.png', width = 10, height = 8, dpi = 300)


