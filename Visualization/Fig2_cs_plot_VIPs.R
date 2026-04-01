### This script plots the variable importance of structural and molecular 
### predictor maps after running PLSR on MEG power data 
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

##### set paths
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/'
res1 = 'Powadjusted_allselected_zscore_10ncomp_pls_noleak'
res2 = 'maps_allselected_Powadjusted_zscore_10ncomp_noleak'

##### get PLS results
p = readMat(paste(root, 'results/PLS_allmaps/', res1, '.mat', sep=""))

r2 = as.numeric(p$average.R2)
rank = as.numeric(p$average.rank)
vip = as.numeric(p$average.vip)
weights = as.data.frame(array(unlist(p$average.weights), dim = c(55, 60)))
rmse_parcel = as.data.frame(array(unlist(p$rmse.parcel), dim = c(200, 50)))
rmse_parcel = rowMeans(rmse_parcel)

##### get labels/classes (ranked and sorted accordingly)
file = paste(root, 'results/PLS_allmaps/', res2, '.xlsx', sep="")
labels = read_xlsx(file)
labels = as.data.frame(labels)
idx = (labels$class == 'microstructural - summary')
labels$class[idx == TRUE] = 'gene expression - summary' # relabel gene expression

##### get labels/classes in original order
file = paste(root, 'maps_classes_selected_plotting.csv', sep="")
labels_orig = read.table(file, header = TRUE, sep=";")
labels_orig = as.data.frame(labels_orig)

# reorder plotting label according to the results
labels_plot <- labels_orig[match(labels$annotation, labels_orig$annotation), ]


################### VIP scores per map
df_vip <- data.frame(
  map = labels_plot$tags_plotting,
  class = labels$class,  # Adjust if the column name is different
  vip_score = labels$scores
)


df_vip %>%
  arrange(vip_score) %>%
  mutate(map=factor(map, levels=map)) %>% 
  ggplot(aes(x=map, y=vip_score)) +
  geom_segment( aes(x=map,xend=map, y=0.09, yend=abs(vip_score)), color="grey") +
#  scale_y_continuous(breaks=seq(0.1,0.7,0.1)) +
#  geom_line( color="grey") +
  geom_point(shape=21, aes(fill = factor(class)), size=4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", size = 1) +  # Horizontal line at VIP = 1
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 5, color = "black", stroke = 0.4)
  )) +
  scale_fill_manual(name = NULL, values=c("#a269b3", "darkorange", "#69b3a2", "#e0ae09", "#6baa58", "#8a7a46", "#697ab3", "#699fb3", "#926575", "#b37d69", "#4c3939")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 15), 
    axis.text.y = element_text(size = 15), 
    axis.title.y = element_text(size = 15, 
    margin = margin(t = 0, r = 8, b = 0, l = 0)), 
    axis.title.x = element_text(size = 14), 
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt"),
    #plot.margin = margin(t=1, r=1, b=1, l=2, unit = "cm"),
    legend.position = "bottom",
    #legend.title = element_text(size = 12),
    legend.text = element_text(size = 15)) +
  ylab('Variable Importance') +
  xlab('Predictor Map') 

ggsave('VIP_colored_power_z.png', dpi = 300, width = 12, height = 8)

