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

################################# Results power
setwd('~/sc/Aging_multiv/')
root   <- '~/sc/Aging_multiv/'
resdir <- paste0(root, 'results/PLS_allmaps_Oct26/PLS_allmaps_power/Mainanalyses_robustness/')
res1   <- 'Powadjusted_allselected_zscore_10ncomp_pls_noleak'
res2   <- 'maps_allselected_Powadjusted_zscore_10ncomp_noleak'

# main CV results (as in your scripts)
p       <- readMat(paste0(resdir, res1, '.mat'))
vip_cv  <- as.numeric(p$average.vip)
weights <- array(unlist(p$average.weights), dim = c(55, 60))

# new results
b <- readMat(paste0(resdir, 'Powadjusted_allselected_zscore_10ncomp_bootstrap_forR.mat'))
d <- readMat(paste0(resdir, 'Powadjusted_distanceCV_spin_forR.mat'))

##### get labels/classes (ranked and sorted accordingly)
file = paste(root, 'results/PLS_allmaps_Oct26/PLS_allmaps_power/Mainanalyses_robustness/', res2, '.xlsx', sep="")
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


vip_cv  <- as.numeric(p$average.vip)
ord_vip <- order(vip_cv, decreasing = TRUE)   # = MATLAB's I: xlsx row k is map ord_vip[k]

# check: xlsx and .mat come from the same run
stopifnot(nrow(labels) == length(vip_cv),
          isTRUE(all.equal(as.numeric(labels$scores), vip_cv[ord_vip], tolerance = 1e-6)))

# undo the VIP sort -> one row per map in ORIGINAL (allmaps column) order
maps_df <- labels[order(ord_vip), c("annotation", "class")]
rownames(maps_df) <- NULL

# plotting names via annotation (order of labels_orig irrelevant)
maps_df$map <- labels_orig$tags_plotting[match(maps_df$annotation, labels_orig$annotation)]
stopifnot(nrow(maps_df) == 55, !any(is.na(maps_df$map)), !any(is.na(maps_df$class)))

# round trip: sorting maps_df by VIP must reproduce your existing labels_plot
stopifnot(all(maps_df$map[ord_vip] == labels_plot$tags_plotting))

# ONE palette for all figures (alphabetical = your Fig. 2 colours)
class_levels <- sort(unique(maps_df$class))
class_cols <- setNames(c("#a269b3","darkorange","#69b3a2","#e0ae09","#6baa58",
                         "#8a7a46","#697ab3","#699fb3","#926575","#b37d69","#4c3939"),
                       class_levels)

df_boot <- maps_df %>%
  mutate(vip     = as.numeric(b$vip.ref),
         lo      = b$vip.CI[, 1],
         hi      = b$vip.CI[, 2],
         p_vip1  = as.numeric(b$p.vip1),
         p_top10 = as.numeric(b$p.topN),
         stable  = p_vip1 >= 0.95)

cat(sprintf("Full-data vs. CV-averaged VIP: r = %.3f\n", cor(df_boot$vip, vip_cv)))

df_boot %>%
  arrange(vip) %>%
  mutate(map = factor(map, levels = map)) %>%
  ggplot(aes(x = map, y = vip)) +
  geom_linerange(aes(ymin = lo, ymax = hi), color = "grey60", linewidth = 0.9) +
  geom_point(aes(fill = class, alpha = stable), shape = 21, size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = 1) +
  scale_fill_manual(name = NULL, values = class_cols) +
  scale_alpha_manual(name = NULL, values = c(`TRUE` = 1, `FALSE` = 0.35),
                     labels = c(`TRUE` = "VIP > 1 in \u2265 95% of resamples",
                                `FALSE` = "VIP > 1 in < 95% of resamples")) +
  guides(fill  = guide_legend(override.aes = list(shape = 21, size = 5, color = "black", stroke = 0.4)),
         alpha = guide_legend(override.aes = list(shape = 21, size = 5, fill = "grey40"))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x  = element_text(angle = 60, hjust = 1, size = 15),
        axis.text.y  = element_text(size = 15),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        axis.title.x = element_text(size = 14),
        plot.margin  = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt"),
        legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 15)) +
  ylab("Variable Importance (95% bootstrap CI)") +
  xlab("Predictor Map")

ggsave('VIP_bootstrap_power_z.png', dpi = 300, width = 12, height = 9)


################################# Results AC
rm(list = ls())
setwd('~/sciebo/Aging_multiv/')
root   <- '~/sciebo/Aging_multiv/'
resdir <- paste0(root, 'results/PLS_allmaps_Oct26/PLS_allmaps_AC/Mainanalyses_robustness/')
res1   <- 'ACall_allselected_zscore_10ncomp_pls_noleak'
res2   <- 'maps_allselected_ACall_zscore_10ncomp_noleak'

# main CV results (as in your scripts)
p       <- readMat(paste0(resdir, res1, '.mat'))
vip_cv  <- as.numeric(p$average.vip)
weights <- array(unlist(p$average.weights), dim = c(55, 60)) 

# new results
b <- readMat(paste0(resdir, 'ACall_allselected_zscore_10ncomp_bootstrap_forR.mat'))
d <- readMat(paste0(resdir, 'ACall_distanceCV_spin_forR.mat'))

##### get labels/classes (ranked and sorted accordingly)
file = paste(root, 'results/PLS_allmaps_Oct26/PLS_allmaps_AC/Mainanalyses_robustness/', res2, '.xlsx', sep="")
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


vip_cv  <- as.numeric(p$average.vip)
ord_vip <- order(vip_cv, decreasing = TRUE)   # = MATLAB's I: xlsx row k is map ord_vip[k]

# check: xlsx and .mat come from the same run
stopifnot(nrow(labels) == length(vip_cv),
          isTRUE(all.equal(as.numeric(labels$scores), vip_cv[ord_vip], tolerance = 1e-6)))

# undo the VIP sort -> one row per map in ORIGINAL (allmaps column) order
maps_df <- labels[order(ord_vip), c("annotation", "class")]
rownames(maps_df) <- NULL

# plotting names via annotation (order of labels_orig irrelevant)
maps_df$map <- labels_orig$tags_plotting[match(maps_df$annotation, labels_orig$annotation)]
stopifnot(nrow(maps_df) == 55, !any(is.na(maps_df$map)), !any(is.na(maps_df$class)))

# round trip: sorting maps_df by VIP must reproduce your existing labels_plot
stopifnot(all(maps_df$map[ord_vip] == labels_plot$tags_plotting))

# ONE palette for all figures (alphabetical = your Fig. 2 colours)
class_levels <- sort(unique(maps_df$class))
class_cols <- setNames(c("#a269b3","darkorange","#69b3a2","#e0ae09","#6baa58",
                         "#8a7a46","#697ab3","#699fb3","#926575","#b37d69","#4c3939"),
                       class_levels)

df_boot <- maps_df %>%
  mutate(vip     = as.numeric(b$vip.ref),
         lo      = b$vip.CI[, 1],
         hi      = b$vip.CI[, 2],
         p_vip1  = as.numeric(b$p.vip1),
         p_top10 = as.numeric(b$p.topN),
         stable  = p_vip1 >= 0.95)

cat(sprintf("Full-data vs. CV-averaged VIP: r = %.3f\n", cor(df_boot$vip, vip_cv)))

df_boot %>%
  arrange(vip) %>%
  mutate(map = factor(map, levels = map)) %>%
  ggplot(aes(x = map, y = vip)) +
  geom_linerange(aes(ymin = lo, ymax = hi), color = "grey60", linewidth = 0.9) +
  geom_point(aes(fill = class, alpha = stable), shape = 21, size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = 1) +
  scale_fill_manual(name = NULL, values = class_cols) +
  scale_alpha_manual(name = NULL, values = c(`TRUE` = 1, `FALSE` = 0.35),
                     labels = c(`TRUE` = "VIP > 1 in \u2265 95% of resamples",
                                `FALSE` = "VIP > 1 in < 95% of resamples")) +
  guides(fill  = guide_legend(override.aes = list(shape = 21, size = 5, color = "black", stroke = 0.4)),
         alpha = guide_legend(override.aes = list(shape = 21, size = 5, fill = "grey40"))) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x  = element_text(angle = 60, hjust = 1, size = 15),
        axis.text.y  = element_text(size = 15),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        axis.title.x = element_text(size = 14),
        plot.margin  = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt"),
        legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 15)) +
  ylab("Variable Importance (95% bootstrap CI)") +
  xlab("Predictor Map")

ggsave('VIP_bootstrap_AC_z.png', dpi = 300, width = 12, height = 9)





