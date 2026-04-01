### This script plots the PLSR beta weights for each frequency bin (MEG power)
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
library(tidyverse)    # Load tidyverse for dplyr and tidyr
library(tidymodels)   # For ML mastery
library(readxl)
library(hrbrthemes)
library(R.matlab)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

############################### Plot with color-coding of the variable categories
##### set path 
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

##################################################
#### Prepare weights_long 
weights_sorted <- weights[order(vip, decreasing = TRUE), ]
weights_sorted$map <- labels_plot$tags_plotting

weights_long <- weights_sorted %>%
  pivot_longer(cols = -map, names_to = "frequency", values_to = "weight")

# convert "V1".."V60" -> numeric 1..60
weights_long$frequency <- as.numeric(gsub("V", "", weights_long$frequency))

# Calculate max abs for color scale
max_abs <- max(abs(weights_long$weight), na.rm = TRUE)

# Ensure ordering exactly as before: top map is first in map_levels
map_levels <- rev(unique(weights_sorted$map))
weights_long$map <- factor(weights_long$map, levels = map_levels)

# Build class mapping by join (robust)
class_map_df <- data.frame(map = labels_plot$tags_plotting, class = labels$class, stringsAsFactors = FALSE)
class_df <- data.frame(map = map_levels, stringsAsFactors = FALSE) %>% left_join(class_map_df, by = "map")
if (any(is.na(class_df$class))) warning("Some maps have no class mapping — check names.")

# palette and named colours
palette_vec <- c("#a269b3","darkorange","#69b3a2","#e0ae09","#6baa58",
                 "#8a7a46","#697ab3","#699fb3","#926575","#b37d69","#4c3939")
class_levels <- unique(class_df$class)
class_colors_named <- setNames(palette_vec[seq_along(class_levels)], class_levels)

# compute label anchor and right-align labels so they don't enter heatmap
minf <- min(weights_long$frequency, na.rm = TRUE)
maxf <- max(weights_long$frequency, na.rm = TRUE)
gap <- (maxf - minf) * 0.04    # make this a bit larger if labels still too close
label_x <- minf - gap          # anchor point to the left of tiles

# text dataframe with same factor levels as weights_long
text_df <- class_df %>%
  mutate(map = factor(map, levels = map_levels),
         x = label_x)

# compute max abs for colorbar (no limits that drop data)
max_abs <- max(abs(weights_long$weight), na.rm = TRUE)

# plot:
p_combined <- ggplot() +
  # heatmap: discrete y using factor 'map'
  geom_tile(data = weights_long, aes(x = frequency, y = map, fill = weight), height = 1) +
  scale_fill_viridis_c(name = "Beta weights", limits = c(-max_abs, max_abs)) +
  
  # coloured labels: right-align them (hjust = 1) and allow this layer to produce the Class legend
  geom_text(data = text_df, aes(x = x, y = map, label = map, color = class),
            hjust = 1, size = 5.1, show.legend = FALSE) +
  
  geom_point(
    data = data.frame(class = class_levels),
    aes(x = NA, y = NA, color = class),   # not plotted
    shape = 15,                          # square
    size = 5,
    inherit.aes = FALSE
  ) +
  
  # color scale + legend styling: make class legend keys appear as filled squares
  scale_color_manual(name = "Categories", values = class_colors_named,
                     guide = guide_legend(override.aes = list(shape = 15, size = 5))) +
  
  # x axis and y axis layout
  scale_x_continuous(breaks = seq(5, 55, by = 5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  
  # keep all tiles (do not set scale limits that drop rows) and allow labels outside panel
  coord_cartesian(xlim = c(label_x - 0.5, maxf + 0.5), clip = "off") +
  
  # nice theme and margins so left labels aren't clipped by canvas edges
  theme_minimal(base_size = 18) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5,
                                margin = margin(r = 140)),
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    # give extra room on the left so labels are visible (tweak if needed)
    plot.margin = margin(t = 5, r = 4, b = 5, l = 5)
  ) +
  labs(x = "Frequency (Hz)", y = "Predictor Map")

# Save (same size as your original)
outfile <- paste(root, 'results/PLS_allmaps/', res2, '_heatmap_all_fixed.png', sep = "")
ggsave(outfile, p_combined, width = 10, height = 11, dpi = 300)

p_combined



