### This script plots the PLSR beta weights for time lage (MEG autocorrelation)
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
library(viridis)

##### set path
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/'
res1 = 'ACall_allselected_10ncomp_pls_noleak'
res2 = 'maps_allselected_ACall_10ncomp_noleak'

############################### Plot with color-coding of the variable categories
##### get PLS results
p = readMat(paste(root, 'results/PLS_AC/ACall/', res1, '.mat', sep=""))

r2 = as.numeric(p$average.R2)
rank = as.numeric(p$average.rank)
vip = as.numeric(p$average.vip)
weights = as.data.frame(array(unlist(p$average.weights), dim = c(55, 40)))
rmse_parcel = as.data.frame(array(unlist(p$rmse.parcel), dim = c(200, 50)))
rmse_parcel = rowMeans(rmse_parcel)

##### get labels/classes (ranked and sorted accordingly)
file = paste(root, 'results/PLS_AC/ACall/', res2, '.xlsx', sep="")
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

## Build weights_long 
weights_sorted <- weights[order(vip, decreasing = TRUE), ]

# attach plotting-friendly map names from labels_plot
weights_sorted$map <- labels_plot$tags_plotting

# Reshape to long format
weights_long <- weights_sorted %>%
  pivot_longer(cols = -map, names_to = "AClag", values_to = "weight")

weights_long$AClag <- as.integer(gsub("^V", "", weights_long$AClag))

# convert from lags to time-axis (ms)
v = (1:40)
ms = (v/300)*1000  
ms = round(ms, digits = 1)
time = as.factor(ms)

# Make sure 'map' is a factor with desired order
#weights_long$map <- factor(weights_long$map, levels = rev(unique(weights_sorted$map)))
weights_long$AClag <- factor(weights_long$AClag, levels = unique(weights_long$AClag))
weights_long$AClag <- as.numeric(gsub("V", "", weights_long$AClag))
weights_long$time <- rep(time,55)

# Preserve ordering exactly like your original code (top = first element)
map_levels <- rev(unique(weights_sorted$map))
weights_long$map <- factor(weights_long$map, levels = map_levels)

# compute max abs for colour scale (safe to set as scale limits)
max_abs <- max(abs(weights_long$weight), na.rm = TRUE)

## FIXED named palette 
class_palette <- c(
  "cortical expansion"           = "#69b3a2",
  "gene expression - summary"    = "#4c3939",
  "gradient - summary"           = "#b37d69",
  "metabolism"                   = "#a269b3",
  "microstructural"              = "#6baa58",
  "network-metrics"              = "#926575",
  "neurotransmitter - enzyme"    = "#699fb3",
  "neurotransmitter - protein"   = "#8a7a46",
  "neurotransmitter - receptor"  = "#e0ae09",
  "neurotransmitter - summary"   = "#697ab3",
  "neurotransmitter - transporter" = "darkorange"
)

# build class mapping
class_map_df <- data.frame(map = labels_plot$tags_plotting, class = labels_plot$class, stringsAsFactors = FALSE)

class_df <- data.frame(map = map_levels, stringsAsFactors = FALSE) %>%
  left_join(class_map_df, by = "map")

if (any(is.na(class_df$class))) {
  warning("Some maps do not have a class in labels_plot. Check names and canonicalization.")
}

# Force the class factor to follow the palette order (so legend order is stable)
class_df$class <- factor(class_df$class, levels = names(class_palette))

# text_df for geom_text (label anchor X to the left)
minf <- min(as.numeric(weights_long$time), na.rm = TRUE)
maxf <- max(as.numeric(weights_long$time), na.rm = TRUE)
gap <- (maxf - minf) * 0.06   # adjust if you want more/less left spacing
label_x <- minf - gap

text_df <- class_df %>%
  mutate(map = factor(map, levels = map_levels),
         x = label_x)

# ensure text_df$class is factor with same levels as palette
text_df$class <- factor(as.character(text_df$class), levels = names(class_palette))


## Plot
p_combined <- ggplot() +
  
  # main heatmap (discrete y = factor map)
  geom_tile(data = weights_long,
            aes(x = time, y = map, fill = weight),
            height = 1) +
  
  # continuous colour for weights (we allow symmetric range from data)
  scale_fill_viridis_c(name = "Beta weights", limits = c(-max_abs, max_abs)) +
  
  # coloured labels (right-aligned so they sit left of the tiles)
  geom_text(data = text_df,
            aes(x = x, y = map, label = map, color = class),
            hjust = 1, size = 5.1, show.legend = FALSE) +
  
  # dummy points to build legend keys (will not be plotted)
  geom_point(data = data.frame(class = names(class_palette)),
             aes(x = NA, y = NA, color = class),
             shape = 15, size = 5, inherit.aes = FALSE) +
  
  # use the fixed named palette (enforces color <-> category mapping and order)
  scale_color_manual(
    name = "Categories",
    values = class_palette,
    breaks = names(class_palette),  # explicit order in legend
    drop = FALSE,
    guide = guide_legend(override.aes = list(shape = 15, size = 5))
  ) +
  
  # axes and limits; we use coord_cartesian to avoid dropping data
  #scale_x_continuous(breaks = seq(5, 55, by = 5), expand = c(0,0)) +
  scale_x_discrete(breaks = c(3.3, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130)) +
  scale_y_discrete(expand = c(0,0)) +
  coord_cartesian(xlim = c(label_x - 0.5, maxf + 0.5), clip = "off") +
  
  # theme + margins: tune left margin so labels are visible
  labs(x = "Autocorrelation time lags (ms)", y = "Predictor Map") +
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
    plot.margin = margin(t = 5, r = 4, b = 5, l = 5)  # increase left space if labels get clipped
  )

# Save output (adjust file name if desired)
outfile <- paste0(root, 'results/PLS_AC/ACall/', res2, '_heatmap_all_fixed_palette.png')
ggsave(outfile, p_combined, width = 10, height = 11, dpi = 300)

# also show in the active R session
p_combined
