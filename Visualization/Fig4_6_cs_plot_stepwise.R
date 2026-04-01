### This script plots the results for the forward selection of brain maps 
### (based on incremental PLSR procedure)
### Christina Stier, 2025

## R version 4.2.2 (2022-10-31)
## RStudio 2023.3.0.386 for macOS

# Load required packages
library(ggplot2)
library(readr)
library(ggrepel)
library(dplyr)
library(stringr)
library(scales)

rm(list = ls())

# Read the CSV output
# Power results
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/results/PLS_allmaps/'
basic = '~/sciebo/Aging_multiv/'
res1 = 'pow_zscored_noleak'
results <- read_csv(paste(root, 'forward_selection_', res1, '.csv', sep=""))

# or AC results
setwd('~/sciebo/Aging_multiv/')
root = '~/sciebo/Aging_multiv/results/PLS_AC/ACall/'
basic = '~/sciebo/Aging_multiv/'
res1 = 'ACall_noleak'
results <- read_csv(paste(root, 'forward_selection_results_', res1, '.csv', sep=""))


# Get class labels and sort according to results map
file = paste(basic, 'maps_classes_selected_plotting.csv', sep="")
labels_orig = read.table(file, header = TRUE, sep=";")
labels_orig = as.data.frame(labels_orig)
labels_plot <- labels_orig[match(results$MapLabel, labels_orig$annotation), ]
results$Tag <- labels_plot$tags_plotting

# make sure results exists and is a plain data.frame
results <- as.data.frame(results)

# safe index & short tag
results$idx <- seq_len(nrow(results))
results$Tag_short <- str_trunc(as.character(results$Tag), 18, ellipsis = "...")

# scale factor for secondary axis (safe with na.rm)
scale_factor <- max(results$CV_R2, na.rm = TRUE) / max(results$NumComponents, na.rm = TRUE)

# ensure the left axis includes 0.8
y_max <- 1.05 * max(results$CV_R2, na.rm = TRUE)
breaks_left <- pretty(c(0, y_max), n = 6)
breaks_left <- sort(unique(c(breaks_left, 0.8)))

p1 <- ggplot(results, aes(x = idx)) +
  # main R^2 line + points
  geom_line(aes(y = CV_R2), color = "#2c7fb8", size = 1) +
  geom_point(aes(y = CV_R2), color = "#2c7fb8", size = 2) +
  # number of components (scaled) as dashed line
  geom_line(aes(y = NumComponents * scale_factor), color = "#006400", linetype = "dashed", size = 0.9) +
  geom_point(aes(y = NumComponents * scale_factor), color = "#006400", size = 1.8) +
  # non-overlapping labels (short names)
  geom_text_repel(data = results,
                  aes(y = CV_R2, label = Tag_short),
                  nudge_x = 0.6,
                  size = 3.5,
                  segment.size = 0.25,
                  segment.color = "grey60",
                  max.overlaps = 40,
                  box.padding = 0.35,
                  direction = "y",
                  hjust = 0) +
  # x axis: steps of 5 (5, 10, 15, ...)
  scale_x_continuous(name = "Number of Maps Added",
                     breaks = seq(5, nrow(results), by = 5),
                     limits = c(1, nrow(results)),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(
    name = "Cross-validated R²",
    limits = c(0, y_max),
    breaks = breaks_left,
    sec.axis = sec_axis(~ . / scale_factor, name = "# PLS Components", breaks = pretty_breaks(6))
  ) +
  labs(title = "Forward Selection of Brain Maps (PLS)") +
  theme_light(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y.right = element_text(color = "#006400"),
    axis.title.y.left = element_text(color = "#2c7fb8")
  )

# save + show
ggsave("forward_selection_improved_labels.png", p1, width = 11, height = 8, dpi = 300)
p1

