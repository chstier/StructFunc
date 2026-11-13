### Visualize MEG age effects and their structural/neurochemical correlates
### PLS coefficients (all maps) + VIP with bootstrap CI
### Power and AC
### Christina Stier, 2026
## R version 4.2.2

rm(list = ls())

library(R.matlab)
library(tidyverse)
library(patchwork)
library(scales)
library(readxl)
if ('package:plyr' %in% search()) detach('package:plyr', unload = TRUE)

## ---- settings ----
root       <- '~/sciebo/Aging_multiv/'
resdir     <- paste0(root, 'results/agemap_standalone/')
outdir     <- paste0(resdir, 'figures/')
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
alphaFDR   <- 0.05
nParc      <- 200
metrics    <- c('Pow_zY', 'AC_zY')
showAgeBar <- TRUE     # panel a: % parcels with age-related increase / decrease
showBands  <- FALSE     # Greek band labels above panel a (power only; no band borders)
markVIPsig <- FALSE     # asterisk if VIP spin-FDR q < .05
orderBy    <- 'own'    # 'own' = VIP of this metric; 'mean' = mean VIP (same order in both figures)
shortClass <- FALSE     # 'neurotransmitter - receptor' -> 'NT receptor' etc.

## ---- print size (A4) ----
page_w <- 7.2          # in (183 mm = full A4 text width)
max_h  <- 9.4          # in; leaves room for the legend text on the page
top_h  <- 1.15         # in, panel a
leg_h  <- 1.0          # in, legends at bottom
fs     <- 7            # base font size (pt) at print size
fs_map <- 9.8         # map labels (pt)

cfgs <- list(
  Pow = list(xlab = 'Frequency (Hz)', xbreak = 10, bands = TRUE),
  AC  = list(xlab = 'Lag (ms)',       xbreak = 20, bands = FALSE))

## ---- style ----
div_cols <- c('#053061','#2166AC','#4393C3','#92C5DE','#D1E5F0','#F7F7F7',
              '#FDDBC7','#F4A582','#D6604D','#B2182B','#67001F')
scale_fill_div <- function(name, lims,
                           guide = guide_colourbar(title.position = 'top', title.hjust = 0.5), ...)
  scale_fill_gradientn(colours = div_cols, limits = lims, oob = squish, name = name,
                       guide = guide, ...)
col_pos <- '#B2182B'; col_neg <- '#2166AC'
theme_fig <- function(base = fs) theme_classic(base_size = base) +
  theme(axis.text    = element_text(colour = 'black', size = base),
        axis.title   = element_text(size = base + 0.5),
        axis.line    = element_line(linewidth = 0.3),
        axis.ticks   = element_line(linewidth = 0.3),
        plot.tag     = element_text(face = 'bold', size = base + 3),
        legend.title = element_text(size = base),
        legend.text  = element_text(size = base - 0.5),
        plot.margin  = margin(2, 2, 2, 2))
save_fig <- function(name, plot, w, h) {
  ggsave(paste0(outdir, name, '.png'), plot, width = w, height = h, dpi = 300, bg = 'white')
  try(ggsave(paste0(outdir, name, '.pdf'), plot, width = w, height = h, device = cairo_pdf),
      silent = TRUE)
}
palette_vec <- c("#a269b3","darkorange","#69b3a2","#e0ae09","#6baa58",
                 "#8a7a46","#697ab3","#699fb3","#926575","#b37d69","#4c3939")
getM <- function(x, nr) matrix(as.numeric(x), nrow = nr)

## ---- map labels and classes (MATLAB order) ----
clean <- function(x) trimws(gsub("^'|'$", "", trimws(x)))
ml <- read.csv(paste0(resdir, 'agemap_maplabels.csv'), stringsAsFactors = FALSE)
ml$annotation <- clean(ml$annotation)
nMaps <- nrow(ml)

labels <- as.data.frame(read_xlsx(paste0(root, 'maps_classes_orig_final.xlsx')))
labels[] <- lapply(labels, function(col) if (is.character(col)) clean(col) else col)
labels$class[labels$class == 'microstructural - summary'] <- 'gene expression - summary'
labels_orig <- read.table(paste0(root, 'maps_classes_selected_plotting.csv'),
                          header = TRUE, sep = ";", quote = "\"", fileEncoding = "UTF-8-BOM",
                          stringsAsFactors = FALSE, strip.white = TRUE)
labels_orig <- labels_orig[labels_orig$annotation != "", ]
labels_plot <- labels_orig[match(labels$annotation, labels_orig$annotation), ]
miss <- labels$annotation[is.na(labels_plot$annotation)]
if (length(miss) > 0) stop('No plotting label for: ', paste(miss, collapse = ', '))

maps <- data.frame(annotation    = labels$annotation,
                   tags_plotting = trimws(labels_plot$tags_plotting),
                   class         = sub('proteine', 'protein', labels$class),
                   stringsAsFactors = FALSE)
stopifnot(nrow(maps) == nMaps, setequal(ml$annotation, maps$annotation),
          !anyDuplicated(maps$tags_plotting))
maps <- maps[match(ml$annotation, maps$annotation), ]
rownames(maps) <- NULL
class_levels <- sort(unique(maps$class))
class_cols   <- setNames(palette_vec[seq_along(class_levels)], class_levels)  # same as Fig. 7
class_labs   <- setNames(if (shortClass)
  sub(' - summary$', ' (summary)', sub('^neurotransmitter - ', 'NT ', class_levels)) else
    class_levels, class_levels)

## ---- load results ----
load_res <- function(prefix) {
  p <- readMat(paste0(resdir, prefix, '_agemap_for_R.mat'))
  need <- c('L.vip', 'L.qVIP', 'Lb.ci', 'Lw.beta', 'lin.r', 'lin.q')
  miss <- need[!need %in% names(p)]
  if (length(miss)) stop(prefix, ': missing fields ', paste(miss, collapse = ', '))
  r <- list(xv = as.numeric(p$freqs), beta = getM(p$Lw.beta, nMaps),
            vip = as.numeric(p$L.vip), qVIP = as.numeric(p$L.qVIP), ci = getM(p$Lb.ci, nMaps),
            linr = getM(p$lin.r, nParc), linq = getM(p$lin.q, nParc))
  stopifnot(ncol(r$beta) == length(r$xv), length(r$vip) == nMaps)
  r
}
RES     <- lapply(setNames(metrics, metrics), load_res)
meanVIP <- rowMeans(sapply(RES, `[[`, 'vip'))

## ---- row height so that the figure fits on the page ----
row_h <- min(0.16, (max_h - (if (showAgeBar) top_h else 0) - leg_h - 0.35) / nMaps)
cat(sprintf('Row height: %.1f pt (map labels %.1f pt)\n', row_h * 72, fs_map))

## =====================================================================
## Figure function
## =====================================================================
make_fig <- function(prefix) {
  cf <- cfgs[[sub('_.*$', '', prefix)]]; r <- RES[[prefix]]     # 'Pow_zY' -> 'Pow'
  if (is.null(cf)) stop('No configuration for ', prefix)
  xv <- r$xv; nF <- length(xv)
  
  xstep  <- median(diff(xv)); xspan <- diff(range(xv))
  x_lo   <- min(xv) - xstep / 2; x_hi <- max(xv) + xstep / 2
  x_cls  <- x_lo - 0.025 * xspan
  x_left <- x_cls - 0.02 * xspan
  xbreaks <- seq(cf$xbreak, max(xv), by = cf$xbreak)
  geom_cells <- if (sd(diff(xv)) < 1e-6) geom_raster else geom_tile
  
  ord <- order(if (orderBy == 'mean') meanVIP else r$vip, decreasing = TRUE)
  lv  <- rev(maps$tags_plotting[ord])
  
  ## ---- (b) PLS coefficients, unthresholded ----
  dfW <- expand.grid(m = seq_len(nMaps), fi = seq_len(nF)) %>%
    mutate(b = as.vector(r$beta), x = xv[fi],
           map = factor(maps$tags_plotting[m], levels = lv))
  cls_df <- maps %>% transmute(map = factor(tags_plotting, levels = lv), class)
  lim <- max(abs(dfW$b))
  brk <- signif(lim * 0.9, 1)
  
  p_beta <- ggplot(dfW, aes(x, map)) +
    geom_cells(aes(fill = b)) +
    geom_point(data = cls_df, aes(x = x_cls, y = map, colour = class),
               shape = 15, size = 2, inherit.aes = FALSE) +
    scale_fill_div('PLS coefficient', c(-lim, lim),
                   breaks = c(-brk, 0, brk),
                   labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE),
                   guide = guide_colourbar(title.position = 'top', title.hjust = 0.5, order = 1,
                                           barwidth = unit(2.8, 'cm'), barheight = unit(0.25, 'cm'))) +
    scale_colour_manual(name = 'Map category', values = class_cols,
                        breaks = class_levels, labels = class_labs[class_levels],
                        guide = guide_legend(ncol = 3, title.position = 'top', order = 2,
                                             byrow = FALSE,
                                             override.aes = list(size = 2.5),
                                             keywidth = unit(0.28, 'cm'),
                                             keyheight = unit(0.28, 'cm'))) +
    scale_x_continuous(breaks = xbreaks, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(xlim = c(x_left, x_hi), expand = FALSE, clip = 'off') +
    labs(x = cf$xlab, y = NULL) + theme_fig() +
    theme(axis.line = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = fs_map))
  
  ## ---- (c) VIP with 95% bootstrap CI ----
  df_vip <- maps %>%
    mutate(vip = r$vip, lo = r$ci[, 1], hi = r$ci[, 2], sig = r$qVIP < alphaFDR,
           map = factor(tags_plotting, levels = lv))
  p_vip <- ggplot(df_vip, aes(vip, map)) +
    geom_vline(xintercept = 1, linetype = 'dashed', colour = 'grey55', linewidth = 0.3) +
    geom_linerange(aes(xmin = lo, xmax = hi), colour = 'grey55', linewidth = 0.35) +
    geom_point(aes(fill = class), shape = 21, colour = 'grey20', size = 1.5, stroke = 0.25) +
    scale_fill_manual(values = class_cols, guide = 'none') +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    expand_limits(x = 0) +
    labs(x = 'VIP\n(95% bootstrap CI)', y = NULL) + theme_fig() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_line(colour = 'grey93', linewidth = 0.25))
  if (markVIPsig && any(df_vip$sig))
    p_vip <- p_vip + geom_text(data = filter(df_vip, sig), aes(x = hi + 0.04),
                               label = '*', size = 1.4 * fs / .pt, vjust = 0.75)
  
  ## ---- (a) age-effect profile ----
  if (showAgeBar) {
    summ <- tibble(rr = as.vector(r$linr), q = as.vector(r$linq), x = rep(xv, each = nParc)) %>%
      group_by(x) %>%
      dplyr::summarise(pos = 100 * mean(q < alphaFDR & rr > 0),
                       neg = -100 * mean(q < alphaFDR & rr < 0), .groups = 'drop')
    bands_df <- if (showBands && cf$bands)
      tibble(lab = c('δ', 'θ', 'α', 'β', 'γ'), lo = c(x_lo, 4, 8, 13, 30),
             hi = c(4, 8, 13, 30, x_hi)) %>% filter(lo < x_hi) %>%
      mutate(hi = pmin(hi, x_hi), mid = (lo + hi) / 2) else NULL
    
    p_bar <- ggplot(summ, aes(x)) +
      geom_col(aes(y = pos), fill = col_pos, width = 0.9 * xstep) +
      geom_col(aes(y = neg), fill = col_neg, width = 0.9 * xstep) +
      geom_hline(yintercept = 0, colour = 'grey30', linewidth = 0.3) +
      annotate('text', x = x_hi + 0.015 * xspan, y = c(50, -50), hjust = 0, lineheight = 0.9,
               label = c('increase\nwith age', 'decrease\nwith age'),
               colour = c(col_pos, col_neg), size = fs / .pt, fontface = 'bold') +
      scale_x_continuous(breaks = xbreaks, expand = c(0, 0),
                         sec.axis = if (!is.null(bands_df))
                           dup_axis(breaks = bands_df$mid, labels = bands_df$lab, name = NULL) else waiver()) +
      scale_y_continuous(labels = abs, breaks = seq(-100, 100, 50)) +
      coord_cartesian(xlim = c(x_left, x_hi), ylim = c(-100, 100), expand = FALSE, clip = 'off') +
      labs(x = NULL, y = '% parcels\n(FDR < .05)') + theme_fig() +
      theme(axis.line.x.top = element_blank(), axis.ticks.x.top = element_blank(),
            axis.text.x.top = element_text(size = fs + 1))
    
    design <- c(area(1, 1), area(2, 1), area(2, 2))
    fig <- p_bar + p_beta + p_vip +
      plot_layout(design = design, widths = c(4, 1),
                  heights = c(top_h, row_h * nMaps), guides = 'collect')
    h <- top_h + row_h * nMaps + leg_h + 0.35
  } else {
    fig <- p_beta + p_vip + plot_layout(widths = c(4, 1), guides = 'collect')
    h <- row_h * nMaps + leg_h + 0.35
  }
  fig <- fig + plot_annotation(tag_levels = 'A') &
    theme(legend.position = 'bottom', legend.box = 'horizontal',
          legend.justification = 'center', legend.box.just = 'top',
          legend.box.spacing = unit(2, 'mm'),
          legend.key.size = unit(0.28, 'cm'),
          legend.spacing.x = unit(0.8, 'mm'),
          legend.text = element_text(size = fs - 0.5, margin = margin(l = 1, r = 4)),
          legend.margin = margin(0, 2.5, 0, 2.5, unit = 'mm'))   # smaller gap left/right of each legend
  save_fig(paste0(prefix, '_FigS_PLSbeta_VIP_linear'), fig, page_w, h)
  
  cat(sprintf('[%s] saved %.1f x %.1f in | coef range %.3f to %.3f | VIP > 1: %d | CI lower > 1: %d\n',
              prefix, page_w, h, min(r$beta), max(r$beta), sum(r$vip > 1), sum(r$ci[, 1] > 1)))
  invisible(fig)
}

for (k in metrics) make_fig(k)
