# Make the main figures

library(dplyr)
library(ggplot2)

d <- readRDS("generated-data/cv_sim_long.rds") %>%
  filter(type == "mean") %>%
  filter(method != "dummy") %>%
  as.data.frame()

# d2 <- d %>% group_by(test_iter, method) %>%
#   summarise(
#     mare = median(abs(re)),
#     mre = median(re),
#     corr = cor(bbmsy_true_trans, bbmsy_est_trans, method = "spearman",
#       use = "pairwise.complete.obs"))
#
# ggplot(d2, aes(method, corr)) + geom_violin()
# ggplot(d2, aes(method, mare)) + geom_violin()
# ggplot(d2, aes(method, mre)) + geom_violin()
#
# d %>% filter(re < 4, re > -4) %>%
#   ggplot(aes(re)) + geom_histogram() + facet_wrap(~method) +
#   geom_vline(xintercept = 0)
#
# d %>% filter(re < 4, re > -4) %>%
#   ggplot(aes(method, re)) + geom_boxplot() +
#   geom_hline(yintercept = 0)
#
# p <- d %>% filter(method != "dummy") %>%
#   ggplot(aes(bbmsy_true, bbmsy_est)) +
#     geom_point(alpha = 0.01) +
#     facet_wrap(~method, ncol = 4) + ylim(0, 3) + xlim(0, 3) +
#     coord_fixed() +
#     geom_abline(intercept = 0, slope = 1, lty = 3)
# suppressWarnings(print(p))

hexagon <- function (x, y, unitcell = 1, ...) {
  polygon(
    hexbin::hexcoords(unitcell)$x + x,
    hexbin::hexcoords(unitcell)$y + y, ...)
}

# testing hexagon binning:
# x <- rnorm(2000)
# y <- rnorm(2000)
# xbins <- 10
# bin <- hexbin::hexbin(x, y, xbnds = c(-5, 5), ybnds = c(-5, 5), xbins = xbins, IDs = TRUE)
# plot(1, 1, xlim = c(-5, 5), ylim = c(-5, 5), type = "n", asp = 1)
# dx <- hcell2xy(bin)$x
# dy <- hcell2xy(bin)$y
# dxy <- data.frame(x = dx, y = dy)
# pal <- grey(seq(1, 0, length.out = max(bin@count)))
# for (i in 1:nrow(dxy)) {
#   hexagon(dxy[i, "x"], dxy[i, "y"], col = pal[bin@count[i]], unitcell = 0.5, border = NA)
# }

#' @param xfrac The fraction over from the left side.
#' @param yfrac The fraction down from the top.
#' @param label The text to label with.
#' @param pos Position to pass to text()
#' @param ... Anything extra to pass to text(), e.g. cex, col.
add_label <- function(xfrac, yfrac, label, pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

clean_names <- dplyr::data_frame(
  method = c("CMSY", "COMSIR", "Costello", "SSCOM",
    "gbm_ensemble", "rf_ensemble", "lm_ensemble", "mean_ensemble"),
  clean_method = c("CMSY", "COM-SIR", "mPRM", "SSCOM",
    "GBM Ensemble", "RF Ensemble", "LM Ensemble", "Mean Ensemble"))

d <- suppressWarnings(inner_join(d, clean_names))

pdf("figs/fig2.pdf", width = 8, height = 4)
par(mfrow = c(2, 4), mgp = c(1.5, 0.5, 0), las = 1, tck = -0.03, oma = c(3.5, 3.5, .5, .5),
  cex = 0.8, mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", col.axis = "grey50",
  col.lab = "grey50")
plyr::l_ply(c(1:4, 8, 7, 6, 5), function(m) {
  xbins <- 75
  xlim <- c(0, max(d$bbmsy_est))
  ylim <- c(0, max(d$bbmsy_est))
  bin <- hexbin::hexbin(
    filter(d, method == unique(d$method)[m])$bbmsy_true,
    filter(d, method == unique(d$method)[m])$bbmsy_est,
    xbnds = xlim, ybnds = ylim, xbins = xbins)
  dx <- hexbin::hcell2xy(bin)$x
  dy <- hexbin::hcell2xy(bin)$y
  dxy <- data.frame(x = dx, y = dy)
  counts <- bin@count
  counts <- round(log(counts*1.5))
  #counts <- round((counts)^0.30)
  if (m %in% 1)
    pal_function <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  if (m %in% 2)
    pal_function <- colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))
  if (m %in% 3)
    pal_function <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))
  if (m %in% 4)
    pal_function <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))
  if (m %in% 5:8)
    pal_function <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))
  #pal_function <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Spectral")))
  pal <- pal_function(max(counts))
  #pal[1:10] <- paste0(pal[1:10], round(seq(50, 99, length.out = 10)))
  plot(1, 1, xlim = c(0, 2.5), ylim = c(0, 2.7), type = "n", asp = 1, xlab = "", ylab = "",
    xaxt = "n", yaxt = "n")
  for (i in 1:nrow(dxy)) {
    hexagon(dxy[i, "x"], dxy[i, "y"], col = pal[counts[i]], unitcell = diff(xlim)/xbins/2,
      border = NA)
  }
  abline(v = 1, lty = "22", col = "#33333350", lwd = 1.5)
  abline(h = 1, lty = "22", col = "#33333350", lwd = 1.5)
  abline(a = 0, b = 1, lty = "22", col = "#33333350", lwd = 1.5)
  box(col = "grey50")
  add_label(-0.01, 0.08, unique(d$clean_method)[m], col = "grey20")
  if (m %in% c(1, 8)) axis(2, at = c(0, 1, 2), col = "grey50")
  if (m %in% 5:8) axis(1, at = c(0, 1, 2), col = "grey50")
    })
 mtext(expression(B/B[MSY]), side = 1, line = 2.1, outer = TRUE, col = "grey20")
 mtext(expression(widehat(B/B[MSY])), side = 2, line = 1.6, outer = TRUE, las = 0, col = "grey20")
dev.off()
