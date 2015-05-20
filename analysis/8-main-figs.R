# Make the main figures
library(dplyr)
library(ggplot2)

d <- readRDS("generated-data/cv_sim_long.rds") %>%
  filter(method != "dummy") %>%
  as.data.frame()

clean_names <- dplyr::data_frame(
  method = c("CMSY", "COMSIR", "Costello", "SSCOM",
    "gbm_ensemble", "rf_ensemble", "lm_ensemble", "mean_ensemble"),
  clean_method = c("CMSY", "COM-SIR", "mPRM", "SSCOM",
    "GBM Ensemble", "RF Ensemble", "LM Ensemble", "Mean Ensemble"),
  order = c(1, 2, 4, 3, 8, 7, 6, 5),
  label_fudge_x = c(
    0, -0.06,0,-0.06,
    0, 0,0,0),
  label_fudge_y = c(
    0, 0,0.03,0,
    0, 0,0,0))

d <- suppressWarnings(inner_join(d, clean_names))

d2 <- d %>% group_by(test_iter, type, clean_method, order, label_fudge_x, label_fudge_y) %>%
  summarise(
    mare = median(abs(re)),
    mre = median(re),
    corr = cor(bbmsy_true_trans, bbmsy_est_trans, method = "spearman",
      use = "pairwise.complete.obs"))

d2$clean_method <- as.factor(d2$clean_method)
d2$clean_method <- reorder(d2$clean_method, d2$order)

d2_slope <- filter(d2, type == "slope")
d2 <- filter(d2, type == "mean")

# ------------------------------------
# Performance distributions simulation
# ------------------------------------

# ggplot(d2, aes(clean_method, corr)) + geom_violin()
# ggplot(d2, aes(clean_method, mare)) + geom_violin()
# ggplot(d2, aes(clean_method, mre)) + geom_violin()
# d %>% filter(re < 4, re > -4) %>%
#   ggplot(aes(re)) + geom_histogram() + facet_wrap(~clean_method) +
#   geom_vline(xintercept = 0)

#bg_plot <- function(colour = "#00000009") {
  #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
    #col = colour, border = FALSE)
#}

performance_panel <- function(dat, column_id, ylim_adj = c(-0.1, 0.1), xaxis = FALSE,
  flip_yaxis = FALSE, yticks = NULL, add_zero_line = TRUE) {

  axis_col <- "grey50"
  label_col <- "grey50"
  cols_ii <- RColorBrewer::brewer.pal(8, "Blues")

  ylim <- range(dat[,column_id][[1]]) + ylim_adj
  if (flip_yaxis) ylim <- rev(ylim)

  plot(1, 1, xlim = c(.8, nrow(clean_names)+0.2),
    ylim = ylim, type = "n", axes = FALSE, ann = FALSE, yaxs = "i")
  rect(4.5, ylim[1], 9.5, ylim[2], col = "#00000010", border = FALSE)
  rect(5.5, ylim[1], 9.5, ylim[2], col = "#00000010", border = FALSE)
  rect(6.5, ylim[1], 9.5, ylim[2], col = "#00000010", border = FALSE)
  rect(7.5, ylim[1], 9.5, ylim[2], col = "#00000010", border = FALSE)
  if(add_zero_line)
    abline(h = 0, col = "grey70", lty = 2, lwd = 1)
  beanplot::beanplot(as.formula(paste(column_id, "~ clean_method")),
    data = dat, add = TRUE,
    border = NA, axes = FALSE, col = cols_ii[c(5, 3, 4, 8)], what =
    c(0, 1, 1, 0), log = FALSE)
  points(jitter(as.numeric(dat$clean_method), amount = 0.09),
    dat[,column_id][[1]], col = "#D0D0D0", pch = 20, cex = 0.37)
  if (is.null(yticks))
    axis(2, las = 1, col = axis_col, col.axis = axis_col)
  else
    axis(2, las = 1, at = yticks, col = axis_col, col.axis = axis_col)

  if (xaxis)
    axis(1, col = axis_col, col.axis = label_col, at = 1:length(unique(dat$clean_method)),
      labels = levels(dat$clean_method), las = 3)
  box(col = axis_col)

}

pdf("../figs/performance-beanplots-sim.pdf", width = 5, height = 5)
par(mfrow = c(3, 1))
par(mar = c(0, 3, 0, 0), cex = 0.8, oma = c(8, 3, 1.5, .5), tck = -0.03, mgp = c(2, 0.5, 0))
performance_panel(d2, "corr", c(-0.1, 0.1), flip_yaxis = FALSE, yticks = c(0, 0.2, 0.4))
mtext("Correlation\nacross populations", side = 2, line = 3, col = "grey40")
performance_panel(d2, "mare", c(-0.05, 0.05), xaxis = FALSE, flip_yaxis = FALSE,
  yticks = c(0.3, 0.4, 0.5, 0.6))
mtext("Inaccuracy\n(MARE)", side = 2, line = 3, col = "grey40")
performance_panel(d2, "mre", c(-0.1, 0.1), flip_yaxis = FALSE, xaxis = TRUE,
  yticks = c(0, 0.2, 0.4))
mtext("Bias and\nprecision (MRE)", side = 2, line = 3, col = "grey40")
dev.off()

d2_long <- reshape2::melt(d2, id.vars = c("test_iter", "clean_method", "order",
    "label_fudge_x", "label_fudge_y"))

d3 <- d2_long %>% group_by(clean_method, order, variable, label_fudge_x, label_fudge_y) %>%
  summarise(
    m = median(value),
    l = quantile(value, 0.25),
    u = quantile(value, 0.75)) %>%
  as.data.frame()

d3$text_col <- ifelse(grepl("Ensemble", d3$clean_method), "grey20", "grey50")

l <- reshape2::dcast(d3, clean_method ~ variable, value.var = "l")
u <- reshape2::dcast(d3, clean_method ~ variable, value.var = "u")
m <- reshape2::dcast(d3, clean_method + text_col ~ variable, value.var = "m")

fudge_x <- reshape2::dcast(d3, clean_method ~ variable, value.var = "label_fudge_x")
fudge_y <- reshape2::dcast(d3, clean_method ~ variable, value.var = "label_fudge_y")

pal <- data_frame(
  col = RColorBrewer::brewer.pal(11, "RdBu"),
  mre = seq(-max(m$mre)*1.05, max(m$mre)*1.05, length.out = 11))
m$col <- pal$col[findInterval(m$mre, pal$mre)]

xlim <- filter(d3, variable == "mare") %>% select(l, u) %>% range
ylim <- filter(d3, variable == "corr") %>% select(l, u) %>% range

pdf("../figs/fig3.pdf", width = 4, height = 3.1)
par(mfrow = c(1, 1), mgp = c(1.5, 0.4, 0), las = 1, tck = -0.012,
  oma = c(3, 3.5, .5, .5), cex = 0.8, mar = c(0, 0, 0, 0),
  xaxs = "i", yaxs = "i", col.axis = "grey50", col.lab = "grey50")
par(xpd = NA)
plot(m$mare, m$corr, xlim = xlim, ylim = ylim, pch = 21, bg = m$col, col = "grey50",
  axes = FALSE, xlab = "", ylab = "", type = "n")
segments(m$mare, l$corr, m$mare, u$corr, col = "#00000050", lwd = 1.4)
segments(l$mare, m$corr, u$mare, m$corr, col = "#00000050", lwd = 1.4)
points(m$mare, m$corr, pch = 21, bg = m$col, col = "grey50", cex = 1.5)
text(
  m$mare + fudge_x$mare,
  m$corr + fudge_y$corr - 0.015,
  m$clean_method, cex = 0.90, pos = 4, col = m$text_col)
box(col = "grey50")
axis(2, col = "grey50", cex.axis = par()[["cex"]], at = seq(0.1, 0.5, 0.1))
par(mgp = par()[["mgp"]] + c(0, -0.25, 0))
axis(1, col = "grey50", cex.axis = par()[["cex"]], at = seq(0.3, 0.6, 0.1))
mtext("Within population inaccuracy (MARE)", side = 1, line = 1.7, col = "grey40", cex = 0.8)
mtext("Across population correlation", side = 2, line = 2.2, col = "grey40", las = 0, cex = 0.8)
#legend("bottomright", bty = "n", pch = c(21, 21, 21), bg = c("red", "red", "blue"),
  #legend = c(4, 0, -1))
dev.off()

#mtext(LETTERS[ii], adj = 0.05, line = -1.5, col = "grey40", cex = 0.8)
#text(0.5, 0.48, panel_labs[ii-4], col = label_col, pos = 4, cex = 1.05)


# ------------------------
# Scatter plots simulation
# ------------------------
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

plot_hex_fig <- function(dat, xbins = 100L, xlab = expression(B/B[MSY]),
  ylab = expression(widehat(B/B[MSY])), lims_hex = c(0, max(dat$bbmsy_est)),
  xlim_plot = c(0, 2.5), ylim_plot = c(0, 2.7), axis_ticks = c(0, 1, 2),
  add_hex = TRUE, alpha = 50) {
  par(mfrow = c(2, 4), mgp = c(1.5, 0.5, 0), las = 1, tck = -0.03,
    oma = c(3.5, 3.5, .5, .5), cex = 0.8, mar = c(0, 0, 0, 0),
    xaxs = "i", yaxs = "i", col.axis = "grey50", col.lab = "grey50")

  #if (max(dat$order) == 8)
    #panels <- c(1:4, 8, 7, 6, 5)
  #else
    panels <- seq_len(length(unique(dat$order)))

  plyr::l_ply(panels, function(m) {

    if (add_hex) {
      xlim <- lims_hex
      ylim <- lims_hex
      bin <- hexbin::hexbin(
        filter(dat, order == m)$bbmsy_true,
        filter(dat, order == m)$bbmsy_est,
        xbnds = xlim, ybnds = ylim, xbins = xbins)
      dx <- hexbin::hcell2xy(bin)$x
      dy <- hexbin::hcell2xy(bin)$y
      dxy <- data.frame(x = dx, y = dy)
      counts <- bin@count
      counts <- round(log(counts*1.5))
    }
    #alternative power transformation:
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
    if (add_hex) {
      pal <- pal_function(max(counts))
    } else {
      pal <- pal_function(10)
    }
    #add transparency to de-emphasize fist colour bins:
    #pal[1:10] <- paste0(pal[1:10], round(seq(50, 99, length.out = 10)))
    plot(1, 1, xlim = xlim_plot, ylim = ylim_plot, type = "n", asp = 1,
      xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    if (add_hex) {
      for (i in 1:nrow(dxy)) {
        hexagon(dxy[i, "x"], dxy[i, "y"], col = pal[counts[i]],
          unitcell = diff(xlim)/xbins/2, border = NA)
      }
    } else {
      dd <- dat %>% filter(order == m)
      points(dd$bbmsy_true, dd$bbmsy_est, col = paste0(pal[8], alpha), pch = 21,
        bg = paste0(pal[5], alpha))
    }
    abline(v = 1, lty = "22", col = "#33333350", lwd = 1.5)
    abline(h = 1, lty = "22", col = "#33333350", lwd = 1.5)
    abline(a = 0, b = 1, lty = "22", col = "#33333350", lwd = 1.5)
    box(col = "grey50")
    add_label(-0.01, 0.08, unique(filter(dat, order == m)$clean_method), col = "grey20")
    if (m %in% c(1, 5)) axis(2, at = axis_ticks, col = "grey50")
    if (m %in% 5:8) axis(1, at = axis_ticks, col = "grey50")
    })
  mtext(xlab, side = 1, line = 2.1, outer = TRUE, col = "grey20")
  mtext(ylab, side = 2, line = 1.6, outer = TRUE,
    las = 0, col = "grey20")
}


d_mean_plot <- filter(d, type == "mean")
d_mean_plot$bbmsy_est[d_mean_plot$bbmsy_est > 10] <- NA
d_mean_plot <- na.omit(d_mean_plot)

pdf("../figs/fig2.pdf", width = 8, height = 4)
plot_hex_fig(d_mean_plot, xbins = 100L)
dev.off()

d_slope_plot <- filter(d, type == "slope")
d_slope_plot$bbmsy_est[d_slope_plot$bbmsy_est > 10] <- NA
d_slope_plot <- na.omit(d_slope_plot)

pdf("../figs/hex-slope-sim.pdf", width = 8, height = 4)
plot_hex_fig(d_slope_plot, xbins = 80L, lims_hex = range(d_slope_plot$bbmsy_est),
  xlim_plot = c(-0.5, 0.5), ylim_plot = c(-0.5, 0.5), axis_ticks = c(-0.4, 0, 0.4),
  xlab = expression(B/B[BMSY]~slope), ylab = expression(widehat(B/B[BMSY])~slope))
dev.off()

# plot the ram extrapolations:
d_ram <- readRDS("generated-data/ram-ensemble-predicted.rds")
d_ram <- suppressWarnings(inner_join(d_ram, clean_names))

pdf("../figs/hex-mean-ram.pdf", width = 8, height = 4)
plot_hex_fig(d_ram, add_hex = FALSE, alpha = 80)
dev.off()

re_ram <- d_ram %>% mutate(
  sq_er = (bbmsy_est - bbmsy_true)^2,
  re    = (bbmsy_est - bbmsy_true) / bbmsy_true)

re_ram %>% group_by(clean_method) %>%
  summarise(mare = median(abs(re)),
  mre = median(re),
  corr = cor(bbmsy_true, bbmsy_est, method = "spearman"))

# ----------------------------------------------
# Example time series plot to motivate the study
# ----------------------------------------------

# try: SBWHITACIR
ram <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE) %>%
  filter(stockid == "SBWHITACIR") %>%
  select(tsyear, Bbmsy_toUse) %>%
  rename(yr = tsyear, b_bmsy = Bbmsy_toUse)

load("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/CMSY/cmsy_rlegacy_output_v0.RData")
cmsy <- cmsy.rlegacy.df0 %>%
  filter(stock_id == "SBWHITACIR") %>%
  select(yr, b_bmsyiq25, b_bmsy, b_bmsyiq75) %>%
  mutate(method = "CMSY")

load("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/COMSIR/comsir_rlegacy_results_table_v0.RData")
comsir <- comsir.rlegacy.df0 %>%
  filter(stock_id == "SBWHITACIR") %>%
  select(yr, b_bmsyiq25, b_bmsy, b_bmsyiq75) %>%
  mutate(method = "COM-SIR")

load("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/SSCOM/sscom_rlegacy_results_table_v0.RData")
sscom <- sscom.rlegacy.df0 %>%
  filter(stock_id == "") %>%
  select(yr, b_bmsy_iq25, b_bmsy, b_bmsy_iq75) %>%
  rename(b_bmsyiq25 = b_bmsy_iq25, b_bmsyiq75 = b_bmsy_iq75) %>%
  mutate(method = "SSCOM")

costello <- read.csv("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/Costello/Costello_rlegacy_results.csv")
cdat <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE)
costello <- merge(costello, unique(cdat[,c("stockid","stocklong")]), by="stocklong", all.x=TRUE)
costello <- costello %>% filter(stockid == "SBWHITACIR") %>%
  select(year, BvBmsy, BvBmsy_LogSD) %>%
  mutate(year, b_bmsyiq25 = exp(log(BvBmsy) - 0.5 * BvBmsy_LogSD),
    b_bmsy = BvBmsy,
    b_bmsyiq75 = exp(log(BvBmsy) + 0.5 * BvBmsy_LogSD)) %>%
  select(-BvBmsy, -BvBmsy_LogSD) %>%
  rename(yr = year) %>%
  mutate(method = "mPRM")

dl_dat <- bind_rows(cmsy, comsir, sscom, costello) %>%
  arrange(method, yr)

#library(ggplot2)
#ggplot(dl_dat, aes(yr, b_bmsy, colour = method)) + geom_line() + xlim(1990, 2011) +
  #geom_ribbon(aes(ymax = b_bmsyiq75, ymin = b_bmsyiq25, fill = method), alpha = 0.2) +
  #geom_line(data = ram, aes(yr, b_bmsy), colour = "black", fill = "black")


col_df <- data_frame(col = paste0(RColorBrewer::brewer.pal(4, "Dark2")),
  method = c("CMSY", "COM-SIR", "SSCOM", "mPRM"))
dl_dat$col <- NULL # for re-running
dl_dat <- inner_join(dl_dat, col_df)

plot_method <- function(dat) {
  polygon(c(dat$yr, rev(dat$yr)), c(dat$b_bmsyiq25, rev(dat$b_bmsyiq75)), border = NA,
    col = paste0(dat$col, 70))
  lines(dat$yr, dat$b_bmsy, col = dat$col, lwd = 2.5)
}

xlim <- c(1990, 2010)
ylim <- range(c(dl_dat$b_bmsyiq25, dl_dat$b_bmsyiq75))

pdf("../figs/motivate.pdf", width = 5, height = 3.5)
par(mar = c(3, 3.4, .5, 4.5), cex = 0.8, oma = c(0, 0, 0, 0), tck = -0.015,
  mgp = c(2, 0.5, 0), col.axis = "grey40", las = 1)
plot(1, 1, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)
abline(h = 1, lty = 2, col = "grey40", lwd = 2)
plyr::d_ply(dl_dat, "method", plot_method)
lines(ram$yr, ram$b_bmsy, lwd = 3.5, col = "grey20")
axis(2, col = "grey40", cex.axis = 0.9)
par(mgp = par()[["mgp"]] + c(0, -0.25, 0))
axis(1, col = "grey40", cex.axis = 0.9)
box(col = "grey40")
mtext("Year", side = 1, line = 1.7, col = "grey20", cex = 0.8)
mtext(expression(widehat(B/B[MSY])), side = 2, line = 1.7, col = "grey20", cex = 0.8,
  las = 0)

lab <- filter(dl_dat, yr == max(dl_dat$yr))
par(xpd = NA)
text(lab$yr+0, lab$b_bmsy, lab$method, pos = 4, col = lab$col, cex = 0.9)
text(unique(lab$yr)+0, 0.05 + ram[nrow(ram), "b_bmsy"][[1]], "Assessed", pos = 4, cex = 0.9)
text(xlim[1]-0.5, ylim[2]-0.1,
  paste0("Example stock:\n", unique(filter(cdat, stockid == "SBWHITACIR")$stocklong)),
  cex = 0.9, pos = 4, col = "grey20")

dev.off()
