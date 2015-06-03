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
    0, -0.03,0,-0.04,
    0, 0,0,-0.05),
  label_fudge_y = c(
    0, -0.02,0.03,0.04,
    0, 0.02,0,-0.02))

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

d2_long <- reshape2::melt(d2, id.vars = c("test_iter", "clean_method", "type", "order",
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

xlim <- filter(d3, variable == "mare") %>% select(l, u) %>% range +
  c(-0.01, 0.01)
ylim <- filter(d3, variable == "corr") %>% select(l, u) %>% range +
  c(-0.02, 0.01)

pdf("../figs/fig3.pdf", width = 4, height = 3.1)
par(mfrow = c(1, 1), mgp = c(1.5, 0.4, 0), las = 1, tck = -0.012,
  oma = c(2.7, 3.5, .5, .5), cex = 0.8, mar = c(0, 0, 0, 0),
  xaxs = "i", yaxs = "i", col.axis = "grey50", col.lab = "grey50")
par(xpd = NA)
plot(m$mare, m$corr, xlim = xlim, ylim = ylim, pch = 21, bg = m$col, col = "grey50",
  axes = FALSE, xlab = "", ylab = "", type = "n")
segments(m$mare, l$corr, m$mare, u$corr, col = "#00000050", lwd = 1.4)
segments(l$mare, m$corr, u$mare, m$corr, col = "#00000050", lwd = 1.4)
points(m$mare, m$corr, pch = 21, bg = m$col, col = "grey50", cex = 1.5,
  lwd = 1.4)
xtext <- m$mare + fudge_x$mare
ytext <- m$corr + fudge_y$corr - 0.015
text(xtext, ytext, m$clean_method, cex = 0.90, pos = 4, col = m$text_col,
  bg = "red")
bg_cols <- c(rep(NA, 4), rep("#00000020", 4))
rect(xtext + 0.007, ytext -0.015, xtext + strwidth(m$clean_method),
  ytext + 0.02, col = bg_cols, border = NA)
box(col = "grey50")
axis(2, col = "grey50", cex.axis = 1, at = seq(0.1, 0.6, 0.1))
par(mgp = par()[["mgp"]] + c(0, -0.25, 0))
axis(1, col = "grey50", cex.axis = 1, at = seq(0.3, 0.6, 0.1))
mtext("Within population inaccuracy (MAPE)", side = 1, line = 1.4, col = "grey40", cex = 0.8)
mtext("Across population correlation", side = 2, line = 2.2, col = "grey40", las = 0, cex = 0.8)
legend("topright", bty = "n",
  fill =
  c(pal[pal$mre < 0.21 & pal$mre > 0.19,]$col,
    "white",
    pal[pal$mre > -0.21 & pal$mre < -0.19,]$col), border = rep("grey50", 3),
  legend = c("  0.2", "  0", "-0.2"), text.col = "grey40", title = "Bias (MPE)")
#text(0.54, 0.5, "Bias (MPE)", col = "grey30", pos = 4)
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
add_label <- function(xfrac, yfrac, label, pos = 4, add_bg = FALSE, ...) {

  u_ <- par("usr")
  width <- u_[2] - u_[1]
  height <- u_[4] - u_[3]
  x <- u_[1] + xfrac * width
  y <- u_[4] - yfrac * height
  text(x, y, label, pos = pos, ...)

  if (add_bg) {
    bg_col <- "#00000015"
    rect(x + width*0.04,                   y - height*0.04,
         x + strwidth(label) + width*0.07, y + strheight(label) - height * 0.01,
         col = bg_col, border = NA)
  }
}


plot_hex_fig <- function(dat, xbins = 100L, xlab = expression(B/B[MSY]),
  ylab = expression(widehat(B/B[MSY])), lims_hex = c(0, max(dat$bbmsy_est)),
  xlim_plot = c(0, 3.4), ylim_plot = c(0, 3.4), axis_ticks = c(0, 1, 2, 3),
  add_hex = TRUE, alpha = 50) {
  par(mfrow = c(2, 4), mgp = c(1.5, 0.5, 0), las = 1, tck = -0.03,
    oma = c(3.5, 3.5, .5, .5), cex = 0.8, mar = c(0, 0, 0, 0),
    xaxs = "i", yaxs = "i", col.axis = "grey50", col.lab = "grey50")

  hexcol1 <- RColorBrewer::brewer.pal(9, "Blues")
  #hexcol1 <- c("#FFFFFF", "#6A4A3C")
  hexcol2 <- RColorBrewer::brewer.pal(9, "Oranges")
  hexcol2 <- RColorBrewer::brewer.pal(9, "YlOrRd")
  #hexcol2 <- c("#FFFFFF", "#CC333F")
  hexcol3 <- RColorBrewer::brewer.pal(9, "Greens")
  hexcol4 <- RColorBrewer::brewer.pal(9, "Purples")
  hexcol_ensemble <- RColorBrewer::brewer.pal(9, "Greys")

  hexcol1 <- hexcol2
  hexcol3 <- hexcol2
  hexcol4 <- hexcol2
  hexcol_ensemble <- RColorBrewer::brewer.pal(9, "YlGnBu")

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
      pal_function <- colorRampPalette(hexcol1)
    if (m %in% 2)
      pal_function <- colorRampPalette(hexcol2)
    if (m %in% 3)
      pal_function <- colorRampPalette(hexcol3)
    if (m %in% 4)
      pal_function <- colorRampPalette(hexcol4)
    if (m %in% 5:8)
      pal_function <- colorRampPalette(hexcol_ensemble)
    if (add_hex) {
      pal <- pal_function(max(counts))
    } else {
      pal <- pal_function(10)
    }
    #add transparency to de-emphasize fist colour bins:
    pal[1:2] <- paste0(pal[1:2], round(seq(80, 99, length.out = 2)))
    plot(1, 1, xlim = xlim_plot, ylim = ylim_plot, type = "n", asp = 1,
      xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i")
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
    add_label(-0.01, 0.08, unique(filter(dat, order == m)$clean_method),
      col = "grey20", add_bg = ifelse(m %in% 1:4, FALSE, TRUE))
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

pdf("../figs/fig2.pdf", width = 7, height = 3.9)
plot_hex_fig(d_mean_plot, xbins = 100L)
dev.off()

d_slope_plot <- filter(d, type == "slope")
d_slope_plot$bbmsy_est[d_slope_plot$bbmsy_est > 10] <- NA
d_slope_plot <- na.omit(d_slope_plot)

pdf("../figs/hex-slope-sim.pdf", width = 7, height = 3.9)
plot_hex_fig(d_slope_plot, xbins = 80L,
  lims_hex = range(c(d_slope_plot$bbmsy_est, d_slope_plot$bbmsy_true)),
  xlim_plot = c(-0.5, 0.5), ylim_plot = c(-0.5, 0.5), axis_ticks = c(-0.4, 0, 0.4),
  xlab = expression(B/B[BMSY]~slope), ylab = expression(widehat(B/B[BMSY])~slope))
dev.off()

# plot the ram extrapolations:
d_ram <- readRDS("generated-data/ram-ensemble-predicted.rds")
clean_names_ram <- clean_names
# clean_names_ram$clean_method <- sub("LM Ensemble", "GAM Ensemble", clean_names_ram$clean_method)
# clean_names_ram$method <- sub("lm_ensemble", "gam_ensemble", clean_names_ram$method)
# clean_names_ram$clean_method <- sub("GBM Ensemble", "GAM Ensemble", clean_names_ram$clean_method)
# clean_names_ram$method <- sub("gbm_ensemble", "gam_ensemble", clean_names_ram$method)
d_ram <- suppressWarnings(inner_join(d_ram, clean_names_ram))
d_ram$bbmsy_est <- as.numeric(as.character(d_ram$bbmsy_est))

d_ram <- filter(d_ram, bbmsy_true < 4, bbmsy_est < 4)

pdf("../figs/hex-mean-ram-cv.pdf", width = 8, height = 4)
plot_hex_fig(d_ram, add_hex = TRUE, alpha = 80,
  lims_hex = range(c(d_ram$bbmsy_est, d_ram$bbmsy_true)), xbins = 25L)
dev.off()

re_ram <- d_ram %>% mutate(
  sq_er = (bbmsy_est - bbmsy_true)^2,
  re    = (bbmsy_est - bbmsy_true) / bbmsy_true)

re_ram_sum <- re_ram %>% group_by(clean_method) %>%
  summarise(mare = median(abs(re)),
  msqe = mean(sq_er),
  mre = median(re),
  corr = cor(bbmsy_true, bbmsy_est, method = "spearman")) %>%
  as.data.frame()

p <- re_ram_sum %>% #filter(!clean_method %in% "LM Ensemble") %>%
  ggplot(aes(mare, corr)) + geom_point(aes(colour = mre), size = 6) +
  geom_text(aes(label = clean_method), size = 3) +
  scale_colour_gradient2(guide = guide_legend(title = "Bias\n(MRE)")) +
  theme_bw() + xlab("Inaccuracy (MARE)") + ylab("Rank-order correlation") +
  xlim(range(re_ram_sum$mare) + c(-0.03, 0.01))
ggsave("../figs/ram-ensemble-performance-cv.pdf", width = 7, height = 5)

# ----------------------------------------------
# Example time series plot to motivate the study
# ----------------------------------------------

# try: SBWHITACIR
ram <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE) %>%
  filter(stockid == "SBWHITACIR") %>%
  select(tsyear, Bbmsy_toUse) %>%
  rename(yr = tsyear, b_bmsy = Bbmsy_toUse)

load("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/CMSY/cmsy_rlegacy_results_table_v0.RData")
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
  filter(stock_id == "SBWHITACIR") %>%
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

####################
# now look at above/below b/bmsy=1 classification success
####################

# # continuous case:
# abovebelow_sim <- d %>% group_by(clean_method, order) %>%
#   summarise(
#     true_above =  sum(bbmsy_est >  1 & bbmsy_true >  1) / length(bbmsy_est),
#     true_below =  sum(bbmsy_est <= 1 & bbmsy_true <= 1) / length(bbmsy_est),
#     false_below = sum(bbmsy_est <= 1 & bbmsy_true >  1) / length(bbmsy_est),
#     false_above = sum(bbmsy_est >  1 & bbmsy_true <= 1) / length(bbmsy_est)) %>%
#   as.data.frame()

# we're in trouble if all rows don't add to 1:
# all(rowSums(abovebelow_sim[,-c(1, 2)]) == 1) %>% stopifnot()
#
# abovebelow_sim$Ensemble <- ifelse(grepl("Ensemble", abovebelow_sim$clean_method),
#   TRUE, FALSE)
# abovebelow_sim <- abovebelow_sim %>%
#   mutate(clean_method = paste(order, clean_method, sep = "-"))
#
# p <- ggplot(abovebelow_sim) +
#   geom_segment(aes(x = true_above, xend = 0 , y = 0, yend = 0), lwd = 1.5, col = "darkgreen") +
#   geom_segment(aes(x = 0, xend = -false_above , y = 0, yend = 0), lwd = 1.5, col = "red") +
#   geom_segment(aes(x = 0, xend = 0 , y = true_below, yend = 0), lwd = 1.5, col = "darkgreen") +
#   geom_segment(aes(x = 0, xend = 0 , y = 0, yend = -false_below), lwd = 1.5, col = "red") +
#   facet_wrap(~clean_method, nrow = 2) + xlab("Above") + ylab("Below") + theme_bw()
# print(p)

get_roc <- function(true, est) {
  library("pROC") # required or pROC:: will generate errors
  y <- pROC::roc(response = true, predictor = est, 5)
  data.frame(sens = y$sensitivities, spec = y$specificities,
    auc = as.numeric(y$auc))
}

dbin <- readRDS("generated-data/cv_sim_binary.rds")
clean_names_bin <- clean_names
clean_names_bin$clean_method <- sub("LM Ensemble", "GLM Ensemble",
  clean_names_bin$clean_method)
clean_names_bin$method <- sub("lm_ensemble", "glm_ensemble", clean_names_bin$method)
dbin <- suppressWarnings(inner_join(dbin, clean_names_bin))

rocs_sim <- dbin %>% group_by(clean_method, order) %>%
  do({get_roc(true = .$bbmsy_true, est = .$bbmsy_est)})
rocs_sim <- rocs_sim %>% as.data.frame() %>%
  mutate(clean_method = paste(order, clean_method, sep = "-"))
rocs_sim$Ensemble <- ifelse(grepl("Ensemble", rocs_sim$clean_method),
  "Ensemble", "Individual")

p <- ggplot(rocs_sim, aes(spec, sens, colour = clean_method, group = clean_method)) +
  geom_line() + coord_equal() +
  geom_abline(intercept = 1, slope = 1, lty = 2, col = "darkgrey") +
  xlim(1, 0) +
  xlab("Specificity") + ylab("Sensitivity") + theme_bw() +
  theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + facet_wrap(~Ensemble) +
  scale_colour_brewer(palette="Spectral", guide = guide_legend(title = "Model"))
  print(p)
ggsave("../figs/roc-sim.pdf", width = 8, height = 5)
