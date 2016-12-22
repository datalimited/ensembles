# General purpose functions to be used throughout

train_spec_mat <- function(x, freq_vec = 1/c(2, 5, 10, 20)) {
  # using AR as smoother, empirical didn't seem to confer much more benefit
  sp <- spec.ar(x/max(x), plot = FALSE)
  # approximate at fixed frequencies - necessary as series of different length
  approx(x = sp$freq, y = sp$spec, xout = freq_vec) %>% as.data.frame
}

mean_slope_bbmsy <- function(dat, years_window = 5L) {
  # chunk of data must have columns: b_bmsy_true, b_bmsy_est
  # message(paste(unique(dat$stockid), unique(dat$iter), sep = "-"))
  if (ncol(dat) > 0) { # not sure what's happening here, but some chunks can have zero columns
    if (nrow(dat) > years_window) { # some have 3 years??
      .n <- nrow(dat)
      i <- seq(.n-(years_window-1), .n)
      bbmsy_true_mean = mean(dat$b_bmsy_true[i])
      bbmsy_est_mean = mean(dat$b_bmsy_est[i])
      ytrue <- dat$b_bmsy_true[i]
      yest <- dat$b_bmsy_est[i]
      bbmsy_true_slope <- coef(mblm::mblm(ytrue ~ i))[[2]]
      bbmsy_est_slope <- coef(mblm::mblm(yest ~ i))[[2]]
      data.frame(bbmsy_true_mean, bbmsy_est_mean, bbmsy_true_slope, bbmsy_est_slope)
    }
  }
}

# A general function for cross-validation testing ensemble models:
cross_val_ensembles <- function(.n, dat, fraction_train = 0.5,
  gbm_formula = "log(bbmsy_true_mean) ~ CMSY + COMSIR + mPRM + SSCOM + LH",
  geo_mean = TRUE,
  individual_models = c("CMSY", "COMSIR", "mPRM", "SSCOM"),
  id = "mean", cache_folder = "generated-data/cv-sim/", distribution = "gaussian",
  lm_formula = "", glm_formula = "", nfold = 3L, weighted = FALSE) {

  ids <- unique(dat$cv_id)
  ids_scrambled <- base::sample(ids)
  chunk_size <- floor(length(ids) / nfold)
  chunk_starts <- seq(1L, chunk_size * nfold, chunk_size)
  # excluding those that aren't used in this sorting due to rounding:
  ids_scrambled_cut <- ids_scrambled[seq_len(chunk_size * nfold)]

  out <- plyr::ldply(seq_along(chunk_starts), function(i) {

    test_ids  <- ids_scrambled_cut[seq(chunk_starts[i], chunk_starts[i] + chunk_size)]
    train_ids <- ids_scrambled_cut[!ids_scrambled_cut %in% test_ids]
    train_dat <- dplyr::filter(dat, cv_id %in% train_ids)
    test_dat  <- dplyr::filter(dat, cv_id %in% test_ids)

    if (weighted) {
      # with weights inverse to true B/BMSY to focus on low values:
      m_gbm <- gbm::gbm(as.formula(gbm_formula),
        data = train_dat, n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01,
        distribution = distribution, weights = 1/train_dat$bbmsy_true_mean)
    } else {
      m_gbm <- gbm::gbm(as.formula(gbm_formula),
        data = train_dat, n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01,
        distribution = distribution)
    }
    test_dat$gbm_ensemble <- tryCatch({gbm::predict.gbm(m_gbm,
      n.trees = m_gbm$n.trees, newdata = test_dat, type = "response")},
      error = function(e) rep(NA, nrow(test_dat)))

    m_rf <- randomForest::randomForest(as.formula(gbm_formula), data = train_dat,
      ntree = 1000L)
    test_dat$rf_ensemble <- tryCatch({predict(m_rf, newdata = test_dat)},
      error = function(e) rep(NA, nrow(test_dat)))

    if (lm_formula != "") {
      m_lm <- lm(as.formula(lm_formula), data = train_dat)
      test_dat$lm_ensemble <- tryCatch({predict(m_lm, newdata = test_dat)},
        error = function(e) rep(NA, nrow(test_dat)))
    }

    # if (lm_formula != "") {
    # library("mgcv")
    # m_gam <- mgcv::gam(as.formula(lm_formula), data = train_dat)
    # test_dat$lm_ensemble <- tryCatch({predict(m_lm, newdata = test_dat)},
    # error = function(e) rep(NA, nrow(test_dat)))
    # }

    if (glm_formula != "") {
      m_glm <- glm(as.formula(glm_formula), data = train_dat, family = binomial(link = logit))
      test_dat$glm_ensemble <- tryCatch({predict(m_glm, newdata = test_dat, type = "response")},
        error = function(e) rep(NA, nrow(test_dat)))
    }

    if (geo_mean) {
      test_dat$mean_ensemble <- exp(rowMeans(log(test_dat[, individual_models])))
    } else {
      test_dat$mean_ensemble <- rowMeans(test_dat[, individual_models])
    }

    test_dat$test_iter <- .n
    test_dat
  })
  out
}

add_label <- function(xfrac, yfrac, label, pos = 4, add_bg = FALSE, ...) {
  # xfrac The fraction over from the left side.
  # yfrac The fraction down from the top.
  # label The text to label with.
  # pos Position to pass to text()
  # ... Anything extra to pass to text(), e.g. cex, col.
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
  xlim_plot = c(0, 3.9), ylim_plot = c(0, 3.3), axis_ticks = c(0, 1, 2, 3),
  add_hex = TRUE, alpha = 30, xbins3 = xbins, lims_hex3 = lims_hex,
  count_transform = 1.5, count_transform3 = 30, oma = c(3.5, 3.5, .5, .5),
  bias = c(rep(3, 12))) {

  rows <- max(dat$order) / 4
  par(mfrow = c(rows, 4), mgp = c(1.5, 0.5, 0), las = 1, tck = -0.03,
    oma = oma, cex = 0.8, mar = c(0, 0, 0.27, 0),
    xaxs = "i", yaxs = "i", col.axis = "grey60", col.lab = "grey50")

  hexcol1 <- RColorBrewer::brewer.pal(9, "Blues")
  hexcol2 <- RColorBrewer::brewer.pal(9, "Oranges")
  hexcol2 <- RColorBrewer::brewer.pal(9, "YlOrRd")
  hexcol3 <- RColorBrewer::brewer.pal(9, "Greens")
  hexcol4 <- RColorBrewer::brewer.pal(9, "Purples")
  hexcol_ensemble <- RColorBrewer::brewer.pal(9, "Greys")
  hexcol_ensemble <- RColorBrewer::brewer.pal(9, "YlGnBu")

  hexcol1 <- hexcol2
  hexcol3 <- hexcol2
  hexcol4 <- hexcol2
  library("RColorBrewer")
  hexcol_third_row <- RColorBrewer::brewer.pal(9, "YlGnBu")

  panels <- seq_len(length(unique(dat$order)))

  plyr::l_ply(panels, function(m) {

    if (add_hex) {
      if (m <= 8) {
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
        xbins_plot <- xbins
      } else { # the optional third row:
        xlim <- lims_hex3
        ylim <- lims_hex3
        bin <- hexbin::hexbin(
          filter(dat, order == m)$bbmsy_true,
          filter(dat, order == m)$bbmsy_est,
          xbnds = xlim, ybnds = ylim, xbins = xbins3)
        dx <- hexbin::hcell2xy(bin)$x
        dy <- hexbin::hcell2xy(bin)$y
        dxy <- data.frame(x = dx, y = dy)
        counts <- bin@count
        xbins_plot <- xbins3
      }
    }
    space = "Lab"
    if (m %in% 1)
      pal_function <- colorRampPalette(hexcol1, space = space, bias = bias[1])
    if (m %in% 2)
      pal_function <- colorRampPalette(hexcol2, space = space, bias = bias[2])
    if (m %in% 3)
      pal_function <- colorRampPalette(hexcol3, space = space, bias = bias[3])
    if (m %in% 4)
      pal_function <- colorRampPalette(hexcol4, space = space, bias = bias[4])
    if (m %in% 5:8)
      pal_function <- colorRampPalette(hexcol_ensemble, space = space, bias = bias[5])
    if (m %in% 9:12) {
      pal_function <- colorRampPalette(hexcol_third_row, space = space, bias = bias[9])
    }
    if (add_hex) {
      pal <- pal_function(max(counts))
    } else {
      pal <- pal_function(10)
    }
    plot(1, 1, xlim = xlim_plot, ylim = ylim_plot, type = "n", asp = 1,
      xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i")
    if (add_hex) {
      # add transparency to de-emphasize fist colour bins:
      pal[1:5] <- paste0(pal[1:5], round(seq(80, 99, length.out = 5)))
      for (i in 1:nrow(dxy)) {
        hexagon(dxy[i, "x"], dxy[i, "y"], col = pal[counts[i]],
          unitcell = diff(xlim)/xbins_plot/2, border = NA)
      }
    } else {
      dd <- dat %>% filter(order == m)
      points(dd$bbmsy_true, dd$bbmsy_est, col = paste0(pal[8], alpha), pch = 21,
        bg = paste0(pal[4], alpha), cex = 0.7)
    }
    line_col <- "#58585840"
    segments(1, 0, 1, 2.8, lty = "22", col = line_col, lwd = 1.3)
    abline(h = 1, lty = "22", col = line_col, lwd = 1.3)
    segments(-1, -1, 2.8, 2.8, lty = "22", col = line_col, lwd = 1.3)
    box(col = "grey60")
    add_label(-0.01, 0.08, paste0("(", letters[m], ") ", unique(filter(dat, order == m)$clean_method)),
      col = "grey20", add_bg = FALSE, cex = 0.95)
    if (m %in% c(1, 5, 9)) axis(2, at = axis_ticks[-length(axis_ticks)], col = "grey60")
    if (m %in% (rows * 4 - 3):(rows * 4)) axis(1, at = axis_ticks, col = "grey60")
  })
  mtext(xlab, side = 1, line = 2.1, outer = TRUE, col = "grey20")
  mtext(ylab, side = 2, line = 1.6, outer = TRUE,
    las = 0, col = "grey20")
}

hexagon <- function (x, y, unitcell = 1, ...) {
  polygon(
    hexbin::hexcoords(unitcell)$x + x,
    hexbin::hexcoords(unitcell)$y + y, ...)
}
