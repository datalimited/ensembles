# General purpose functions to be used throughout

train_spec_mat <- function(x, freq_vec = 1/c(2, 5, 10, 20)) {
  # using AR as smoother, empirical didn't seem to confer much more benefit
  sp <- spec.ar(x/max(x), plot = FALSE)
  # approximate at fixed frequencies - necessary as series of different length
  approx(x = sp$freq, y = sp$spec, xout = freq_vec) %>% as.data.frame
}

mean_slope_bbmsy <- function(dat, years_window = 3) {
  # chunk of data must have columns: b_bmsy_true, b_bmsy_est
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

# A general function for cross-validation testing ensemble models:
cross_val_ensembles <- function(.n, dat, fraction_train = 0.5,
  gbm_formula = "log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM + LH",
  geo_mean = TRUE,
  individual_models = c("CMSY", "COMSIR", "Costello", "SSCOM"),
  id = "mean", cache_folder = "generated-data/cv-sim/", distribution = "gaussian",
  lm_formula = "", glm_formula = "", nfold = 3L) {

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

    m_gbm <- gbm::gbm(as.formula(gbm_formula),
      data = train_dat, n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01,
      distribution = distribution)
    #saveRDS(m_gbm, file = paste0(cache_folder, id, "-", .n, "-gbm.rds"))

    test_dat$gbm_ensemble <- tryCatch({gbm::predict.gbm(m_gbm,
      n.trees = m_gbm$n.trees, newdata = test_dat, type = "response")},
      error = function(e) rep(NA, nrow(test_dat)))

    m_rf <- randomForest::randomForest(as.formula(gbm_formula), data = train_dat, ntree = 1000L)
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
