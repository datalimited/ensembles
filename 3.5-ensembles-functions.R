# General purpose functions to be used throughout

# A general function for cross-validation testing ensemble models:
cross_val_ensembles <- function(.n, dat, fraction_train = 0.5,
  gbm_formula = "log(bbmsy_true_mean) ~ CMSY + COM.SIR + Costello + SSCOM + LH",
  geo_mean = TRUE,
  individual_models = c("CMSY", "COM.SIR", "Costello", "SSCOM"),
  id = "mean", cache_folder = "generated-data/cv-sim/", distribution = "gaussian",
  lm_formula = "", glm_formula = "") {

  cv_ids_set <- FALSE # legacy code, leaving in in case needed for RAM stocks
  while(!cv_ids_set) {
    nstocks <- length(unique(dat$stock_id))
    train_ids <- sample(nstocks, round(nstocks * fraction_train))
    test_ids <- seq_len(nstocks)[-train_ids]

    train_stock_ids <- unique(dat$stock_id)[train_ids]
    test_stock_ids <- unique(dat$stock_id)[test_ids]

    train_dat <- filter(dat, stock_id %in% train_stock_ids)
    test_dat <- filter(dat, stock_id %in% test_stock_ids)

    cv_ids_set <- TRUE
  }

  m_gbm <- gbm::gbm(as.formula(gbm_formula),
    data = train_dat, n.trees = 3000L, interaction.depth = 3, shrinkage = 0.001,
    distribution = distribution)
  saveRDS(m_gbm, file = paste0(cache_folder, id, "-", .n, "-gbm.rds"))
  test_dat$gbm_ensemble <- tryCatch({gbm::predict.gbm(m_gbm,
    n.trees = m_gbm$n.trees, newdata = test_dat, type = "response")},
    error = function(e) rep(NA, nrow(test_dat)))

  if (lm_formula != "") {
    m_lm <- lm(as.formula(lm_formula), data = train_dat)
    test_dat$lm_ensemble <- tryCatch({predict(m_lm, newdata = test_dat)},
      error = function(e) rep(NA, nrow(test_dat)))
  }

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
}
