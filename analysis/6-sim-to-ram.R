# apply the simulation-trained ensemble model to the RAM dataset

library("randomForest")
library("dplyr")

source("4-ensemble-functions.R")

# prep the RAM data:
ram <- readRDS("generated-data/ram_fits.rds")

ram_sum <- ram %>%
  group_by(stockid, method) %>%
  do(mean_slope_bbmsy(.)) %>%
  as.data.frame()

ram_meta <- ram %>%
  group_by(stockid) %>%
  summarise(
    scientificname = scientificname[1],
    habitat = habitat[1],
    max_catch = max(catch),
    total_catch = sum(catch))
ram_sum <- inner_join(ram_sum, ram_meta)

spec_ram <- ram %>%
  arrange(stockid, tsyear) %>%
  filter(method == "CMSY") %>% # pick one
  group_by(stockid) %>%
  do(train_spec_mat(.$catch)) %>%
  rename(spec_freq = x, spec_dens = y) %>%
  as.data.frame()

spec_ram_wide <- spec_ram %>%
  mutate(spec_freq = paste0("spec_freq_", spec_freq)) %>%
  reshape2::dcast(stockid ~ spec_freq,
    value.var = "spec_dens")

# save a data frame of 'true' operating model values to merge in:
trues <- select(ram_sum, stockid, bbmsy_true_mean, bbmsy_true_slope)
trues <- trues[!duplicated(trues), ] # one value per operating model stockid

# switch from long to wide format for modelling:
d_mean <- reshape2::dcast(ram_sum, stockid + scientificname + habitat + max_catch ~ method,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, stockid, bbmsy_true_mean)) %>%
  inner_join(spec_ram_wide)
d_slope <- reshape2::dcast(ram_sum, stockid + scientificname + habitat + max_catch ~ method,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, stockid, bbmsy_true_slope)) %>%
  inner_join(spec_ram_wide)

# bring in the simulation formatted data to build the simulation-trained models:
d_mean_sim <- readRDS("generated-data/sim-mean-dat.rds")
m_rf <- randomForest::randomForest(
  log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM +
    spec_freq_0.05 + spec_freq_0.2,
  data = d_mean_sim, ntree = 1000L)

m_gbm <- gbm::gbm(
  log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM +
    spec_freq_0.05 + spec_freq_0.2,
  data = d_mean_sim, distribution = "gaussian",
  n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01)

m_lm <- lm(
  log(bbmsy_true_mean) ~ (CMSY + COMSIR + Costello + SSCOM)^2,
      # max_catch + spec_freq_0.05 + spec_freq_0.2)^2,
  data = d_mean_sim)

library("mgcv")
m_gam <- mgcv::gam(
  log(bbmsy_true_mean) ~ s(CMSY) + s(COMSIR) + s(Costello) + s(SSCOM),
    # + s(CMSY:COMSIR) + s(CMSY:Costello) + s(CMSY:SSCOM) + s(COMSIR:Costello) +
    # s(COMSIR:SSCOM) + s(Costello:SSCOM) +
      # s(max_catch) + s(spec_freq_0.05) + s(spec_freq_0.2),
  data = d_mean_sim)

# join in life-history data for Costello method
# spp_categories is in the datalimited package as data:
ts_dat <- dplyr::left_join(datalimited::ram_ts, datalimited::spp_categories)
# format before to save time, can be done once before cross-validation:
ts_dat <- plyr::ddply(ts_dat, "stockid", function(x) {
  datalimited::format_prm(year = x$year, catch = x$catch, bbmsy = x$bbmsy_ram,
    species_cat = x$spp_category[1L])
})

# Apply the simulation ensembles to a cross-validated version of the RAM dataset:
cv_ensemble_ram <- function(nfold = 3L, .n = 1L) {
  # To be fair to all methods, we will repeatedly rebuild the mPRM/Costello method
  # in cross-validation splits. Each time, we will:
  # 1. rebuild mPRM on the training split
  # 2. predict the testing split with the just-fitted mPRM model
  # 3. to the 'testing' portion: apply the simulation-trained ensemble model
  # 4. repeat on the other splits
  # 5. repeat the whole process N times

  ids <- unique(d_mean$stockid)
  ids_scrambled <- base::sample(ids)
  chunk_size <- floor(length(ids) / nfold)
  chunk_starts <- seq(1L, chunk_size * nfold, chunk_size)
  # note that we exclude those (a small number) that aren't used in this sorting
  # due to rounding; they'll be captured in the next iteration of CV:
  ids_scrambled_cut <- ids_scrambled[seq_len(chunk_size * nfold)]

  out <- plyr::ldply(seq_along(chunk_starts), function(i) {
    test_ids  <- ids_scrambled_cut[
      seq(chunk_starts[i], chunk_starts[i] + chunk_size)] %>% na.omit()
    train_ids <- ids_scrambled_cut[!ids_scrambled_cut %in% test_ids]
    train_dat <- dplyr::filter(ts_dat, stockid %in% train_ids)
    test_dat  <- dplyr::filter(ts_dat, stockid %in% test_ids)

    mprm <- datalimited::fit_prm(train_dat)
    test_dat$b_bmsy_est <- datalimited::predict_prm(test_dat, model = mprm)

    # not all years available because of lagged catches:
    test_dat <- test_dat %>% dplyr::filter(!is.na(b_bmsy_est))
    test_dat <- test_dat %>% dplyr::filter(!is.na(bbmsy))

    test_dat_sum <- test_dat %>%
      rename(b_bmsy_true = bbmsy) %>%
      group_by(stockid) %>%
      do(mean_slope_bbmsy(.)) %>%
      as.data.frame()

    # check that these match, or we're in trouble:
    aa <- dplyr::inner_join(
      select(test_dat_sum, stockid, bbmsy_true_mean),
      select(d_mean, stockid, bbmsy_true_mean) %>%
        rename(bbmsy_true_mean_original = bbmsy_true_mean), by = "stockid")
    identical(as.numeric(aa$bbmsy_true_mean),
      as.numeric(aa$bbmsy_true_mean_original)) %>%
      stopifnot()

    # now substitute these in for the Costello values in the original d_mean:
    # note that due to the inner_join this also accomplishes subsetting
    # the full dataset to the cross-validation chunk level:
    d_test <- test_dat_sum %>%
      rename(Costello = bbmsy_est_mean) %>%
      select(stockid, Costello) %>%
      dplyr::inner_join(select(d_mean, -Costello), by = "stockid")

    # we're in trouble if these don't match:
    stopifnot(identical(length(test_ids), nrow(d_test)))

    # now extrapolate with the ensemble models and return the whole data frame
    d_test$rf_ensemble <- exp(predict(m_rf, newdata = d_test))
    d_test$gbm_ensemble <- exp(predict(m_gbm, newdata = d_test, n.trees = m_gbm$n.trees))
    geo_mean <- TRUE
    individual_models <- c("CMSY", "COMSIR", "Costello", "SSCOM")
    if (geo_mean) {
      d_test$mean_ensemble <- exp(rowMeans(log(d_test[, individual_models])))
    } else {
      d_test$mean_ensemble <- rowMeans(d_test[, individual_models])
    }
    d_test$lm_ensemble <- exp(predict(m_lm, newdata = d_test))
    d_test$gam_ensemble <- exp(predict(m_gam, newdata = d_test))
    d_test$.n = .n # for identification purposes

    d_test
  })
}

library("doParallel")
registerDoParallel(cores = 4L)
qq <- plyr::ldply(seq_len(8L), function(i) cv_ensemble_ram(nfold = 3L, .n = i),
  .parallel = TRUE)

d_mean_long <- qq %>%
  select(-max_catch, -spec_freq_0.05, -spec_freq_0.2, -habitat) %>%
  reshape2::melt(id.vars = c("stockid", "bbmsy_true_mean"),
    variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean)

saveRDS(d_mean_long, file = "generated-data/ram-ensemble-predicted.rds")
