# apply the simulation-trained ensemble model to the RAM dataset

set.seed(123)
library("randomForest")
library("dplyr")

source("0-ensemble-functions.R")

# prep the RAM data:
ram <- readRDS("generated-data/ram_fits.rds") %>%
  as.data.frame() %>%
  filter(!is.na(b_bmsy_true))

ram$method <- sub("Costello", "mPRM", ram$method)

ram_sum <- ram %>%
  group_by(stockid, method) %>%
  do(mean_slope_bbmsy(.)) %>%
  as.data.frame()

ram_meta <- ram %>%
  group_by(stockid) %>%
  summarise(scientificname = scientificname[1])
ram_sum <- inner_join(ram_sum, ram_meta)

spec_ram <- ram %>%
  arrange(stockid, tsyear) %>%
  filter(method == "CMSY") %>% # pick one
  group_by(stockid) %>%
  do(train_spec_mat(.$catch, freq_vec = 1/c(5, 20))) %>%
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
d_mean <- reshape2::dcast(ram_sum, stockid + scientificname  ~ method,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, stockid, bbmsy_true_mean)) %>%
  inner_join(spec_ram_wide)
d_slope <- reshape2::dcast(ram_sum, stockid + scientificname ~ method,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, stockid, bbmsy_true_slope)) %>%
  inner_join(spec_ram_wide)

# bring in the simulation formatted data to build the simulation-trained models:
d_mean_sim <- readRDS("generated-data/sim-mean-dat.rds")
m_rf <- randomForest::randomForest(
  log(bbmsy_true_mean) ~ CMSY + COMSIR + mPRM + SSCOM,
    # spec_freq_0.05 + spec_freq_0.2,
  data = d_mean_sim, ntree = 1000L)

m_gbm <- gbm::gbm(
  log(bbmsy_true_mean) ~ CMSY + COMSIR + mPRM + SSCOM,
    # spec_freq_0.05 + spec_freq_0.2,
  data = d_mean_sim, distribution = "gaussian",
  n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01)

m_lm <- lm(
  log(bbmsy_true_mean) ~ (CMSY + COMSIR + mPRM + SSCOM)^2,
      # spec_freq_0.05 + spec_freq_0.2)^2,
  data = d_mean_sim)

pdf("../figs/lm-coefs.pdf", width = 5, height = 6)
par(mfrow = c(1, 1))
par(mar = c(4, 12, 1, 3), cex = 0.9)
par(xpd = FALSE)
arm::coefplot(arm::standardize(m_lm), main = "")
par(xpd = NA)
mtext("Standardized regression coefficient", side = 3, line = 3)
dev.off()

# load the RAM data formatted for mPRM:
ram_prm_dat <- readRDS("generated-data/ram_prm_dat.rds")

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
    train_dat <- dplyr::filter(ram_prm_dat, stockid %in% train_ids)
    test_dat  <- dplyr::filter(ram_prm_dat, stockid %in% test_ids)

    mprm <- datalimited::fit_prm(train_dat,
      eqn = log(bbmsy) ~
        mean_scaled_catch +
        scaled_catch +
        scaled_catch1 +
        scaled_catch2 +
        scaled_catch3 +
        scaled_catch4 +
        species_cat +
        catch_to_rolling_max +
        time_to_max +
        years_back +
        initial_slope - 1)

    # on rare occassions the testing dataset will have new factor levels on the
    # species category
    # if that's the case, we'll bypass this chunk and move on:
    test_dat$b_bmsy_est <- tryCatch({datalimited::predict_prm(test_dat, model = mprm)},
      error = function(e) rep(NA, nrow(test_dat)))

    if (sum(is.na(test_dat$b_bmsy_est)) < nrow(test_dat)) {
      # i.e. if there are all NA values then the mPRM prediction failed, move on...
      # this happens extremely rarely

      test_dat$gbm_ensemble <- tryCatch({gbm::predict.gbm(m_gbm,
        n.trees = m_gbm$n.trees, newdata = test_dat, type = "response")},
        error = function(e) rep(NA, nrow(test_dat)))

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

      # now substitute these in for the mPRM values in the original d_mean:
      # note that due to the inner_join this also accomplishes subsetting
      # the full dataset to the cross-validation chunk level:

       d_test <- test_dat_sum %>%
        rename(mPRM = bbmsy_est_mean) %>%
        select(stockid, mPRM) %>%
        dplyr::inner_join(select(d_mean, -mPRM), by = "stockid")

      # we're in trouble if these don't match:
      stopifnot(identical(length(test_ids), nrow(d_test)))

      # now extrapolate with the ensemble models and return the whole data frame
      d_test$rf_ensemble <- exp(predict(m_rf, newdata = d_test))
      d_test$gbm_ensemble <- exp(predict(m_gbm, newdata = d_test, n.trees = m_gbm$n.trees))
      geo_mean <- TRUE
      individual_models <- c("CMSY", "COMSIR", "mPRM", "SSCOM")
      if (geo_mean) {
        d_test$mean_ensemble <- exp(rowMeans(log(d_test[, individual_models])))
      } else {
        d_test$mean_ensemble <- rowMeans(d_test[, individual_models])
      }
      d_test$lm_ensemble <- exp(predict(m_lm, newdata = d_test))
      d_test$.n = .n # for identification purposes
      d_test
    } else {
      print("skipped iteration", .n)
    }
  })
}

library("doParallel")
registerDoParallel(cores = parallel::detectCores())
.parallel <- ifelse(Sys.info()[["sysname"]] == "Windows", FALSE, TRUE)
qq <- plyr::ldply(seq_len(100L), function(i) cv_ensemble_ram(nfold = 3L, .n = i),
  .parallel = .parallel)

d_mean_long <- qq %>%
  select(-spec_freq_0.05, -spec_freq_0.2) %>%
  reshape2::melt(id.vars = c("stockid", "bbmsy_true_mean"),
    variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean)

saveRDS(d_mean_long, file = "generated-data/ram-ensemble-predicted.rds")
