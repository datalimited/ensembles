library(dplyr)
library(ggplot2)
if (!file.exists("generated-data/cv-gbms/")) system("mkdir generated-data/cv-gbms")

dsim <- readRDS("raw-data/batch1-results.rds")

dsim <- dsim %>%
  transform(stock_id =
  paste0(stock_id, "_sigmaC_", sigmaC, "_sigmaR_", sigmaR, "_lh", LH, "_it_", iter)) %>%
  arrange(stock_id, iter, year)

if(!file.exists("generated-data/dsim_sum.rds")) {
  dsim_sum <- dsim %>%
    group_by(stock_id, method_id, LH) %>%
    do({

      years_mean <- 3
      .n <- nrow(.)
      i <- seq(.n-(years_mean-1), .n)

      bbmsy_true_mean = mean(.$b_bmsy_true[i])
      bbmsy_est_mean = mean(.$b_bmsy_est[i])

      ytrue <- .$b_bmsy_true[i]
      yest <- .$b_bmsy_est[i]
      bbmsy_true_slope <- coef(mblm::mblm(ytrue ~ i))[[2]]
      bbmsy_est_slope <- coef(mblm::mblm(yest ~ i))[[2]]

      data.frame(bbmsy_true_mean, bbmsy_est_mean, bbmsy_true_slope, bbmsy_est_slope)
    }) %>% as.data.frame
    saveRDS(dsim_sum, file = "generated-data/dsim_sum.rds")
} else {
  dsim_sum <- readRDS("generated-data/dsim_sum.rds")
}

dsim_meta <- dsim %>%
  group_by(stock_id) %>%
  summarise(
    max_catch = max(catch),
    total_catch = sum(catch))

dsim_sum <- inner_join(dsim_sum, dsim_meta)

trues <- select(dsim_sum, stock_id, method_id, bbmsy_true_mean, bbmsy_true_slope)
trues <- trues[!duplicated(trues), ]

d_mean <- reshape2::dcast(dsim_sum, stock_id + LH ~ method_id,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, -bbmsy_true_slope))

d_slope <- reshape2::dcast(dsim_sum, stock_id + LH ~ method_id,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, -bbmsy_true_mean))

m <- gbm::gbm(log(bbmsy_true_mean) ~ CMSY + COM.SIR + Costello + SSCOM + LH,
#m <- gbm::gbm(bbmsy_true_mean ~ CMSY + COM.SIR + Costello + SSCOM + LH,
  data = d_mean, n.trees = 3000L, interaction.depth = 3, shrinkage = 0.001)
partial <- plyr::ldply(seq_len(5), function(i) {
  dd <- gbm::plot.gbm(m, i.var = i, return.grid = TRUE)
  dd$predictor <- names(dd)[1]
  names(dd)[1] <- "predictor_value"
  dd$y <- exp(dd$y)
  dd
})

ggplot(partial, aes(predictor_value, y)) + geom_line() + facet_wrap(~predictor) + xlim(0, 3)

# a general function for cross-validation testing ensemble models:

cross_val_ensembles <- function(.n, dat, fraction_train = 0.5,
  formula = "log(bbmsy_true_mean) ~ CMSY + COM.SIR + Costello + SSCOM + LH",
  id = "mean", geomean = TRUE) {

  cv_ids_set <- FALSE # legacy code, leaving in in case needed for RAM stocks
  while(!cv_ids_set) {
    nstocks <- length(unique(dat$stock_id))

    train_ids <- sample(nstocks, round(nstocks * fraction_train))
    test_ids <- seq_len(nstocks)[-train_ids]

    train_stock_ids <- unique(dat$stock_id)[train_ids]
    test_stock_ids <- unique(dat$stock_id)[test_ids]

    train_dat <- filter(dat, stock_id %in% train_stock_ids)
    test_dat <- filter(dat, stock_id %in% test_stock_ids)

    # make sure we have all 'test' species categories in our 'train' data
    # or the lm() prediction will fail
    cv_ids_set <- TRUE
  }

  m_gbm <- gbm::gbm(as.formula(formula),
    data = train_dat, n.trees = 2000L, interaction.depth = 3, shrinkage = 0.001,
    distribution = "gaussian")

  saveRDS(m_gbm, file = paste0("generated-data/cv-gbms/", id, "-", .n, "-gbm.rds"))

  test_dat$gbm_ensemble <- tryCatch({gbm::predict.gbm(m_gbm,
    n.trees = m_gbm$n.trees, newdata = test_dat)},
    error = function(e) rep(NA, length(test_ids)))

  if (geomean) {
  test_dat$mean_ensemble <- exp(c(log(test_dat$CMSY), log(test_dat$COM.SIR),
    log(test_dat$Costello), log(test_dat$SSCOM)))
  } else {
  test_dat$mean_ensemble <- c(test_dat$CMSY, test_dat$COM.SIR,
    test_dat$Costello, test_dat$SSCOM)
  }

  test_dat$test_iter <- .n
  test_dat
}

library("doParallel")
registerDoParallel(cores = 4)

cv_sim_mean <- plyr::ldply(seq_len(8), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, id = "sim-mean",
      formula = "log(bbmsy_true_mean) ~ CMSY + COM.SIR + Costello + SSCOM + LH"))
cv_sim_mean$gbm_ensemble <- exp(cv_sim_mean$gbm_ensemble)

cv_sim_slope <- plyr::ldply(seq_len(8), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles( .n = .n, dat = d_slope, id = "sim-slope",
      geomean = FALSE,
      formula = "bbmsy_true_slope ~ CMSY + COM.SIR + Costello + SSCOM + LH"))

############## fix from here:
# now melt all the methods:
cv_sim_mean_long <- reshape2::melt(cv_sim_mean,
  id.vars = c("stock_id", "test_iter", "LH", "bbmsy_true_mean"),
  variable.name = "method_id", value.name = "bbmsy_est_mean")

cv_sim_slope_long <- reshape2::melt(cv_sim_slope,
  id.vars = c("stock_id", "test_iter", "LH", "bbmsy_true_slope"),
  variable.name = "method_id", value.name = "bbmsy_est_slope")

ggplot(cv_sim_mean_long, aes(bbmsy_true_mean, bbmsy_est_mean)) +
  geom_point(alpha = 0.03) +
    facet_wrap(~method_id) + ylim(0, 3) + xlim(0, 3)

ggplot(cv_sim_slope_long, aes(bbmsy_true_slope, bbmsy_est_slope)) +
  geom_point(alpha = 0.03) +
    facet_wrap(~method_id) + xlim(-.5, .5) + ylim(-.5, .5)

cv_sim_mean_long %>% group_by(method_id) %>%
  summarise(spearman = cor(bbmsy_true_mean, bbmsy_est_mean, method = "spearman",
    use = "pairwise.complete.obs"))

cv_sim_slope_long %>% group_by(method_id) %>%
  summarise(spearman = cor(bbmsy_true_slope, bbmsy_est_slope, method = "spearman",
    use = "pairwise.complete.obs"))

re <- cv_sim_mean_long %>% mutate(
   sq_er = (log(bbmsy_est_mean) - log(bbmsy_true_mean))^2,
  #sq_er = ((bbmsy_est_mean) - (bbmsy_true_mean))^2,
   re = (log(bbmsy_est_mean) - log(bbmsy_true_mean)) / log(bbmsy_true_mean)) %>%
  #re = ((bbmsy_est_mean) - (bbmsy_true_mean)) / (bbmsy_true_mean)) %>%
  group_by(method_id, test_iter) %>%
  summarize(
    mare = median(abs(re), na.rm = TRUE),
    rmse = sqrt(mean(sq_er, na.rm = TRUE)))

ggplot(re, aes(method_id, mare)) + geom_boxplot()
ggplot(re, aes(method_id, rmse)) + geom_boxplot()
