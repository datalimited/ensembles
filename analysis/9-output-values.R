# This file generates 'magic' numbers for input into the paper. These are,
# for example, ranges of MARE or correlation or AUC values. They get saved to
# `../text/values.rda` and are input through knitr in `../text/ms.Rmd`.

get_performance_stats <- function(dat, digits_fold = 1, digits_raw = 2,
  collapse = "--", mre_digits = 2) {

  dat <- dat %>%
    mutate(ensemble = ifelse(grepl("ensemble", clean_method), TRUE, FALSE)) %>%
    mutate(machine = ifelse(grepl("GBM|RF", clean_method), TRUE, FALSE))

  x <- list()

  mach <- dat %>% filter(ensemble == TRUE, machine == TRUE)
  ensemble <- dat %>% filter(ensemble == TRUE)
  ind <- dat %>% filter(ensemble == FALSE)

  x$mach_vs_ind_mare_fold <- c(
    round((min(ind$mare) - max(mach$mare)) / max(mach$mare), digits_fold),
    round((max(ind$mare) - min(mach$mare)) / min(mach$mare), digits_fold)
  ) * 100

#   x$mach_vs_ind_bias <- c(
#     (ind$mre %>% abs %>% max) / (mach$mre %>% abs %>% min) %>% round(digits_fold),
#     (ind$mre %>% abs %>% min) / (mach$mre %>% abs %>% max) %>% round(digits_fold)
#   ) * 100

  x$ind_corr_range <- round(range(ind$corr), digits_raw) %>% formatC(digits = digits_raw, format = "f")
  x$mach_corr_range <- round(range(mach$corr), digits_raw) %>% formatC(digits = digits_raw, format = "f")

  x$ind_mare_range <- round(range(ind$mare), digits_raw) %>% formatC(digits = digits_raw, format = "f")
  x$mach_mare_range <- round(range(mach$mare), digits_raw) %>% formatC(digits = digits_raw, format = "f")

  x$ensemble_mre_range <- ensemble$mre %>% range %>% round(mre_digits) %>% formatC(digits = mre_digits, format = "f")
  x$mach_mre_range <- mach$mre %>% range %>% round(mre_digits) %>% formatC(digits = mre_digits, format = "f")
  x$ind_mre_range <- ind$mre %>% range %>% round(mre_digits) %>% formatC(digits = mre_digits, format = "f")

  lapply(x, paste, collapse = collapse)
}

d_sim_perf_wide <- readRDS("generated-data/d_sim_perf_wide.rds")
mean_sim <- get_performance_stats(d_sim_perf_wide)

re_ram_sum <- readRDS("generated-data/re_ram_sum.rds")
mean_ram <- get_performance_stats(re_ram_sum)

d_slope_error <- readRDS("generated-data/d_slope_error.rds")
slope_sim <- get_performance_stats(d_slope_error, mre_digits = 3L)

auc_sim_mean <- readRDS("generated-data/auc_sim_mean.rds")
auc_sim <- list()

auc_sim$mach_range <- auc_sim_mean %>% filter(grepl("GBM|RF", clean_method)) %>%
  select(auc) %>% range %>% round(2)

auc_sim$ind_range <- auc_sim_mean %>% filter(grepl("Ensemble", clean_method) == FALSE) %>%
  select(auc) %>% range %>% as.numeric() %>% round(2)

d_ram <- readRDS("generated-data/ram-ensemble-predicted.rds")
ram_stocks_n <- length(unique(d_ram$stockid))

save(mean_sim, mean_ram, slope_sim, auc_sim, ram_stocks_n, file = "../text/values.rda")
