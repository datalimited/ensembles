library("dplyr")
# library("datalimited")
load_all("../datalimited")

ramts <- readRDS("generated-data/ramts.rds")
ramts <- ramts %>% filter(!is.na(c_touse)) # TODO won't this create gaps in years?

message("Running CMSY fits to RAM database...")
ram_cmsy <- ramts %>%
  filter(stockid %in% unique(ramts$stockid)[1:6]) %>%
  group_by(stockid) %>%
  do({
    message(.$stocklong[1])
    cmsy_out <- cmsy(
      yr             = .$year,
      ct             = .$c_touse,
      prior_log_mean = .$log_mean[1],
      prior_log_sd   = .$log_sd[1],
      start_r        = resilience(.$res[1]),
      sig_r          = 0.05,
      reps           = 5e4)
    bbmsy <- cmsy_out$biomass[, -1] / cmsy_out$quantities$bmsy
    bbmsy[is.infinite(bbmsy)] <- NA  # TODO investigate the -Inf values
    # TODO switch log to TRUE:
    bbmsy_out <- summarize_bbmsy(bbmsy, na.rm = TRUE, log = FALSE) # TODO investigate the negative values
    data.frame(year = .$year, c_touse = .$c_touse, bbmsy_out)
  }) %>% as.data.frame

message("Running COM-SIR fits to RAM database...")
ram_comsir <- ramts %>%
  filter(stockid %in% unique(ramts$stockid)[1:2]) %>%
  group_by(stockid) %>%
  do({
    message(.$stocklong[1])
    comsir_out <- comsir(
      yr             = .$year,
      ct             = .$c_touse,
      x              = 0.5,
      k              = 800,
      r              = 0.6,
      a              = 0.8,
      start_r        = resilience(.$res[1]),
      norm_k         = TRUE,
      logk           = TRUE,
      norm_r         = FALSE,
      norm_a         = FALSE,
      norm_x         = FALSE,
      nsim           = 1e6, # was 1e6
      cv             = 0.4,
      logistic_model = TRUE,
      obs            = FALSE,
      n_posterior    = 5e3) # was 5e3
    bbmsy <- reshape2::dcast(comsir_out$quantities, yr ~ sample_id,
      value.var = "bbmsy")[,-1]
    # TODO switch log to TRUE:
    bbmsy_out <- summarize_bbmsy(bbmsy, log = FALSE)
  })
    data.frame(year = .$year, c_touse = .$c_touse, bbmsy_out)

library("ggplot2")
# ggplot(ram_cmsy, aes(year, c_touse)) + geom_point() + facet_wrap(~stockid)
ggplot(ram_cmsy, aes(year, bbmsy_q50)) + geom_point() + facet_wrap(~stockid)
