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
    data.frame(year = .$yr, c_touse = .$c_touse, bbmsy_out)
  }) %>% as.data.frame
