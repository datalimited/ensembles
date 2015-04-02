library("dplyr")
# library("datalimited")
load_all("../datalimited")

ramts <- readRDS("generated-data/ramts.rds")
ramts <- ramts %>% filter(!is.na(c_touse)) # TODO won't this create gaps in years?

message("Running CMSY fits to RAM database...")
ram_cmsy <- ramts %>%
  filter(stockid %in% unique(ramts$stockid)[c(1L,4L)]) %>%
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
      reps           = 1e5)
    bbmsy <- cmsy_out$biomass[, -1] / cmsy_out$quantities$bmsy
    bbmsy[is.infinite(bbmsy)] <- NA  # TODO investigate the -Inf values
    # TODO switch log to TRUE:
    bbmsy_out <- summarize_bbmsy(bbmsy, na.rm = TRUE, log = FALSE) # TODO investigate the negative values
    data.frame(year = .$year, c_touse = .$c_touse, bbmsy_out)
  }) %>% as.data.frame

message("Running COM-SIR fits to RAM database...")
ram_comsir <- ramts %>%
  filter(stockid %in% unique(ramts$stockid)[c(1:6)]) %>%
  plyr::ddply("stockid", function(.) {})
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
      # start_r        = c(0.2, 1.0),
      norm_k         = FALSE,
      logk           = TRUE,
      norm_r         = FALSE,
      norm_a         = FALSE,
      norm_x         = FALSE,
      normal_like    = FALSE,
      nsim           = 1e6, # was 1e6
      cv             = 0.4,
      logistic_model = TRUE,
      obs            = FALSE,
      n_posterior    = 5e3) # was 5e3
    if(!is.null(comsir_out)) {
      bbmsy <- reshape2::dcast(comsir_out$quantities, sample_id ~ yr,
        value.var = "bbmsy")[,-1] # convert long to wide format
      # TODO switch log to TRUE:
      bbmsy_out <- summarize_bbmsy(bbmsy, log = FALSE, na.rm = TRUE)
      data.frame(year = .$year, c_touse = .$c_touse, bbmsy_out)
    } else {
      data.frame(year = .$year, c_touse = .$c_touse, bbmsy_out = rep(NA, nrow(.)))
    }
  }) %>% as.data.frame


library("ggplot2")
ggplot(ram_cmsy, aes(year, c_touse)) + geom_line() + facet_wrap(~stockid, scales = "free_x")

ram_cmsy$bbmsy_q25[ram_cmsy$bbmsy_q25 < 0] <- 0
ram_cmsy$bbmsy_q75[ram_cmsy$bbmsy_q75 > 2.5] <- 2.5
ram_cmsy$bbmsy_q2.5[ram_cmsy$bbmsy_q2.5 < 0] <- 0
ram_cmsy$bbmsy_q97.5[ram_cmsy$bbmsy_q97.5 > 2.5] <- 2.5

left_join(ram_cmsy, ramts[,c("stockid", "stocklong"), ]) %>%
ggplot(aes(year, bbmsy_q50)) + geom_line()  +
  geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
  geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), alpha = 0.1) +
  facet_wrap(~stocklong, scales = "free_x") + geom_hline(yintercept = 1, lty = 2) +
  scale_y_continuous(limits = c(0, 2.5))

ggplot(ram_comsir, aes(year, bbmsy_q50)) + geom_line() + geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), alpha = 0.2) + facet_wrap(~stockid, scales = "free_x") + geom_hline(yintercept = 1, lty = 2)

ggplot(ram_comsir, aes(year, bbmsy_q50)) + geom_point() + facet_wrap(~stockid)

