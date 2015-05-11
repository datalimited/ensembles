library("dplyr")
# library("datalimited")
load_all("../datalimited")

ramts <- readRDS("generated-data/ramts.rds")

# unknowns <- filter(ramts, res == "")
# unknowns <- unique(unknowns$stockid)

# Re-fit catch-MSY because there was an issue with the stochastic sigmaR:
message("Running CMSY fits to RAM database...")
ram_cmsy <- ramts %>%
  filter(stockid %in% unique(ramts$stockid)) %>%
# filter(stockid %in% unknowns) %>%
  group_by(stockid) %>%
  do({
    cmsy_out <- cmsy(
      yr             = .$year,
      ct             = .$catch,
      prior_log_mean = .$log_mean[1],
      prior_log_sd   = .$log_sd[1],
      start_r        = resilience(.$res[1], unknown_equals_medium = TRUE),
      sig_r          = 0.05,
      reps           = 10000L)
    if(!is.null(cmsy_out)) {
      bbmsy <- cmsy_out$biomass[, -1] / cmsy_out$bmsy
      bbmsy_out <- summarize_bbmsy(bbmsy, na.rm = TRUE, log = TRUE)
      data.frame(year = .$year, c_touse = .$catch, bbmsy_out)
    } else {
      data.frame(year = .$year, c_touse = .$catch,
        bbmsy_out = rep(NA, nrow(.)))
    }
  }) %>% as.data.frame
ram_cmsy <- select(ram_cmsy, -bbmsy_out)
# saveRDS(ram_cmsy, file = "~/Desktop/ram_cmsy.rds")

# ram_cmsy_unknown_equal_medium <- ramts %>%
#   #filter(stockid %in% unique(ramts$stockid)) %>%
#   filter(stockid %in% unknowns) %>%
#   group_by(stockid) %>%
#   do({
#     cmsy_out <- cmsy(
#       yr             = .$year,
#       ct             = .$catch,
#       prior_log_mean = .$log_mean[1],
#       prior_log_sd   = .$log_sd[1],
#       start_r        = resilience(.$res[1], unknown_equals_medium = TRUE),
#       sig_r          = 0.05,
#       reps           = 5000L)
#     if(!is.null(cmsy_out)) {
#       bbmsy <- cmsy_out$biomass[, -1] / cmsy_out$bmsy
#       bbmsy_out <- summarize_bbmsy(bbmsy, na.rm = TRUE, log = TRUE)
#       data.frame(year = .$year, c_touse = .$catch, bbmsy_out)
#     } else {
#       data.frame(year = .$year, c_touse = .$catch,
#         bbmsy_out = rep(NA, nrow(.)))
#     }
#   }) %>% as.data.frame
# ram_cmsy_unknown_equal_medium <- select(ram_cmsy_unknown_equal_medium,
#   -bbmsy_out)
# saveRDS(ram_cmsy_unknown_equal_medium,
#   file = "~/Desktop/ram_cmsy_unknown_equal_medium.rds")
#
# #message("Running COM-SIR fits to RAM database...")
# #ram_comsir <- ramts %>%
# #  filter(stockid %in% unique(ramts$stockid)[c(1:6)]) %>%
# #  plyr::ddply("stockid", function(.) {})
# #group_by(stockid) %>%
# #  do({
# `.` <- subset(ramts, stockid == unique(ramts$stockid)[2])
# #    # message(.$stocklong[1])
#    comsir_out <- comsir(
#      yr             = .$year,
#      ct             = .$catch,
#      x              = 0.8,
#      k              = 800,
#      r              = 0.6,
#      a              = 0.8,
#      start_r        = resilience(.$res[1]),
#      norm_k         = FALSE,
#      logk           = TRUE,
#      norm_r         = FALSE,
#      norm_a         = FALSE,
#      norm_x         = FALSE,
#      normal_like    = FALSE,
#      nsim           = 2e6, # was 1e6
#      cv             = 0.4,
#      logistic_model = TRUE,
#      obs            = FALSE,
#      n_posterior    = 5e3) # was 5e3
# #    if(!is.null(comsir_out)) {
#      bbmsy <- reshape2::dcast(comsir_out$quantities, sample_id ~ yr,
#        value.var = "bbmsy")[,-1] # convert long to wide format
# #      # TODO switch log to TRUE:
#      bbmsy_out <- summarize_bbmsy(bbmsy, log = FALSE, na.rm = TRUE)
#      data.frame(year = .$year, catch = .$catch, bbmsy_out)
# #    } else {
# #      data.frame(year = .$year, catch = .$catch, bbmsy_out = rep(NA, nrow(.)))
# #    }
# #  }) %>% as.data.frame
#
# library("ggplot2")
# # ggplot(ram_cmsy, aes(year, c_touse)) + geom_line() +
# #   facet_wrap(~stockid, scales = "free_x")
#
# # ram_cmsy$bbmsy_q25[ram_cmsy$bbmsy_q25 < 0] <- 0
# # ram_cmsy$bbmsy_q75[ram_cmsy$bbmsy_q75 > 2.5] <- 2.5
# # ram_cmsy$bbmsy_q2.5[ram_cmsy$bbmsy_q2.5 < 0] <- 0
# # ram_cmsy$bbmsy_q97.5[ram_cmsy$bbmsy_q97.5 > 2.5] <- 2.5

library("ggplot2")
p <- left_join(ram_cmsy, ramts[,c("year", "stockid", "stocklong", "bbmsy_ram"), ])
saveRDS(p, file = "generated-data/ram-cmsy-2015.rds")

pdf("figs/cmsy-ram.pdf", width = 5, height = 4)
plyr::d_ply(p, "stocklong", function(x) {
  message(x$stocklong[1])
  q <- ggplot(x, aes(year, bbmsy_q50)) + geom_line()  +
    geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
    geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), alpha = 0.1) +
    facet_wrap(~stocklong, scales = "free_x") +
    geom_hline(yintercept = 1, lty = 2) +
    scale_y_continuous(limits = c(0, 3)) +
    geom_line(aes(year, bbmsy_ram), colour = "red") +
    theme_bw() +
    theme(plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  print(q)
})
dev.off()

# ggplot(aes(year, bbmsy_q50)) + geom_line()  +
#   geom_ribbon(aes(ymin = bbmsy_q25, ymax = bbmsy_q75), alpha = 0.2) +
#   facet_wrap(~stocklong, scales = "free_x") +
#   geom_hline(yintercept = 1, lty = 2) +
#   scale_y_continuous(limits = c(0, 2.5))
# ggsave("figs/cmsy-ram.pdf", width = 20, height = 20)
#
#
# ggplot(ram_comsir, aes(year, bbmsy_q50)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), alpha = 0.2) +
#   facet_wrap(~stockid, scales = "free_x") +
#   geom_hline(yintercept = 1, lty = 2)
#
# ggplot(ram_comsir, aes(year, bbmsy_q50)) + geom_point() +
#   facet_wrap(~stockid)

# names(ram_cmsy)[3:length(names(ram_cmsy))] <- paste0(names(ram_cmsy)[3:length(names(ram_cmsy))], "_notequal")
#
# x <- ram_cmsy %>%
#   select(stockid, year, bbmsy_q50_notequal, bbmsy_q2.5_notequal, bbmsy_q97.5_notequal) %>%
#   inner_join(select(ram_cmsy_unknown_equal_medium, stockid, year, bbmsy_q50, bbmsy_q2.5, bbmsy_q97.5))
#
#
# xx <- reshape2::melt(x, id.vars = c("stockid", "year"))
#
# library(ggplot2)
# p <- ggplot(xx, aes(year, value, colour = variable)) + facet_wrap(~stockid) + geom_line() + geom_hline(yintercept = 1, lty = 2) + theme_bw()
# ggsave("~/Desktop/resilience-check.pdf", width = 13, height = 8)
