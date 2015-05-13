library("ggplot2")
library("doParallel")
library("datalimited") # install first, e.g.: devtools::install("path/to/datalimited")
registerDoParallel(cores = 4)

ramts <- readRDS("generated-data/ramts.rds")
sig_r_vector <- c(0.05) # so we can iterate across multiple values

# Re-fit catch-MSY because there was an issue with the stochastic sigmaR:
for (i in seq_along(sig_r_vector)) {
  message(paste("Fitting CMSY with sigR of", sig_r_vector[i]))
  ram_cmsy <- plyr::ddply(ramts, "stockid", .parallel = TRUE, .fun = function(.) {
    cmsy_out <- cmsy(yr      = .$year,
                     ct      = .$catch,
                     start_r = resilience(.$res[1]),
                     sig_r   = sig_r_vector[i],
                     reps    = 10000L)
    if (!is.null(cmsy_out)) {
      bbmsy <- cmsy_out$biomass[, -1] / cmsy_out$bmsy
      bbmsy_out <- summarize_bbmsy(bbmsy, na.rm = TRUE, log = TRUE)
      data.frame(year = .$year, c_touse = .$catch, bbmsy_out)
    } else {
      data.frame(year = .$year, c_touse = .$catch, bbmsy_out = rep(NA, nrow(.)))
    }
  })
  ram_cmsy <- select(ram_cmsy, -bbmsy_out)

  p <- dplyr::left_join(ram_cmsy, ramts[,c("year", "stockid", "stocklong", "bbmsy_ram"), ])
  saveRDS(p, file = paste0("generated-data/ram-cmsy-sigr-", sig_r_vector[i], ".rds"))

  pdf(paste0("figs/cmsy-ram-sigr-", sig_r_vector[i], ".pdf"), width = 5, height = 4)
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
}

# Fit to the simulation time series:
# TODO: need to get the data from someone
