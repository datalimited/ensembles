# Save as 3-weights file, but use the simulation output to build the ensemble model

library("dplyr")

dsim <- readRDS("raw-data/batch1-results.rds")
# note that there were duplicate rows here:
# dsim <- dsim[!duplicated(dsim), ] # takes a long time
# saveRDS(dsim, file = "raw-data/batch1-results.rds")

# combine over iterations to make the model fitting tractable:
dsim2 <- dsim %>%
  group_by(stock_id, year, method_id, sigmaC, sigmaR, ED, SEL, TS, AR, UR) %>%
  summarise(
    b_bmsy_true = median(b_bmsy_true, na.rm = TRUE),
    b_bmsy_est = median(b_bmsy_est, na.rm = TRUE)) %>%
  arrange(method_id, stock_id, year) %>%
  as.data.frame

# now drop all but last 10 years:
dsim3 <- plyr::ddply(dsim2,
  c("stock_id", "method_id", "sigmaC", "sigmaR", "ED", "SEL", "TS", "AR", "UR"),
  function(x) {
    if (nrow(x) >= 10) { # note that some are shorter - ask why
      out <- x[(nrow(x)-9):nrow(x),]
      out$year_before_end <- seq_len(10)
      out
    }
  })

dsim_wide <- reshape2::dcast(dsim3,
  stock_id + year + year_before_end + sigmaC + AR + UR + ED + SEL + TS + sigmaR +
    b_bmsy_true ~ method_id,
  value.var = "b_bmsy_est")
# lots of NAs, not sure why:
dsim_wide <- na.omit(dsim_wide)

m_gbm_sim <- gbm::gbm(b_bmsy_true ~
    COM.SIR + Costello + SSCOM + CMSY,
  data = dsim_wide, distribution = "gaussian", n.trees = 2500, interaction.depth = 10,
  shrinkage = 0.05)

pdf("figs/gbm-partial-sim.pdf", width = 8, height = 8)
par(cex = 0.6)
par(mfrow = c(2, 2));for(i in 1:4) gbm::plot.gbm(m_gbm_sim, i.var = i, ylim = c(0.7, 1.5))
dev.off()

saveRDS(m_gbm_sim, file = "generated-data/m_gbm_sim.rds")
