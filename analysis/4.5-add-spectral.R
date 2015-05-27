# bring in spectral values based on Coilin's work:

library("dplyr")

source("4-ensemble-functions.R")
# raw data from:
# ~/Dropbox/FisheriesWorkingGroupPhaseII/Decision_trees/harvest_dynamics_classification/StochasticSimFullFactorial_2013-06-21_no_status.RData

load("raw-data/StochasticSimFullFactorial_2013-06-21_no_status.RData")

spec <- simsStoch %>%
  arrange(stock_id, year) %>%
  group_by(stock_id, sigmaC, sigmaR, LH, iter, ED) %>%
  do(train_spec_mat(.$catch)) %>%
  rename(spec_freq = x, spec_dens = y) %>%
  as.data.frame()

library("ggplot2")
p <- ggplot(spec, aes(spec_freq, log(spec_dens), group = spec_freq)) +
  geom_boxplot() + facet_wrap(~ED)
ggsave("../figs/spectral-distributions.pdf", width = 8, height = 7)

spec_wide <- spec %>%
  mutate(spec_freq = paste0("spec_freq_", spec_freq)) %>%
  reshape2::dcast(stock_id + sigmaC + sigmaR + LH + iter + ED ~ spec_freq,
    value.var = "spec_dens")

saveRDS(spec_wide, file = "generated-data/spec-frequencies-sim.rds")

