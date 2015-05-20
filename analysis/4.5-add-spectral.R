# bring in spectral values based on Coilin's work:

library("dplyr")

train_spec_mat <- function(x, freq_vec = 1/c(2, 5, 10, 20)) {
  # using AR as smoother, empirical didn't seem to confer much more benefit
  sp <- spec.ar(x, plot = FALSE)
  # approximate at fixed frequencies - necessary as series of different length
  approx(x = sp$freq, y = sp$spec, xout = freq_vec) %>% as.data.frame
}

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

