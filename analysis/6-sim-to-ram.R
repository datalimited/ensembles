# apply the simulation-trained ensemble model to the RAM dataset

library("randomForest")
library("dplyr")

source("4-ensemble-functions.R")

# prep the RAM data:
ram <- readRDS("generated-data/ram_fits.rds")

ram_sum <- ram %>%
  group_by(stockid, method) %>%
  do(mean_slope_bbmsy(.)) %>%
  as.data.frame()

ram_meta <- ram %>%
  group_by(stockid) %>%
  summarise(
    habitat = habitat[1],
    max_catch = max(catch),
    total_catch = sum(catch))
ram_sum <- inner_join(ram_sum, ram_meta)

# join in spectral values:
train_spec_mat <- function(x, freq_vec = 1/c(5, 20)) {
  # using AR as smoother, empirical didn't seem to confer much more benefit
  sp <- spec.ar(x, plot = FALSE)
  # approximate at fixed frequencies - necessary as series of different length
  approx(x = sp$freq, y = sp$spec, xout = freq_vec) %>% as.data.frame
}
spec_ram <- ram %>%
  arrange(stockid, tsyear) %>%
  group_by(stockid) %>%
  do(train_spec_mat(.$catch)) %>%
  rename(spec_freq = x, spec_dens = y) %>%
  as.data.frame()

spec_ram_wide <- spec_ram %>%
  mutate(spec_freq = paste0("spec_freq_", spec_freq)) %>%
  reshape2::dcast(stockid ~ spec_freq,
    value.var = "spec_dens")

# save a data frame of 'true' operating model values to merge in:
trues <- select(ram_sum, stockid, bbmsy_true_mean, bbmsy_true_slope)
trues <- trues[!duplicated(trues), ] # one value per operating model stockid

# switch from long to wide format for modelling:
d_mean <- reshape2::dcast(ram_sum, stockid + habitat + max_catch ~ method,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, stockid, bbmsy_true_mean)) %>%
  inner_join(spec_ram_wide)
d_slope <- reshape2::dcast(ram_sum, stockid + habitat + max_catch ~ method,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, stockid, bbmsy_true_slope)) %>%
  inner_join(spec_ram_wide)

# build the simulation-trained model:
d_mean_sim <- readRDS("generated-data/sim-mean-dat.rds")

m_rf <- randomForest::randomForest(
  log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM +
  max_catch + spec_freq_0.05 + spec_freq_0.2, data = d_mean_sim)

m_gbm <- gbm::gbm(
  log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM +
  max_catch + spec_freq_0.05 + spec_freq_0.2, data = d_mean_sim,
  n.trees = 10000L, interaction.depth = 2, shrinkage = 0.001)

m_lm <- lm(
  log(bbmsy_true_mean) ~ (CMSY + COMSIR + Costello + SSCOM +
  max_catch + spec_freq_0.05 + spec_freq_0.2)^2, data = d_mean_sim)

d_mean$rf_ensemble <- predict(m_rf, newdata = d_mean) %>% exp()
d_mean$gbm_ensemble <- predict(m_gbm, newdata = d_mean, n.trees = m_gbm$n.trees) %>% exp()
d_mean$mean_ensemble <- exp(log(d_mean$COMSIR) + log(d_mean$CMSY) + log(d_mean$Costello) + log(d_mean$SSCOM))
d_mean$lm_ensemble <- predict(m_lm, newdata = d_mean) %>% exp()

d_mean_long <- d_mean %>%
  select(-max_catch, -spec_freq_0.05, -spec_freq_0.2, -habitat) %>%
  reshape2::melt(id.vars = c("stockid", "bbmsy_true_mean"),
  variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean)

saveRDS(d_mean_long, file = "generated-data/ram-ensemble-predicted.rds")
