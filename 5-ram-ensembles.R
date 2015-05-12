library("dplyr")
library("ggplot2")
source("3.5-ensembles-functions.R")

if (!file.exists("generated-data/cv-ram/")) system("mkdir generated-data/cv-ram")

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

# add a binary above/below B/Bmsy = 1 column for a potential response:
ram_sum <- ram_sum %>% mutate(
  above_bbmsy1_true = ifelse(bbmsy_true_mean > 1, 1, 0))

# save a data frame of 'true' operating model values to merge in:
trues <- select(ram_sum, stockid, above_bbmsy1_true,
  bbmsy_true_mean, bbmsy_true_slope)
trues <- trues[!duplicated(trues), ] # one value per operating model stockid

# switch from long to wide format for modelling:
d_mean <- reshape2::dcast(ram_sum, stockid + habitat ~ method,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, -bbmsy_true_slope))
d_slope <- reshape2::dcast(ram_sum, stockid + habitat ~ method,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, -bbmsy_true_mean))

# run a model on all the data to generate data for partial dependence plots:
m <- gbm::gbm(log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM + habitat,
  data = d_mean, n.trees = 3000L, interaction.depth = 3, shrinkage = 0.001)
partial <- plyr::ldply(seq_len(5), function(i) {
  dd <- gbm::plot.gbm(m, i.var = i, return.grid = TRUE)
  dd$predictor <- names(dd)[1]
  names(dd)[1] <- "predictor_value"
  dd$y <- exp(dd$y)
  dd
})

ggplot(partial, aes(predictor_value, y)) + geom_line() + facet_wrap(~predictor) + xlim(0, 3)

# work through cross validation of ensemble models:
library("doParallel")
registerDoParallel(cores = 4)

cv_ram_mean <- plyr::ldply(seq_len(8), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE, id = "sim-mean",
      gbm_formula = "log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM + habitat",
      lm_formula = "log(bbmsy_true_mean) ~ (CMSY + COMSIR + Costello + SSCOM + habitat)^2"))
cv_sim_mean$gbm_ensemble <- exp(cv_sim_mean$gbm_ensemble)
cv_sim_mean$lm_ensemble <- exp(cv_sim_mean$lm_ensemble)

cv_ram_slope <- plyr::ldply(seq_len(8), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_slope, geo_mean = FALSE, id = "sim-slope",
     gbm_formula = "bbmsy_true_slope ~ CMSY + COMSIR + Costello + SSCOM + habitat",
     lm_formula = "bbmsy_true_slope ~ (CMSY + COMSIR + Costello + SSCOM + habitat)^2"))

cv_ram_binary <- plyr::ldply(seq_len(4), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE,
      id = "sim-mean", distribution = "bernoulli",
      gbm_formula = "above_bbmsy1_true ~ CMSY + COMSIR + Costello + SSCOM + habitat",
      glm_formula = "above_bbmsy1_true ~ (CMSY + COMSIR + Costello + SSCOM + habitat)^2"))
saveRDS(cv_ram_binary, file = "generated-data/cv_ram_binary.rds") # used in 5-roc.R

# -------------------------------------------------------
# now switch to long format data, summarize, and compare:
cv_ram_mean_long <- reshape2::melt(select(cv_ram_mean, -above_bbmsy1_true),
  id.vars = c("stockid", "test_iter", "habitat", "bbmsy_true_mean"),
  variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean) %>%
  mutate(
    type = "mean",
    bbmsy_true_trans = log(bbmsy_true),
    bbmsy_est_trans = log(bbmsy_est))

cv_ram_slope_long <- reshape2::melt(select(cv_ram_slope, -above_bbmsy1_true),
  id.vars = c("stockid",  "test_iter", "habitat", "bbmsy_true_slope"),
  variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_slope) %>%
  mutate(
    type = "slope",
    bbmsy_true_trans = bbmsy_true,
    bbmsy_est_trans = bbmsy_est)

ggplot(cv_ram_mean_long, aes(bbmsy_true, bbmsy_est)) +
  geom_point(alpha = 0.07) +
    facet_wrap(~method) + ylim(0, 3) + xlim(0, 3)

ggplot(cv_ram_slope_long, aes(bbmsy_true, bbmsy_est)) +
  geom_point(alpha = 0.07) +
    facet_wrap(~method) + xlim(-.5, .5) + ylim(-.5, .5)

cv_ram_long <- suppressWarnings(
  dplyr::bind_rows(cv_ram_mean_long, cv_ram_slope_long))
saveRDS(cv_ram_long, "generated-data/cv_ram_long.rds")

cors <- cv_ram_long %>% group_by(method, type) %>%
  summarise(spearman = cor(bbmsy_true_trans, bbmsy_est_trans, method = "spearman",
    use = "pairwise.complete.obs")) %>% as.data.frame %>%
  mutate(spearman = round(spearman, 4)) %>%
  arrange(type, -spearman)

ggplot(cors, aes(x = method, xend = method, yend = spearman)) +
  geom_segment(y = 0, lwd = 1.2) + facet_wrap(~type)

cors %>% tidyr::spread(type, spearman) %>%
  ggplot(aes(mean, slope)) + geom_point() +
  geom_text(aes(label = method), hjust = 0.5, vjust = -0.5) +
  xlab("Spearman correlation of B/BMSY mean") +
  ylab("Spearman correlation of B/BMSY slope")

re <- cv_ram_long %>% mutate(
  sq_er = (bbmsy_est_trans - bbmsy_true_trans)^2,
  re    = (bbmsy_est_trans - bbmsy_true_trans) / bbmsy_true_trans,
  pe    = bbmsy_est_trans / bbmsy_true_trans)

# make a data frame that summarizes all these performance metrics at once:
re2 <- re %>% group_by(type, method) %>%
  summarize(
    performance = median(abs(re), na.rm = TRUE),
    performance_l = quantile(abs(re), 0.25, na.rm = TRUE),
    performance_u = quantile(abs(re), 0.75, na.rm = TRUE)) %>%
  mutate(summary = "MARE")

cors2 <- cors %>% mutate(
  summary = "spearman",
  performance = spearman,
  performance_l = spearman,
  performance_u = spearman) %>%
  select(-spearman)

performance <- bind_rows(re2, cors2) %>% as.data.frame

performance %>% reshape2::dcast(method + type ~ summary, value.var = c("performance")) %>%
  ggplot(aes(MARE, spearman)) + geom_point() + geom_text(aes(label = method), hjust = 0.2) +
  facet_wrap(~type, scales = "free") +
  xlab("MARE within stocks") + ylab("Spearman correlation across stocks")
