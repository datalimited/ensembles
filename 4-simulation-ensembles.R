library("dplyr")
library("ggplot2")
source("3.5-ensembles-functions.R")

if (!file.exists("generated-data/cv-sim/")) system("mkdir generated-data/cv-sim")

dsim <- readRDS("raw-data/batch1-results.rds")

# extend the stock_id to be truly unique per scenario:
# (note that there are still multiple iterations with the same stock_id)
dsim <- dsim %>%
  transform(stock_id =
  paste0(stock_id, "_sigmaC_", sigmaC, "_sigmaR_", sigmaR, "_lh", LH)) %>%
  arrange(stock_id, iter, year) # critical since not all in order

# summarise the mean and slope in the last N years:
if(!file.exists("generated-data/dsim_sum.rds")) {
  dsim_sum <- dsim %>%
    group_by(stock_id, method_id, iter) %>%
    do({
      years_mean <- 3
      .n <- nrow(.)
      i <- seq(.n-(years_mean-1), .n)
      bbmsy_true_mean = mean(.$b_bmsy_true[i])
      bbmsy_est_mean = mean(.$b_bmsy_est[i])
      ytrue <- .$b_bmsy_true[i]
      yest <- .$b_bmsy_est[i]
      bbmsy_true_slope <- coef(mblm::mblm(ytrue ~ i))[[2]]
      bbmsy_est_slope <- coef(mblm::mblm(yest ~ i))[[2]]
      data.frame(bbmsy_true_mean, bbmsy_est_mean, bbmsy_true_slope, bbmsy_est_slope)
    }) %>% as.data.frame
    saveRDS(dsim_sum, file = "generated-data/dsim_sum.rds")
} else {
  dsim_sum <- readRDS("generated-data/dsim_sum.rds")
}

# join in some characteristics that we'll use in models:
dsim_meta <- dsim %>%
  group_by(stock_id, iter) %>%
  summarise(
    LH = LH[1],
    max_catch = max(catch),
    total_catch = sum(catch))
dsim_sum <- inner_join(dsim_sum, dsim_meta)

# add a binary above/below B/Bmsy = 1 column for a potential response:
dsim_sum <- dsim_sum %>% mutate(
  above_bbmsy1_true = ifelse(bbmsy_true_mean > 1, 1, 0))

# save a data frame of 'true' operating model values to merge in:
trues <- select(dsim_sum, stock_id, iter, above_bbmsy1_true,
  bbmsy_true_mean, bbmsy_true_slope)
trues <- trues[!duplicated(trues), ] # one value per operating model stock_id

# switch from long to wide format for modelling:
d_mean <- reshape2::dcast(dsim_sum, stock_id + iter + LH ~ method_id,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, -bbmsy_true_slope))
d_slope <- reshape2::dcast(dsim_sum, stock_id + iter + LH ~ method_id,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, -bbmsy_true_mean))

# run a model on all the data to generate data for partial dependence plots:
m <- gbm::gbm(log(bbmsy_true_mean) ~ CMSY + COM.SIR + Costello + SSCOM + LH,
  data = d_mean, n.trees = 3000L, interaction.depth = 3, shrinkage = 0.001)
partial <- plyr::ldply(seq_len(5), function(i) {
  dd <- gbm::plot.gbm(m, i.var = i, return.grid = TRUE)
  dd$predictor <- names(dd)[1]
  names(dd)[1] <- "predictor_value"
  dd$y <- exp(dd$y)
  dd
})

# partial dependence plot:
ggplot(partial, aes(predictor_value, y)) + geom_line() + facet_wrap(~predictor) + xlim(0, 3)

# work through cross validation of ensemble models:
library("doParallel")
registerDoParallel(cores = 4)

cv_sim_mean <- plyr::ldply(seq_len(8), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE, id = "sim-mean",
      gbm_formula = "log(bbmsy_true_mean) ~ CMSY + COM.SIR + Costello + SSCOM + LH",
      lm_formula = "log(bbmsy_true_mean) ~ (CMSY + COM.SIR + Costello + SSCOM + LH)^2"))
cv_sim_mean$gbm_ensemble <- exp(cv_sim_mean$gbm_ensemble)
cv_sim_mean$lm_ensemble <- exp(cv_sim_mean$lm_ensemble)

cv_sim_slope <- plyr::ldply(seq_len(8), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_slope, geo_mean = FALSE, id = "sim-slope",
     gbm_formula = "bbmsy_true_slope ~ CMSY + COM.SIR + Costello + SSCOM + LH",
     lm_formula = "bbmsy_true_slope ~ (CMSY + COM.SIR + Costello + SSCOM + LH)^2"))

cv_sim_binary <- plyr::ldply(seq_len(4), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE,
      id = "sim-mean", distribution = "bernoulli",
      gbm_formula = "above_bbmsy1_true ~ CMSY + COM.SIR + Costello + SSCOM + LH",
      glm_formula = "above_bbmsy1_true ~ (CMSY + COM.SIR + Costello + SSCOM + LH)^2"))
saveRDS(cv_sim_binary, file = "generated-data/cv_sim_binary.rds") # used in 5-roc.R

# -------------------------------------------------------
# now switch to long format data, summarize, and compare:
cv_sim_mean_long <- reshape2::melt(select(cv_sim_mean, -above_bbmsy1_true),
  id.vars = c("stock_id", "iter", "test_iter", "LH", "bbmsy_true_mean"),
  variable.name = "method_id", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean) %>%
  mutate(
    type = "mean",
    bbmsy_true_trans = log(bbmsy_true),
    bbmsy_est_trans = log(bbmsy_est))

cv_sim_slope_long <- reshape2::melt(select(cv_sim_slope, -above_bbmsy1_true),
  id.vars = c("stock_id", "iter", "test_iter", "LH", "bbmsy_true_slope"),
  variable.name = "method_id", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_slope) %>%
  mutate(
    type = "slope",
    bbmsy_true_trans = bbmsy_true,
    bbmsy_est_trans = bbmsy_est)

ggplot(cv_sim_mean_long, aes(bbmsy_true, bbmsy_est)) +
  geom_point(alpha = 0.01) +
    facet_wrap(~method_id) + ylim(0, 3) + xlim(0, 3)

ggplot(cv_sim_slope_long, aes(bbmsy_true, bbmsy_est)) +
  geom_point(alpha = 0.01) +
    facet_wrap(~method_id) + xlim(-.5, .5) + ylim(-.5, .5)

cv_sim_long <- suppressWarnings(
  dplyr::bind_rows(cv_sim_mean_long, cv_sim_slope_long))
saveRDS(cv_sim_long, "generated-data/cv_sim_long.rds")

cors <- cv_sim_long %>% group_by(method_id, type) %>%
  summarise(spearman = cor(bbmsy_true_trans, bbmsy_est_trans, method = "spearman",
    use = "pairwise.complete.obs")) %>% as.data.frame %>%
  mutate(spearman = round(spearman, 4)) %>%
  arrange(type, -spearman)

ggplot(cors, aes(x = method_id, xend = method_id, yend = spearman)) +
  geom_segment(y = 0, lwd = 1.2) + facet_wrap(~type)

cors %>% tidyr::spread(type, spearman) %>%
  ggplot(aes(mean, slope)) + geom_point() +
  geom_text(aes(label = method_id), hjust = 0.5, vjust = -0.5) +
  xlab("Spearman correlation of B/BMSY mean") +
  ylab("Spearman correlation of B/BMSY slope")

re <- cv_sim_long %>% mutate(
  sq_er = (bbmsy_est_trans - bbmsy_true_trans)^2,
  re    = (bbmsy_est_trans - bbmsy_true_trans) / bbmsy_true_trans,
  pe    = bbmsy_est_trans / bbmsy_true_trans)

# make a data frame that summarizes all these performance metrics at once:
re2 <- re %>% group_by(type, method_id) %>%
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

performance %>% reshape2::dcast(method_id + type ~ summary, value.var = c("performance")) %>%
  ggplot(aes(MARE, spearman)) + geom_point() + geom_text(aes(label = method_id), hjust = 0.2) +
  facet_wrap(~type, scales = "free") +
  xlab("MARE within stocks") + ylab("Spearman correlation across stocks")
