# The main file to run ensemble models for the simulation data.

library("dplyr")
library("ggplot2")
source("4-ensemble-functions.R")

if (!file.exists("generated-data/cv-sim/")) system("mkdir generated-data/cv-sim")

dsim <- readRDS("raw-data/batch1-results.rds")
dsim$method_id <- sub("COM.SIR", "COMSIR", dsim$method_id) # to match RAM fits
dsim <- rename(dsim, stockid = stock_id, method = method_id) # to match RAM fits

# extend the stockid to be truly unique per scenario:
# (note that there are still multiple iterations with the same stockid)
# TODO: should 'scenario' be across sigmaC and sigmaR values too?
dsim <- dsim %>%
  mutate(stockid =
    paste0(stockid, "_sigmaC_", sigmaC, "_sigmaR_", sigmaR, "_lh", LH)) %>%
  arrange(stockid, iter, year) # critical since not all in order

# summarise the mean and slope in the last N years:
library("doParallel")
registerDoParallel(cores = 4)

dsim_sum <- plyr::ddply(dsim, c("stockid", "method", "iter"),
  .parallel = TRUE, .fun = mean_slope_bbmsy)
saveRDS(dsim_sum, file = "generated-data/dsim_sum.rds")
# dsim_sum <- readRDS("generated-data/dsim_sum.rds")

# join in some characteristics that we'll use in models:
dsim_meta <- dsim %>%
  group_by(stockid, iter) %>%
  summarise(
    LH = LH[1],
    max_catch = max(catch),
    total_catch = sum(catch))
dsim_sum <- inner_join(dsim_sum, dsim_meta)

# add a binary above/below B/Bmsy = 1 column for a potential response:
dsim_sum <- dsim_sum %>% mutate(
  above_bbmsy1_true = ifelse(bbmsy_true_mean > 1, 1, 0))

# save a data frame of 'true' operating model values to merge in:
trues <- select(dsim_sum, stockid, iter, above_bbmsy1_true,
  bbmsy_true_mean, bbmsy_true_slope)
trues <- trues[!duplicated(trues), ] # one value per operating model stockid

# switch from long to wide format for modelling:
d_mean <- reshape2::dcast(dsim_sum, stockid + iter + LH ~ method,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, -bbmsy_true_slope))
d_slope <- reshape2::dcast(dsim_sum, stockid + iter + LH ~ method,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, -bbmsy_true_mean))

# strip out those with NAs (in CMSY)
# create problems with some models, such as randomForest, otherwise
d_mean <- na.omit(d_mean)
d_slope <- na.omit(d_slope)

# run a model on all the data to generate data for partial dependence plots:
m <- gbm::gbm(log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM + LH,
  data = d_mean, n.trees = 10000L, interaction.depth = 2, shrinkage = 0.001)
partial <- plyr::ldply(seq_len(5), function(i) {
  dd <- gbm::plot.gbm(m, i.var = i, return.grid = TRUE)
  dd$predictor <- names(dd)[1]
  names(dd)[1] <- "predictor_value"
  dd$y <- exp(dd$y)
  dd
})

# partial dependence plot:
p <- ggplot(partial, aes(predictor_value, y)) + geom_line() +
  facet_wrap(~predictor) + xlim(0, 3)
ggsave("figs/partial-sim.pdf", width = 7, height = 5)

partial_2d <- plyr::ldply(1:5, function(x) plyr::ldply(1:5, function(y) {
  dd <- gbm::plot.gbm(m, i.var = c(x, y), return.grid = TRUE, continuous.resolution = 20)
  dd$var1 <- names(dd)[1]
  dd$var2 <- names(dd)[2]
  names(dd)[1] <- "x"
  names(dd)[2] <- "y"
  names(dd)[3] <- "z"
  dd$z <- exp(dd$z)
  dd
}))

# check colour pallete:
zlim <- c(min(partial_2d$z), max(partial_2d$z))
pal <- colorRampPalette(c("blue", "white", "red"))(16)[-c(13:16)]
#plot(seq(min(zlim), max(zlim), length.out = 12), 1:12, col = pal)
pdf("figs/partial-sim-2d.pdf", width = 9, height = 9)
par(mfrow = c(5, 5), mar = c(3,3,1,1), oma = c(4, 4, 1, 1), cex = 0.5)
par(xpd = NA, mgp = c(1.5, 0.5, 0))
plyr::d_ply(partial_2d, c("var1", "var2"), function(x) {
  xx <- reshape2::dcast(x, x ~ y, value.var = "z")
  image(xx[,1], as.numeric(colnames(xx)[-1]), as.matrix(xx[,-1]), main = "",
    xlab = unique(x$var1), ylab = unique(x$var2),
    zlim = zlim,
    col = colorRampPalette(c("blue", "white", "red"))(16)[-c(13:16)])
})
dev.off()
message(paste("zlim were", round(zlim, 2), collapse = " "))

# work through cross validation of ensemble models:
cv_sim_mean <- plyr::ldply(seq_len(4), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE, id = "sim-mean",
      gbm_formula = "log(bbmsy_true_mean) ~ CMSY + COMSIR + Costello + SSCOM + LH",
      lm_formula = "log(bbmsy_true_mean) ~ (CMSY + COMSIR + Costello + SSCOM + LH)^2"))
cv_sim_mean$gbm_ensemble <- exp(cv_sim_mean$gbm_ensemble)
cv_sim_mean$lm_ensemble <- exp(cv_sim_mean$lm_ensemble)

cv_sim_slope <- plyr::ldply(seq_len(4), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_slope, geo_mean = FALSE, id = "sim-slope",
     gbm_formula = "bbmsy_true_slope ~ CMSY + COMSIR + Costello + SSCOM + LH",
     lm_formula = "bbmsy_true_slope ~ (CMSY + COMSIR + Costello + SSCOM + LH)^2"))

cv_sim_binary <- plyr::ldply(seq_len(4), .parallel = TRUE,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE,
      id = "sim-mean", distribution = "bernoulli",
      gbm_formula = "above_bbmsy1_true ~ CMSY + COMSIR + Costello + SSCOM + LH",
      glm_formula = "above_bbmsy1_true ~ (CMSY + COMSIR + Costello + SSCOM + LH)^2"))
saveRDS(cv_sim_binary, file = "generated-data/cv_sim_binary.rds") # used in 5-roc.R

# -------------------------------------------------------
# now switch to long format data, summarize, and compare:
cv_sim_mean_long <- reshape2::melt(select(cv_sim_mean, -above_bbmsy1_true),
  id.vars = c("stockid", "iter", "test_iter", "LH", "bbmsy_true_mean"),
  variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean) %>%
  mutate(
    type = "mean",
    bbmsy_true_trans = log(bbmsy_true),
    bbmsy_est_trans = log(bbmsy_est))

cv_sim_slope_long <- reshape2::melt(select(cv_sim_slope, -above_bbmsy1_true),
  id.vars = c("stockid", "iter", "test_iter", "LH", "bbmsy_true_slope"),
  variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_slope) %>%
  mutate(
    type = "slope",
    bbmsy_true_trans = bbmsy_true,
    bbmsy_est_trans = bbmsy_est)

# ggplot(cv_sim_mean_long, aes(bbmsy_true, bbmsy_est)) +
#   geom_point(alpha = 0.01) +
#     facet_wrap(~method) + ylim(0, 3) + xlim(0, 3)
#
# ggplot(cv_sim_slope_long, aes(bbmsy_true, bbmsy_est)) +
#   geom_point(alpha = 0.01) +
#     facet_wrap(~method) + xlim(-.5, .5) + ylim(-.5, .5)

cv_sim_long <- suppressWarnings(
  dplyr::bind_rows(cv_sim_mean_long, cv_sim_slope_long))
saveRDS(cv_sim_long, "generated-data/cv_sim_long.rds")

# cors <- cv_sim_long %>% group_by(method, type) %>%
#   summarise(spearman = cor(bbmsy_true_trans, bbmsy_est_trans, method = "spearman",
#     use = "pairwise.complete.obs")) %>% as.data.frame %>%
#   mutate(spearman = round(spearman, 4)) %>%
#   arrange(type, -spearman)

# ggplot(cors, aes(x = method, xend = method, yend = spearman)) +
#   geom_segment(y = 0, lwd = 1.2) + facet_wrap(~type)

cors %>% tidyr::spread(type, spearman) %>%
  ggplot(aes(mean, slope)) + geom_point() +
  geom_text(aes(label = method), hjust = 0.5, vjust = -0.5) +
  xlab("Spearman correlation of B/BMSY mean") +
  ylab("Spearman correlation of B/BMSY slope")

re <- cv_sim_long %>% mutate(
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

p <- performance %>%
  reshape2::dcast(method + type ~ summary, value.var = c("performance")) %>%
  ggplot(aes(MARE, spearman)) +
  #geom_point() +
  geom_text(aes(label = method), hjust = 0.2, size = 3) +
  facet_wrap(~type, scales = "free") +
  xlab("MARE within stocks") + ylab("Spearman correlation across stocks")
ggsave("figs/performance-sim-scatter.pdf", width = 8, height = 4)
