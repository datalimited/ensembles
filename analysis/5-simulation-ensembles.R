# The main file to run ensemble models for the simulation data.

library("dplyr")
library("ggplot2")
source("0-ensemble-functions.R")

dsim <- readRDS("generated-data/dsim.rds")
spec <- readRDS("generated-data/spec-frequencies-sim.rds")
dsim <- suppressWarnings(left_join(dsim, spec)) # warnings on character-factor conversions
dsim$method_id <- sub("COM.SIR", "COMSIR", dsim$method_id) # to match RAM fits
dsim$method_id <- sub("Costello", "mPRM", dsim$method_id)
dsim <- rename(dsim, stockid = stock_id, method = method_id) # to match RAM fits
# extend the stockid to be truly unique per scenario:
# (note that there are still multiple iterations with the same stockid)
# cv_id defines groups that should *all* be in the testing or training
# datasets (note that stockid and cv_id are now essentially the same
# these remain due to legacy code)
dsim <- dsim %>%
  mutate(cv_id = paste0(stockid, "_sigC_", sigmaC, "_sigR_", sigmaR)) %>%
  mutate(stockid =
      paste0(stockid, "_sigmaC_", sigmaC, "_sigmaR_", sigmaR, "_lh_", LH)) %>%
  arrange(stockid, iter, year) # critical since not all in order

# summarise the mean and slope in the last N years:
library("doParallel")
registerDoParallel(cores = parallel::detectCores())
.parallel <- ifelse(Sys.info()[["sysname"]] == "Windows", FALSE, TRUE)

dsim_sum <- plyr::ddply(dsim, c("stockid", "method", "iter"),
  .parallel = .parallel, .fun = mean_slope_bbmsy)
saveRDS(dsim_sum, file = "generated-data/dsim_sum.rds")
# dsim_sum <- readRDS("generated-data/dsim_sum.rds")

# join in some characteristics that we'll use in models:
dsim_meta <- dsim %>%
  group_by(stockid, iter) %>%
  summarise(
    cv_id = cv_id[1],
    LH = LH[1],
    spec_freq_0.05 = spec_freq_0.05[1],
    spec_freq_0.1 = spec_freq_0.1[1],
    spec_freq_0.2 = spec_freq_0.2[1],
    spec_freq_0.5 = spec_freq_0.5[1])
dsim_sum <- inner_join(dsim_sum, dsim_meta)
dsim_sum$LH <- as.factor(dsim_sum$LH)

# add a binary above/below B/Bmsy = 1 column for a potential response:
dsim_sum <- dsim_sum %>% mutate(
  above_bbmsy1_true = ifelse(bbmsy_true_mean > 1, 1, 0))

# save a data frame of 'true' operating model values to merge in:
trues <- select(dsim_sum, stockid, iter, above_bbmsy1_true,
  bbmsy_true_mean, bbmsy_true_slope)
trues <- trues[!duplicated(trues), ] # one value per operating model stockid

# switch from long to wide format for modelling:
d_mean <- reshape2::dcast(dsim_sum,
  stockid + cv_id + iter + LH + spec_freq_0.05 + spec_freq_0.2 ~ method,
  value.var = "bbmsy_est_mean")  %>%
  inner_join(select(trues, -bbmsy_true_slope))
d_slope <- reshape2::dcast(dsim_sum,
  stockid + cv_id + iter + LH + spec_freq_0.05 + spec_freq_0.2 ~ method,
  value.var = "bbmsy_est_slope") %>%
  inner_join(select(trues, -bbmsy_true_mean))

# strip out those with NAs (in CMSY)
# create problems with some models, such as randomForest, otherwise
d_mean <- na.omit(d_mean)
d_slope <- na.omit(d_slope)
saveRDS(d_mean, file = "generated-data/sim-mean-dat.rds")

# Number of variables in the ensemble models
# gets used to, for example, figure out how many panels to plot
nvar <- 6L

# run a model on all the data to generate data for partial dependence plots:
# m_rf <- randomForest::randomForest(
#   log(bbmsy_true_mean) ~ CMSY + COMSIR + mPRM + SSCOM + LH +
#   spec_freq_0.05 + spec_freq_0.2, data = d_mean)

# best.iter <- gbm::gbm.perf(m_gbm1, method="cv", oobag.curve = FALSE)
# print(best.iter)
# library("caret")
# x <- dplyr::select(d_mean, CMSY, COMSIR, mPRM, SSCOM, spec_freq_0.05, spec_freq_0.2)
# mm_gbm <- train(x = x , y = log(d_mean_sim$bbmsy_true_mean), method = "gbm",
#   trControl = trainControl(method = "cv", number = 3, repeats = 1),
#   tuneGrid =
#     expand.grid(
#       interaction.depth = c(1, 2, 4, 6),
#       n.trees = c(500, 2000, 4000, 8000),
#       shrinkage = c(0.1, 0.01, 0.005),
#       n.minobsinnode = 10),  # has little effect here, use default
#   verbose = FALSE)
# pdf("../figs/gbm-selection-rmse.pdf", width = 8, height = 6)
# plot(mm_gbm)
# dev.off()
#
# pdf("../figs/gbm-selection-rsq.pdf", width = 8, height = 6)
# plot(mm_gbm, metric = "Rsquared")
# dev.off()

mlm <- lm(log(bbmsy_true_mean) ~ CMSY + COMSIR + mPRM + SSCOM +
  spec_freq_0.05 + spec_freq_0.2, data = d_mean)
mlm2 <- lm(log(bbmsy_true_mean) ~ (CMSY + COMSIR + mPRM + SSCOM +
 spec_freq_0.05 + spec_freq_0.2)^2, data = d_mean)
mlm3 <- lm(log(bbmsy_true_mean) ~ (CMSY + COMSIR + mPRM + SSCOM +
 spec_freq_0.05 + spec_freq_0.2)^3, data = d_mean)
bbmle::AICtab(mlm, mlm2, mlm3)

m <- gbm::gbm(log(bbmsy_true_mean) ~ CMSY + COMSIR + mPRM + SSCOM +
    spec_freq_0.05 + spec_freq_0.2, distribution = "gaussian",
  data = d_mean, n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01)
d_mean$gbm_pred <- exp(gbm::predict.gbm(m, n.trees = 2000L))

m_slope <- gbm::gbm(bbmsy_true_slope ~ CMSY + COMSIR + mPRM + SSCOM +
    spec_freq_0.05 + spec_freq_0.2, distribution = "gaussian",
  data = d_slope, n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01)

m_basic <- gbm::gbm(log(bbmsy_true_mean) ~ CMSY + COMSIR + mPRM + SSCOM,
  distribution = "gaussian",
  data = d_mean, n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01)

m_slope_basic <- gbm::gbm(bbmsy_true_slope ~ CMSY + COMSIR + mPRM + SSCOM,
  distribution = "gaussian",
  data = d_slope, n.trees = 2000L, interaction.depth = 6, shrinkage = 0.01)

p1 <- ggplot(d_mean, aes(CMSY, gbm_pred)) + geom_point()
p2 <- ggplot(d_mean, aes(COMSIR, gbm_pred)) + geom_point()
p3 <- ggplot(d_mean, aes(mPRM, gbm_pred)) + geom_point()
p4 <- ggplot(d_mean, aes(SSCOM, gbm_pred)) + geom_point()
p5 <- ggplot(d_mean, aes(spec_freq_0.05, gbm_pred)) + geom_point()
p6 <- ggplot(d_mean, aes(spec_freq_0.2, gbm_pred)) + geom_point()
pdf("../figs/predictor-vs-predicted-gbm.pdf", width = 10, height = 5)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)
dev.off()

p1 <- ggplot(d_mean, aes(CMSY, bbmsy_true_mean)) + geom_point()
p2 <- ggplot(d_mean, aes(COMSIR, bbmsy_true_mean)) + geom_point()
p3 <- ggplot(d_mean, aes(mPRM, bbmsy_true_mean)) + geom_point()
p4 <- ggplot(d_mean, aes(SSCOM, bbmsy_true_mean)) + geom_point()
p5 <- ggplot(d_mean, aes(spec_freq_0.05, bbmsy_true_mean)) + geom_point()
p6 <- ggplot(d_mean, aes(spec_freq_0.2, bbmsy_true_mean)) + geom_point()
pdf("../figs/predictor-vs-truth-gbm.pdf", width = 10, height = 5)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)
dev.off()

partial <- plyr::ldply(seq_len(4L), function(i) {
  dd <- gbm::plot.gbm(m_basic, i.var = i, return.grid = TRUE)
  dd$predictor <- names(dd)[1]
  names(dd)[1] <- "predictor_value"
  dd$y <- exp(dd$y)
  dd
})

partial_slope <- plyr::ldply(seq_len(4L), function(i) {
  dd <- gbm::plot.gbm(m_slope_basic, i.var = i, return.grid = TRUE)
  dd$predictor <- names(dd)[1]
  names(dd)[1] <- "predictor_value"
  dd
})

partial_names <- data_frame(
  predictor = c("CMSY", "COMSIR", "mPRM", "SSCOM",
    "spec_freq_0.05", "spec_freq_0.2"),
  predictor_label = c("(a) CMSY", "(b) COM-SIR", "(d) mPRM", "(c) SSCOM",
    "(f) Long-term spectral density", "(e) Short-term spectral density"))
partial <- inner_join(partial, partial_names)
partial_slope <- inner_join(partial_slope, partial_names)

# partial dependence plot:
p <- ggplot(partial, aes(predictor_value, y)) + geom_line() +
  facet_wrap(~predictor_label, scales = "free_x") + ylim(0.3, 1.6) + theme_bw() +
  ylab(expression(Average~widehat(B/B[MSY]))) +
  xlab("Predictor value") +
  xlim(0, 2) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"))
ggsave("../figs/partial-sim.pdf", width = 7, height = 5)

p <- ggplot(partial_slope, aes(predictor_value, y)) + geom_line() +
  facet_wrap(~predictor_label, scales = "free_x") + ylim(-0.2, 0.2) + theme_bw() +
  ylab(expression(Average~predicted~slope~of~B/B[MSY])) +
  xlab("Predictor value") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white")) +
      xlim(-0.4, 0.4)
ggsave("../figs/partial-sim-slope.pdf", width = 7, height = 5)

partial_2d <- plyr::ldply(1:4, function(x) plyr::ldply(1:4, function(y) {
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
zlim <- c(min(partial_2d$z)+0.21, max(partial_2d$z))
pal <- rev(colorRampPalette(c("red", "white", "blue"), space = "Lab")(13))[-c(1:1)]
# white should line up with 1: (adjust the above 2 lines as needed)

pdf("../figs/partial-2d-col-check.pdf", width = 6, height = 4)
plot(seq(min(zlim), max(zlim), length.out = 12), 1:12, bg = pal, pch = 21, cex = 3);abline(v = 1)
dev.off()

nvar_basic <- 4L
ii <<- 0
lab <<- 0
panels <- matrix(NA, ncol = nvar_basic, nrow = nvar_basic, byrow = TRUE) %>% upper.tri %>%
  t %>% as.vector

pdf("../figs/partial-sim-2d.pdf", width = 7, height = 5.5)
par(mfrow = c(nvar_basic-1, nvar_basic), mar = c(3.5,3.5,1,1), oma = c(1, 1, 0.5, 0.5), cex = 0.7)
par(xpd = NA, mgp = c(1.7, 0.5, 0), tck = -0.03, las = 1)
plyr::d_ply(partial_2d, c("var1", "var2"), function(x) {
  ii <<- ii + 1
  if (panels[ii]) {
    xx <- reshape2::dcast(x, x ~ y, value.var = "z")
    image(xx[,1], as.numeric(colnames(xx)[-1]), as.matrix(xx[,-1]), main = "",
      xlab = unique(x$var1), ylab = unique(x$var2),
      zlim = zlim,
      col = pal)
    lab <<- lab + 1
    add_label(0, -0.08, label = paste0("(", letters[lab], ")"))
  } else {
    if (ii < nvar_basic^2 - nvar_basic) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
      if (ii == 1) {
        blocks <- seq(0.1, 0.9, length.out = length(pal))
        for(i in seq_along(pal)) {
          rect(xleft = 0.9, ybottom = blocks[i],
               xright = 1.1, ytop = blocks[i+1]+0.001,
               border = NA, col = pal[i])
        }
        text(0.9, blocks[1], round(min(zlim), 1), pos = 2)
        text(0.9, 0.5, 1, pos = 2)
        text(0.9, blocks[length(blocks)], round(max(zlim), 1), pos = 2)
        mtext(expression(widehat(B/B[MSY])), side = 2, las = 0, cex = 0.8)
      }
    }
  }
})
dev.off()
message(paste("zlim were", round(zlim, 2), collapse = " "))

# the base form of the ensemble models:
eq <- paste0("CMSY + COMSIR + mPRM + ",
  "SSCOM + spec_freq_0.05 + spec_freq_0.2")
eq_basic <- paste0("CMSY + COMSIR + mPRM + SSCOM")

# work through cross validation of ensemble models:
cv_sim_mean <- plyr::ldply(seq_len(50), .parallel = .parallel,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE, id = "sim-mean",
      gbm_formula = paste0("log(bbmsy_true_mean) ~ ", eq),
      lm_formula = paste0("log(bbmsy_true_mean) ~ (", eq, ")^2"), weighted = TRUE))
cv_sim_mean <- cv_sim_mean %>% mutate(
  gbm_ensemble = exp(gbm_ensemble),
  rf_ensemble  = exp(rf_ensemble),
  lm_ensemble  = exp(lm_ensemble))
cv_sim_mean$cv_id <- NULL

cv_sim_mean_basic <- plyr::ldply(seq_len(50), .parallel = .parallel,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean, geo_mean = TRUE, id = "sim-mean",
      gbm_formula = paste0("log(bbmsy_true_mean) ~ ", eq_basic),
      lm_formula = paste0("log(bbmsy_true_mean) ~ (", eq_basic, ")^2"), weighted = TRUE))
cv_sim_mean_basic <- cv_sim_mean_basic %>% mutate(
  gbm_ensemble = exp(gbm_ensemble),
  rf_ensemble  = exp(rf_ensemble),
  lm_ensemble  = exp(lm_ensemble))
cv_sim_mean_basic$cv_id <- NULL

cv_sim_slope <- plyr::ldply(seq_len(50), .parallel = .parallel,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_slope, geo_mean = FALSE, id = "sim-slope",
      gbm_formula = paste0("bbmsy_true_slope ~ ", eq_basic),
      lm_formula = paste0("bbmsy_true_slope ~ (", eq_basic, ")^2")))
cv_sim_slope$cv_id <- NULL

cv_sim_binary <- plyr::ldply(seq_len(1), .parallel = .parallel,
  .fun = function(.n)
    cross_val_ensembles(.n = .n, dat = d_mean,
      id = "sim-binary", distribution = "bernoulli",
      gbm_formula = paste0("above_bbmsy1_true ~ ", eq_basic),
      glm_formula = paste0("above_bbmsy1_true ~ (", eq_basic, ")^2")))
cv_sim_binary$cv_id <- NULL

# now switch to long format data, summarize, and compare:
cv_sim_binary_long <- cv_sim_binary %>%
  select(-spec_freq_0.05, -spec_freq_0.2, -bbmsy_true_mean) %>%
  reshape2::melt(id.vars = c("stockid", "iter", "test_iter", "LH", "above_bbmsy1_true"),
    variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = above_bbmsy1_true) %>%
  mutate(type = "binary")
saveRDS(cv_sim_binary_long, file = "generated-data/cv_sim_binary.rds") # used in 5-roc.R

cv_sim_mean_basic_long <- cv_sim_mean_basic %>%
  select(-spec_freq_0.05, -spec_freq_0.2, -above_bbmsy1_true) %>%
  reshape2::melt(id.vars = c("stockid", "iter", "test_iter", "LH", "bbmsy_true_mean"),
    variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean) %>%
  mutate(
    type = "mean",
    bbmsy_true_trans = log(bbmsy_true),
    bbmsy_est_trans = log(bbmsy_est))

cv_sim_mean_long <- cv_sim_mean %>%
  select(-spec_freq_0.05, -spec_freq_0.2, -above_bbmsy1_true) %>%
  reshape2::melt(id.vars = c("stockid", "iter", "test_iter", "LH", "bbmsy_true_mean"),
    variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_mean) %>%
  mutate(
    type = "mean",
    bbmsy_true_trans = log(bbmsy_true),
    bbmsy_est_trans = log(bbmsy_est))

cv_sim_slope_long <- cv_sim_slope %>%
  select(-spec_freq_0.05, -spec_freq_0.2, -above_bbmsy1_true) %>%
  reshape2::melt(id.vars = c("stockid", "iter", "test_iter", "LH", "bbmsy_true_slope"),
    variable.name = "method", value.name = "bbmsy_est") %>%
  rename(bbmsy_true = bbmsy_true_slope) %>%
  mutate(
    type = "slope",
    bbmsy_true_trans = bbmsy_true,
    bbmsy_est_trans = bbmsy_est)

cv_sim_long <- suppressWarnings(
  dplyr::bind_rows(cv_sim_mean_long, cv_sim_slope_long))

re <- cv_sim_long %>% mutate(
  sq_er = (bbmsy_est - bbmsy_true)^2,
  re    = (bbmsy_est - bbmsy_true) / abs(bbmsy_true))

re %>% saveRDS(file = "generated-data/cv_sim_long.rds")

re <- cv_sim_mean_basic_long %>% mutate(
  sq_er = (bbmsy_est - bbmsy_true)^2,
  re    = (bbmsy_est - bbmsy_true) / abs(bbmsy_true))

re %>% saveRDS(file = "generated-data/cv_sim_mean_basic_long.rds")
