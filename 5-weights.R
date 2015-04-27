# Build weighted ensemble models and evaluate performance

library("dplyr")
library("reshape2")
library("ggplot2")
library("gbm")
library("datalimited")

SHRINKAGE <- 0.001
INT.DEPTH <- 4
N.TREES <- 2000
filename_base <- paste(N.TREES, SHRINKAGE, INT.DEPTH, sep = "-")

ram_fits <- readRDS("generated-data/ram_fits.rds") %>%
  rename(sigma_effort = sigma, slope_effort = slope)

# Try fitting a model to the time series B/Bmsy data:
m_gbm <- gbm(log(b2bmsy_true) ~
    year_before_end + COMSIR + Costello + SSCOM + CMSY_new_prior +
    resilience + sigma_effort + slope_effort +
    habitat + log_max_catch + log_total_catch +
    envtemp,
  data = ram_fits, distribution = "gaussian", n.trees = N.TREES,
  interaction.depth = INT.DEPTH, shrinkage = SHRINKAGE)

factor_variables <- c(6:9, 12)
pdf(paste0("figs/gbm-partial-", filename_base, ".pdf"), width = 11, height = 7)
par(mfrow = c(3, 4))
par(oma = c(2, 5, 1, 1))
par(mar = c(4, 6, 1, 1))
par(cex = 0.8)
for(i in 1:12) {
  if(i %in% factor_variables) {
    x <- plot(m_gbm, i.var = i, xlim = c(-0.8, 0), return.grid = TRUE)
    plot(exp(x[,2]), x[,1], yaxt = "n", ylab = "", xlab = "", pch = 19, xlim = exp(c(-0.8, 0)))
    axis(2, at = x[,1], labels = x[,1], las = 1, cex.axis = 1.3)
  } else {
    x <- plot(m_gbm, i.var = i, ylim = c(-0.8, 0), xlab = "", return.grid = TRUE)
    plot(x[,1], exp(x[,2]), ylab = "", xlab = "", type = "l", ylim = exp(c(-0.8, 0)))
  }
    mtext(names(x)[1], side = 1, line = 2.2)
}
dev.off()

pdf("figs/partial-eg-effort-sscom.pdf", width = 5, height = 4)
plot(m_gbm, i.var = c(4, 7))
dev.off()

pdf("figs/partial-eg-comsir-costello.pdf", width = 5, height = 4)
plot(m_gbm, i.var = c(2, 3))
dev.off()

ram_fits$gbm_ensemble <- exp(predict(m_gbm, n.trees = m_gbm$n.trees))
ram_fits <- plyr::adply(ram_fits, 1, function(x)
  data.frame(mean_ensemble = exp(mean(c(log(x$CMSY_new_prior),
    log(x$Costello), log(x$SSCOM), log(x$COMSIR))))))

ram_fits$one <- rnorm(nrow(ram_fits), mean = 1, sd = 0.05)

png("figs/pairs-models-ram.png", width = 700, height = 700)
pairs(log(ram_fits[,c("b2bmsy_true", "CMSY_new_prior", "COMSIR", "Costello", "SSCOM", "gbm_ensemble", "mean_ensemble")]), col = "#00000010")
dev.off()

# Now switch back to long format data and calculate goodness of fit stats:
ram <- melt(ram_fits, measure.vars = c("CMSY_new_prior", "COMSIR",
  "Costello", "SSCOM", "gbm_ensemble", "mean_ensemble", "one"),
  value.name = "b2bmsy", variable.name = "method")

# relative error in slope, centered intercept; correlation of 10 years
lm_re <- plyr::ddply(ram, c("stockid", "method"), function(x) {
  centered_year <- x$year_before_end - mean(x$year_before_end)
  m <- lm(log(x$b2bmsy) ~ centered_year)
  slope_est <- coef(m)[[2]]
  intercept_est <- coef(m)[[1]]

  m <- lm(log(x$b2bmsy_true) ~ centered_year)
  slope_true <- coef(m)[[2]]
  intercept_true <- coef(m)[[1]]

  data.frame(slope_b2bmsy = slope_est, int_b2bmsy = intercept_est,
    slope_b2bmsy_true = slope_true, int_b2bmsy_true = intercept_true)
})

lm_re <- lm_re %>% mutate(
  slope_re = (slope_b2bmsy - slope_b2bmsy_true) / slope_b2bmsy_true,
  int_re = (int_b2bmsy - int_b2bmsy_true) / int_b2bmsy_true)

group_by(lm_re, method) %>% summarise(mare_slope = median(abs(slope_re)),
  mare_int = median(abs(int_re)))

# time series plots:
widths <- data.frame(method = c("CMSY_new_prior", "COMSIR", "Costello",
  "SSCOM", "gbm_ensemble", "mean_ensemble", "one"),
  type = c("1-individual", "1-individual", "1-individual",
    "1-individual", "2-ensemble", "2-ensemble", "2-ensemble"))
ram$widths <- NULL
ram <- inner_join(ram, widths)

# pdf("figs/ensemble-ts-ram.pdf", width = 25, height = 25)
# q <- ggplot(subset(ram, method != "one"),
#   aes(year_before_end, b2bmsy, colour = method)) +
#   geom_line(aes(size = type)) + scale_size_manual(values = c(0.5, 1.5)) +
#   geom_line(aes(year_before_end, b2bmsy_true), lwd = 2, colour = "black") +
#   geom_hline(yintercept = 1, lty = 2) +
#   scale_y_continuous(limits = c(0, 3)) +
#   theme_bw() + facet_wrap(~stockid) +
#   theme(plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()) +
#   xlab("Year (last 10 years)") +
#   ylab(expression(B/B[MSY]))
# print(q)
# dev.off()

# check the correlation of estimated vs. "true":
cors <- group_by(ram, method, stockid) %>%
  summarise(cor = cor(b2bmsy, b2bmsy_true)) %>%
  as.data.frame()

p1 <- ggplot(lm_re, aes(method, slope_re)) + geom_boxplot() + ylim(-4, 4)
p2 <- ggplot(lm_re, aes(method, int_re)) + geom_boxplot() + ylim(-4, 4)
p3 <- ggplot(cors, aes(method, cor)) + geom_boxplot() + ylim(-1, 1)
pdf("figs/slope-intercept-re-ram-ensembles.pdf", width = 7, height = 9)
gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
dev.off()

## bring in the ensemble model from the simulation output:

# a temporary version with properly named columns:
ram_fits_temp <- ram_fits %>% dplyr::rename(CMSY = CMSY_new_prior, COM.SIR = COMSIR)

# and predict on RAM data with the simulation-weighted ensemble model:
m_gbm_sim <- readRDS("generated-data/m_gbm_sim.rds")
ram_fits$gbm_sim_ensemble <- gbm::predict.gbm(m_gbm_sim,
  newdata = ram_fits_temp, n.trees = m_gbm_sim$n.trees)

######################
# now do the above, but each time subset the data (cross validate) and record
# the output (relative error value)

# build PRM-formatted data once to save time:
# ram_ts and spp_categories are in the datalimited package data:
prm_dat <- inner_join(ram_ts, spp_categories, by = "scientificname")
prm_dat <- filter(prm_dat, stocklong %in% ram_fits$stocklong)
prm_dat <- plyr::ddply(prm_dat, "stockid", function(x) {
  format_prm(year = x$year, catch = x$catch, bbmsy = x$bbmsy_ram,
    species_cat = x$spp_category[1L])
})

library("doParallel")
registerDoParallel(cores = 4)

if(!file.exists("cv-ensembles")) system("mkdir cv-ensembles")
system("rm -r cv-ensembles/*")
if(!file.exists("cv-gbms")) system("mkdir cv-gbms")
system("rm -r cv-gbms/*")
if(!file.exists("cv-b2bmsy-lm")) system("mkdir cv-b2bmsy-lm")
system("rm -r cv-b2bmsy-lm/*")

set.seed(1)
cv_weights <- plyr::ldply(seq_len(4), .parallel = TRUE, .fun = function(.n) {

  cv_ids_set <- FALSE

  while(!cv_ids_set) {
    nstocks <- length(unique(ram_fits$stockid))

    train_ids <- sample(nstocks, round(nstocks*0.5))
    test_ids <- seq_len(nstocks)[-train_ids]

    train_stock_ids <- unique(ram_fits$stockid)[train_ids]
    test_stock_ids <- unique(ram_fits$stockid)[test_ids]

    train_dat <- filter(ram_fits, stockid %in% train_stock_ids)
    test_dat <- filter(ram_fits, stockid %in% test_stock_ids)

    # need to buid new Costello-style model each time:
    train_prm_dat <- filter(prm_dat, stockid %in% train_stock_ids)
    test_prm_dat <- filter(prm_dat, stockid %in% test_stock_ids)
    train_prm_dat <- filter(train_prm_dat, !is.na(bbmsy)) # or gbm() fails

    # make sure we have all 'test' species categories in our 'train' data
    # or the lm() prediction will fail
    if (all(unique(test_prm_dat$species_cat) %in% unique(train_prm_dat$species_cat)))
      cv_ids_set <- TRUE
  }

  mprm <- fit_prm(train_prm_dat)
  mprm_gbm <- fit_prm(train_prm_dat, type = "gbm", n.trees = N.TREES,
    interaction.depth = INT.DEPTH,
    shrinkage = SHRINKAGE)
  # might fail if factor levels don't match on train/test:
  test_prm_dat$Costello <- tryCatch({predict_prm(test_prm_dat, model = mprm)},
    error = function(e) rep(NA, nrow(test_prm_dat)))
  test_prm_dat$Costello_gbm <- tryCatch({predict_prm(test_prm_dat, model = mprm_gbm)},
    error = function(e) rep(NA, nrow(test_prm_dat)))
  test_prm_dat$tsyear <- test_prm_dat$year

  # now sub back CV Costello results into main data frame:
  test_dat$Costello <- NULL
  test_dat <- inner_join(test_dat,
    test_prm_dat[,c("stockid", "tsyear", "Costello", "Costello_gbm")],
    by = c("stockid", "tsyear"))

  m_gbm <- gbm(log(b2bmsy_true) ~
      year_before_end + COMSIR + Costello + SSCOM + CMSY_new_prior +
      resilience + sigma_effort + slope_effort + habitat +
      log_max_catch + log_total_catch +
      envtemp,
    data = train_dat, distribution = "gaussian", n.trees = N.TREES,
    interaction.depth = INT.DEPTH, shrinkage = SHRINKAGE)
  # save the model to extract form later
  saveRDS(m_gbm, file = paste0("cv-gbms/", .n, "-gbm.rds"))

  library("randomForest")
  m_rf <- randomForest(log(b2bmsy_true) ~
      year_before_end + COMSIR + Costello + SSCOM + CMSY_new_prior +
      resilience + sigma_effort + slope_effort + habitat +
      log_max_catch + log_total_catch +
      envtemp, data = train_dat)

  test_dat$gbm_ensemble <- tryCatch({exp(predict(m_gbm,
    n.trees = N.TREES, newdata = test_dat))},
    error = function(e) rep(NA, length(test_ids)))

  test_dat$rf_ensemble <- tryCatch({exp(predict(m_rf, newdata = test_dat))},
    error = function(e) rep(NA, length(test_ids)))

  test_dat <- plyr::adply(test_dat, 1, function(x)
    data.frame(mean_ensemble = exp(mean(c(log(x$CMSY_new_prior),
      log(x$Costello), log(x$SSCOM), log(x$COMSIR))))))

  test_dat$one <- rnorm(nrow(test_dat), mean = 1, sd = 0.05)

  # save for after:
  test_dat$.n <- .n
  saveRDS(test_dat, file = paste0("cv-ensembles/", .n, ".rds"))

  ram <- melt(test_dat, measure.vars = c("CMSY_new_prior", "COMSIR",
    "Costello", "SSCOM", "gbm_ensemble", "mean_ensemble", "Costello_gbm",
    "gbm_sim_ensemble", "rf_ensemble", "one"),
    value.name = "b2bmsy", variable.name = "method")

  # relative error in slope, centered intercept; correlation of 10 years
  lm_out <- plyr::ddply(ram, c("stockid", "method"), function(x) {
    centered_year <- x$year_before_end - mean(x$year_before_end)
    m <- lm(log(x$b2bmsy) ~ centered_year)
    slope_est <- coef(m)[[2]]
    intercept_est <- coef(m)[[1]]
    m <- lm(log(x$b2bmsy_true) ~ centered_year)
    slope_true <- coef(m)[[2]]
    intercept_true <- coef(m)[[1]]
    data.frame(slope_b2bmsy = slope_est, int_b2bmsy = intercept_est,
      slope_b2bmsy_true = slope_true, int_b2bmsy_true = intercept_true)
  })

  # save the slopes and intercepts to plot later:
  saveRDS(lm_out, file = paste0("cv-b2bmsy-lm/", .n, ".rds"))

  lm_re <- lm_out %>% mutate(
    slope_re = (slope_b2bmsy - slope_b2bmsy_true) / slope_b2bmsy_true,
    int_re = (int_b2bmsy - int_b2bmsy_true) / int_b2bmsy_true,
    slope_se = (slope_b2bmsy - slope_b2bmsy_true)^2,
    int_se = (int_b2bmsy - int_b2bmsy_true)^2)

  # check the correlation of estimated vs. "true":
  cors_within <- group_by(ram, method, stockid) %>%
    summarise(cor_within = cor(b2bmsy, b2bmsy_true)) %>%
    as.data.frame

  cors_across <- group_by(lm_out, method) %>%
    summarise(
      cor_slope_across = cor(slope_b2bmsy, slope_b2bmsy_true),
      cor_int_across = cor(int_b2bmsy, int_b2bmsy_true))

  out <- group_by(lm_re, method) %>%
    summarise(mare_slope = median(abs(slope_re), na.rm = TRUE),
      mare_int = median(abs(int_re), na.rm = TRUE),
      rmse_slope = sqrt(mean(slope_se, na.rm = TRUE)),
      rmse_int = sqrt(mean(int_se, na.rm = TRUE))) %>%
    as.data.frame
  out$.n <- .n
  out <- inner_join(out, cors_within, by = "method")
  out <- inner_join(out, cors_across, by = "method")
  out
})
saveRDS(cv_weights, file = "generated-data/cv_weights.rds") # cache results
cv_weights <- readRDS("generated-data/cv_weights.rds")

cv_weights <- filter(cv_weights, method != "one")
cv_weights$ensemble <- grepl("ensemble", cv_weights$method)

# gather the time series:
cv_ts_files <- list.files("cv-ensembles", full.names = TRUE)
cv_ts_dat <- list()
for (i in seq_along(cv_ts_files)) {
  cv_ts_dat[[i]] <- readRDS(cv_ts_files[i])
}
cv_ts_dat <- dplyr::rbind_all(cv_ts_dat) %>% as.data.frame

# gather the relative influence values:
cv_gbms <- list.files("cv-gbms", full.names = TRUE)
cv_gbms_dat <- list()
for (i in seq_along(cv_gbms)) {
  cv_gbms_dat[[i]] <- readRDS(cv_gbms[i])
}
gbm_ri <- cv_gbms_dat %>% plyr::ldply(function(x) {
  q <- gbm::relative.influence(x, n.trees = x$n.trees)
  q / max(q)
})
gbm_ri_long <- reshape2::melt(gbm_ri)
gbm_ri_order <- group_by(gbm_ri_long, variable) %>%
  summarise(mean_ri = mean(value)) %>%
  as.data.frame() %>%
  mutate(order = rank(-mean_ri))
gbm_ri_long <- inner_join(gbm_ri_long, gbm_ri_order) %>%
  mutate(parameter = reorder(variable, -order))

p <- ggplot(gbm_ri_long, aes(parameter, value)) + geom_boxplot() +
  coord_flip() + xlab("") + ylab("Scaled relative importance")
  #geom_point(alpha = 0.5, position = position_jitter(width = 0.1))
ggsave("figs/relative-importance-cv.pdf", width = 4, height = 4.5)
# end relative importance plots

# now try plotting slopes and intercepts across all cv reps
cv_slope_int_files <- list.files("cv-b2bmsy-lm/", full.names = TRUE)
cv_slope_int <- list()
for (i in seq_along(cv_slope_int_files)) {
  cv_slope_int[[i]] <- readRDS(cv_slope_int_files[i])
}
cv_slope_int <- dplyr::rbind_all(cv_slope_int) %>% as.data.frame %>%
  filter(!method %in% c("one", "Costello_gbm")) %>%
  mutate(int_b2bmsy_true = exp(int_b2bmsy_true), int_b2bmsy = exp(int_b2bmsy))

slope_int_plot <- function(x, y, type = "", lim = c(-0.25, 0.25)) {
 ggplot(cv_slope_int, aes_string(x, y)) +
   geom_point(alpha = 0.1) +
  facet_wrap(~method) +
  xlab(paste("B/Bmsy", type, "assessment")) +
  ylab(paste("B/Bmsy", type, "data-limited")) +
  stat_smooth(method = "loess") + coord_fixed(ratio = 1) +
  geom_abline(slope = 1, intercept = 0) + ylim(lim) + xlim(lim)
}

p1 <-slope_int_plot("slope_b2bmsy_true", "slope_b2bmsy", "slope")
ggsave("figs/cv-slope-pairs.pdf", width = 6, height = 6)

p2 <-slope_int_plot("int_b2bmsy_true", "int_b2bmsy", "intercept", lim = c(0, 3))
ggsave("figs/cv-int-pairs.pdf", width = 6, height = 6)
## end plotting slopes and intercepts from CV

cv_weights_long <- cv_weights %>% select(-.n) %>%
  reshape2::melt(id.vars = c("method", "stockid", "ensemble"),
    variable.name = "metric") %>%
  filter(!metric %in% c("mare_slope", "mare_int"))

boxplot_performance <- function(dat, metrics, ylab = "value") {
  p99 <- filter(cv_weights_long, metric %in% metrics) %>%
    ggplot(aes(method, value, fill = ensemble)) +
    geom_boxplot() + xlab("") + ylab(ylab) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p99)
}

pdf(paste0("figs/cv-ensemble-boxplots-", filename_base, ".pdf"),
  width = 11.5, height = 11)
boxplot_performance(dat,
  metrics = c("rmse_slope", "rmse_int", "cor_slope_across", "cor_int_across"))
dev.off()

## make smaller version for slides
pdf(paste0("figs/cv-ensemble-boxplots-rmse-", filename_base, ".pdf"),
  width = 6, height = 5)
boxplot_performance(dat, ylab = "RMSE within stocks",
  metrics = c("rmse_slope", "rmse_int"))
dev.off()

pdf(paste0("figs/cv-ensemble-boxplots-cor-", filename_base, ".pdf"),
  width = 6, height = 5)
boxplot_performance(dat, ylab = "Correlation across stocks",
  metrics = c("cor_slope_across", "cor_int_across"))
dev.off()

# time series plots with cv-ensembles:
# -----------------------------------

# switch back to long format data:
cv_ts_dat_long <- melt(cv_ts_dat, measure.vars = c("CMSY_new_prior", "COMSIR",
  "Costello", "SSCOM", "Costello_gbm", "gbm_sim_ensemble", "gbm_ensemble",
  "mean_ensemble", "rf_ensemble", "one"),
  value.name = "b2bmsy", variable.name = "method")
cv_ts_dat_long <- cv_ts_dat_long %>% group_by(stockid, year_before_end, method) %>%
  summarise(median_b2bmsy = median(b2bmsy),
    q25_b2bmsy = quantile(b2bmsy, probs = 0.25),
    q75_b2bmsy = quantile(b2bmsy, probs = 0.75)) %>%
  mutate(type = ifelse(grepl("ensemble", method), "2-ensemble", "1-individual"))


big_ts_plot <- function(dat, facet = facet_wrap(~stockid)) {
  suppressWarnings({ # geom_path warnings on empty rows
    # q <- ggplot(subset(cv_ts_dat_long, !method %in% c("one", "Costello_gbm") & type == "2-ensemble"),
    q <- ggplot(dat,
      aes(year_before_end, median_b2bmsy, colour = method)) +
      geom_ribbon(aes(ymin = q25_b2bmsy, ymax = q75_b2bmsy, fill = method), alpha = 0.2, lwd = 0) +
      geom_line(lwd = 1) +
      geom_line(aes(size = type)) + scale_size_manual(values = c(0.3, 1.8)) +
      geom_line(data = filter(ram_fits, stockid %in% unique(dat$stockid)),
    aes(year_before_end, b2bmsy_true), lwd = 2, colour = "black") +
      geom_hline(yintercept = 1, lty = 2) +
      scale_y_continuous(limits = c(0, 3)) +
      theme_bw() + facet +
      theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
      xlab("Year (last 10 years)") +
      ylab(expression(B/B[MSY]))
    print(q)
  })
}

pdf(paste0("figs/cv-ensemble-ts-ram-", filename_base, ".pdf"),
  width = 25, height = 25)
big_ts_plot(subset(cv_ts_dat_long, !method %in% c("one", "Costello_gbm")))
dev.off()

## make a smaller version for slides

stocks <- c("HADROCK", "SEELNSSA1", "YELL3LNO", "PHAKEPCOAST")

pdf(paste0("figs/cv-ensemble-eg-ts-ram-", filename_base, ".pdf"), width = 10, height = 6)
big_ts_plot(filter(cv_ts_dat_long,
    !method %in% c("one", "Costello_gbm"),
    stockid %in% stocks))
dev.off()

pdf(paste0("figs/cv-ensemble-eg-basic-ts-ram-", filename_base, ".pdf"), width = 10, height = 6)
big_ts_plot(filter(cv_ts_dat_long,
    !method %in% c("one", "Costello_gbm"),
    type == "2-ensemble" | method == "Costello",
    stockid %in% stocks))
dev.off()
pdf(paste0("figs/cv-eg-basic-ts-ram-", filename_base, ".pdf"), width = 10, height = 6)
big_ts_plot(filter(cv_ts_dat_long,
  !method %in% c("one", "Costello_gbm"),
  type == "1-individual",
  stockid %in% stocks))
dev.off()


filter(ram_fits, stockid %in% stocks) %>%
  group_by(stockid) %>%
  summarise(name = stocklong[1])

# make a plot to motivate ensembles
p <- cv_ts_dat_long %>%
  filter(stockid == "ALBAIO", type != "2-ensemble",
!method %in% c("one", "Costello_gbm")) %>%
  ggplot(aes(year_before_end, median_b2bmsy, colour = method)) +
  geom_line(lwd = 1.4) +
  theme_bw() +
theme(plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()) +
  geom_line(data = filter(ram_fits, stockid %in% "ALBAIO"),
    aes(year_before_end, b2bmsy_true), lwd = 2.4, lty = 1,
    colour = "black", alpha = 0.5) + ylab(expression(B/B[MSY])) +
      xlab("Year (last 10 years)")
ggsave("figs/ensemble-motivation.pdf", width = 7, height = 4)


## demo gbms:
x <- seq(1, 10, length.out = 50)
y <- c(seq(5, 25, length.out = 25), seq(25, 5, length.out = 25)) + rnorm(50, mean = 0, sd = 2.5)
trees <- c(1, 5, 100, 1000, 20000)
plot(x, y)
o <- plyr::llply(trees, function(i) {
  m <- gbm(y ~ x, n.trees = i, shrinkage = 0.001)
  m
})

par(mfrow = c(1, 5))
for(i in seq_along(o)) {
  qq <- plot(o[[i]], return.grid = TRUE)
  plot(x, y, main = trees[i])
  lines(qq)
}
