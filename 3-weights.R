# Build weighted ensemble models and evaluate performance

library("dplyr")
library("reshape2")
library("ggplot2")
library("gbm")

ram_fits <- readRDS("generated-data/ram_fits.rds")

# Try fitting a model to the time series B/Bmsy data:
m_gbm <- gbm(b2bmsy_true ~
    year_before_end + COMSIR + Costello + SSCOM + CMSY_new_prior +
    resilience + sigma + slope + habitat,
  data = ram_fits, distribution = "gaussian", n.trees = 5000, interaction.depth = 4,
  shrinkage = 0.0005)

pdf("figs/gbm-partial.pdf", width = 8, height = 8)
par(cex = 0.6)
# par(mfrow = c(3, 3));for(i in 1:9) plot(m_gbm, i.var = i, ylim = c(-1, 0))
par(mfrow = c(3, 3));for(i in 1:9) plot(m_gbm, i.var = i, ylim = c(0.8, 1.2))
dev.off()

# plot(m_gbm, i.var = c(3,2), ylim = c(0.8, 1.2))

ram_fits$gbm_ensemble <- (predict(m_gbm, n.trees = 5000))
ram_fits <- plyr::adply(ram_fits, 1, function(x)
  data.frame(mean_ensemble = exp(mean(c(log(x$CMSY_new_prior),
    log(x$Costello), log(x$SSCOM), log(x$COMSIR))))))

ram_fits$one <- rnorm(nrow(ram_fits), mean = 1, sd = 0.05)

pdf("figs/pairs-models-ram.pdf", width = 9, height = 9)
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

pdf("figs/ensemble-ts-ram.pdf", width = 25, height = 25)
q <- ggplot(subset(ram, method != "one"),
  aes(year_before_end, b2bmsy, colour = method)) +
  geom_line(aes(size = type)) + scale_size_manual(values = c(0.5, 1.5)) +
  geom_line(aes(year_before_end, b2bmsy_true), lwd = 2, colour = "black") +
  geom_hline(yintercept = 1, lty = 2) +
  scale_y_continuous(limits = c(0, 3)) +
  theme_bw() + facet_wrap(~stockid) +
  theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Year (last 10 years)") +
  ylab(expression(B/B[MSY]))
print(q)
dev.off()

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

######################
# now do the above, but each time subset the data and record the output (RE value)

library("doParallel")
registerDoParallel(cores = 4)

cv_weights <- plyr::ldply(seq_len(16), .parallel = TRUE, .fun = function(.n) {

  train_ids <- sample(seq_len(nrow(ram_fits)), round(nrow(ram_fits)*0.5))
  test_ids <- seq_len(nrow(ram_fits))[-train_ids]

  train_dat <- ram_fits[train_ids, ]
  test_dat <- ram_fits[test_ids, ]

  m_gbm <- gbm(b2bmsy_true ~
      year_before_end + COMSIR + Costello + SSCOM + CMSY_new_prior +
      resilience + sigma + slope + habitat,
    data = train_dat, distribution = "gaussian", n.trees = 5000,
    interaction.depth = 4, shrinkage = 0.0005)

  test_dat$gbm_ensemble <- tryCatch({predict(m_gbm,
    n.trees = 5000, newdata = test_dat)},
    error = function(e) rep(NA, length(test_ids)))

  test_dat <- plyr::adply(test_dat, 1, function(x)
    data.frame(mean_ensemble = exp(mean(c(log(x$CMSY_new_prior),
      log(x$Costello), log(x$SSCOM), log(x$COMSIR))))))

  test_dat$one <- rnorm(nrow(test_dat), mean = 1, sd = 0.05)

  ram <- melt(test_dat, measure.vars = c("CMSY_new_prior", "COMSIR",
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

  # check the correlation of estimated vs. "true":
  cors <- group_by(ram, method, stockid) %>%
    summarise(cor = cor(b2bmsy, b2bmsy_true)) %>%
    as.data.frame

  out <- group_by(lm_re, method) %>%
    summarise(mare_slope = median(abs(slope_re), na.rm = TRUE),
      mare_int = median(abs(int_re), na.rm = TRUE)) %>% as.data.frame
  out$.n <- .n
  out <- inner_join(out, cors)
  out
})

p1 <- ggplot(cv_weights, aes(method, mare_slope)) + geom_boxplot()
p2 <- ggplot(cv_weights, aes(method, mare_int)) + geom_boxplot()
p3 <- ggplot(cors, aes(method, cor)) + geom_boxplot() + ylim(-1, 1)

pdf("figs/slope-intercept-cv-re-ram-ensembles.pdf", width = 8, height = 9)
gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
dev.off()
