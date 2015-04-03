# start with the existing data based on the simulations

library("dplyr")
lh <- readRDS("raw-data/OldRAMData-for-lifehistory-categories.RDS")
lh <- select(lh, spname, habitat, pricecategory, envtemp, tl,
  agematmax, resilience, LifeHist2)
  # filter(database == "Ram")
lh <- lh[!duplicated(lh), ]

ram <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE)

ram_fits <- readRDS("generated-data/ram-orig-fits.rds")

d <- ram_fits %>% rename(stockid = stock, tsyear = year) %>%
  inner_join(select(ram, tsyear, stockid, stocklong, Bbmsy_toUse, res, CtoUse))

ram_sp <- ram[,c("stockid", "scientificname")]
ram_sp <- ram_sp[!duplicated(ram_sp), ]

# some duplicates on join - no idea why, checking later
d <- d[!duplicated(d), ]

library(ggplot2)
# simdat <- load("rat-data/CMSY_COMSIR_SSCOM_COSTELLO_STO_BATCH_ALL_RESULTS 2013-06-24.Rdata")
# batch1results

# need to widen and lengthen to strip out years where anything is NA
assessed_bbmsy <- select(d, stockid, tsyear, Bbmsy_toUse)
assessed_bbmsy <- assessed_bbmsy[!duplicated(assessed_bbmsy), ]
d2 <- reshape2::dcast(d, stockid + tsyear + stocklong ~ method, value.var = "b2bmsy") %>%
  inner_join(assessed_bbmsy)
d2 <- na.omit(d2)
set.seed(123)
d2$Minto <- rnorm(nrow(d2), mean = 1.0, sd = 0.3)

d2_long <- reshape2::melt(select(d2, -Bbmsy_toUse), id.vars = c("stockid", "stocklong", "tsyear"))
d2_long <- inner_join(d2_long, d2[,c("stockid", "tsyear", "Bbmsy_toUse")]) %>%
  rename(b2bmsy = value, method = variable)

# Get median b/bmsy in last 5 years
d5 <- d2_long[!is.na(d2_long$b2bmsy) & !is.na(d2_long$Bbmsy_toUse), ] %>%
  group_by(stockid, method) %>%
  summarise(
    maxyear = max(tsyear),
    nyear = length(b2bmsy),
    b2bmsy = median(b2bmsy[(nyear-4):nyear]),
    Bbmsy_toUse = median(Bbmsy_toUse[(nyear-4):nyear])) %>%
  as.data.frame

mare_ <- function(xhat, xtrue) {
  median(abs((xhat - xtrue) / xtrue))
}

d5_mare <- d5 %>% group_by(stockid, method) %>%
  mutate(mare = mare_(b2bmsy, Bbmsy_toUse))

d5_mean <- d5 %>% group_by(stockid) %>%
  summarise(b2bmsy = mean(b2bmsy), Bbmsy_toUse = Bbmsy_toUse[1])

d5_mean <- as.data.frame(d5_mean) %>%
  group_by(stockid) %>%
  mutate(mare = mare_(b2bmsy, Bbmsy_toUse)) %>%
  mutate(method = "ensemble-average") %>%
  as.data.frame

cols <- c("stockid", "method", "b2bmsy", "Bbmsy_toUse", "mare")
d5_mare <- rbind(d5_mare[,cols], d5_mean[,cols])

# make a wide dataset for running models:
d5_wide <- reshape2::dcast(d5, stockid ~ method, value.var = "b2bmsy") %>%
  inner_join(d5_mean[,c("stockid", "Bbmsy_toUse")])

# bring in some life-history data etc.
d5_wide <- inner_join(d5_wide, ram_sp) %>% as.data.frame
d5_wide <- lh %>% rename(scientificname = spname) %>%
  inner_join(d5_wide)
d5_wide <- mutate(d5_wide, habitat = as.factor(habitat), pricecategory = as.factor(pricecategory), envtemp = as.factor(envtemp), resilience = as.factor(resilience))

f_lm <- as.formula("Bbmsy_toUse ~ CMSY_new_prior*habitat + CMSY_old_prior*habitat + COMSIR*habitat  + SSCOM*habitat")
f_gbm <- as.formula("Bbmsy_toUse ~ CMSY_new_prior + CMSY_old_prior + COMSIR + SSCOM + Minto + habitat")
# run some in-bag models:
m_lm <- lm(f_lm, data = d5_wide)
m_gbm <- gbm::gbm(f_gbm,
  data = d5_wide, n.trees = 2000, interaction.depth = 5, shrinkage = 0.001,
  distribution = "gaussian")

# first bind lm:
d5_wide$b2bmsy <- predict(m_lm)
d5_wide <- group_by(d5_wide, stockid) %>%
  mutate(mare = mare_(b2bmsy, Bbmsy_toUse)) %>%
  as.data.frame %>%
  mutate(method = "lm-in-bag")
d5_mare <- rbind(d5_mare[,cols], d5_wide[,cols])

# then bind gbm (this is terrible, I know)
d5_wide$mare <- NULL
d5_wide$b2bmsy <- NULL
d5_wide$b2bmsy <- gbm::predict.gbm(m_gbm, n.trees = 2000, type = "response")
d5_wide <- group_by(d5_wide, stockid) %>%
  mutate(mare = mare_(b2bmsy, Bbmsy_toUse)) %>%
  as.data.frame %>%
  mutate(method = "gbm-in-bag")
d5_mare <- rbind(d5_mare[,cols], d5_wide[,cols])
d5_mare <- as.data.frame(d5_mare)

# and first pass at example out-of-bag:
set.seed(1234)

cross_val_weights <- function() {
  train_ids <- base::sample(seq_len(nrow(d5_wide)), round(nrow(d5_wide)*0.5))
  test_ids <- seq_len(nrow(d5_wide))[-train_ids]

  m_lm <- lm(f_lm,
    data = d5_wide[train_ids, ])
  m_gbm <- gbm::gbm(f_gbm,
    data = d5_wide[train_ids, ], n.trees = 1000, interaction.depth = 5, shrinkage = 0.001,
    distribution = "gaussian")

  lm_out <- vector(mode = "numeric", length = nrow(d5_wide))
  gbm_out <- vector(mode = "numeric", length = nrow(d5_wide))

  tryCatch({lm_out[test_ids] <- predict(m_lm, newdata = d5_wide[test_ids,])},
    error = function(e) rep(NA, length(test_ids)))
  lm_out[train_ids] <- NA

  tryCatch({gbm_out[test_ids] <- gbm::predict.gbm(m_gbm, n.trees = 1000,
    type = "response", newdata = d5_wide[test_ids, ])},
    error = function(e) rep(NA, length(test_ids)))
  gbm_out[train_ids] <- NA

  cbind(lm_out, gbm_out)
}

x <- lapply(seq_len(300L), function(yy) cross_val_weights())
cross_out_lm <- matrix(ncol = length(x), nrow = nrow(d5_wide))
cross_out_gbm <- matrix(ncol = length(x), nrow = nrow(d5_wide))
for (i in seq_along(x)) {
  cross_out_lm[, i] <- x[[i]][,1]
  cross_out_gbm[, i] <- x[[i]][,2]
}

# now build back into relative error df:
d5_wide$b2bmsy <- apply(cross_out_lm, 1, median, na.rm = TRUE)
d5_wide <- group_by(d5_wide, stockid) %>%
  mutate(mare = mare_(b2bmsy, Bbmsy_toUse)) %>%
  mutate(method = "lm-cv")
d5_mare <- rbind(d5_mare[,cols], d5_wide[,cols])

d5_wide$b2bmsy <- apply(cross_out_gbm, 1, median, na.rm = TRUE)
d5_wide <- group_by(d5_wide, stockid) %>%
  mutate(mare = mare_(b2bmsy, Bbmsy_toUse)) %>%
  mutate(method = "gbm-cv")
d5_mare <- rbind(d5_mare[,cols], d5_wide[,cols])

## order and plot:
method_order <- d5_mare %>% group_by(method) %>%
  summarise(median_are = median(mare, na.rm = TRUE)) %>%
  arrange(median_are) %>%
  mutate(or = 1:length(method), method_ordered = reorder(method, or)) %>%
  as.data.frame
d5_mare <- inner_join(method_order, d5_mare)

d5_mare <- filter(d5_mare, !method %in% c("Costello", "gbm-in-bag", "lm-in-bag"))

scatter_plot <- function(dat) {
  ggplot(dat, aes(b2bmsy, Bbmsy_toUse, colour = method)) +
    geom_point(alpha = 0.2) +
    theme_bw() +
    ylim(0, 2) +
    xlim(0, 2) + stat_smooth(method = "loess", alpha = 0.15, lwd = 1) +
    geom_abline(aes(intercept = 0, slope = 1), lty = 3) +
    xlab("B/Bmsy estimated") +
    ylab("B/Bmsy assessment") +
    coord_fixed()
}

p <- scatter_plot(filter(d5_mare, !method %in%
    c("CMSY_new_prior", "CMSY_old_prior", "COMSIR", "Costello", "SSCOM")))
ggsave("figs/ensemble-ram-scatter.pdf", width = 8, height = 5)

p <- scatter_plot(filter(d5_mare, method %in%
    c("CMSY_new_prior", "CMSY_old_prior", "COMSIR", "Costello", "SSCOM")))
ggsave("figs/single-ram-scatter.pdf", width = 8, height = 5)

p <- ggplot(d5_mare, aes(method_ordered, log10(mare))) +
  geom_boxplot(notch = TRUE, alpha = 0.7, outlier.colour = NA) +
  geom_point(alpha = 0.1, position = position_jitter(width = 0.1)) +
  ylab("log10(Absolute relative error)") +
  xlab("Method from highest to lowest median absolute relative error") + coord_flip()
ggsave("figs/ensemble-re-ram-boxplots.pdf", width = 7, height = 6)

sp <- group_by(d5_mare, method) %>%
  summarise(sp = cor(b2bmsy, Bbmsy_toUse, method = "spearman")) %>%
  arrange(sp) %>%
  mutate(or = seq_along(method), method = reorder(method, -or))

p <- ggplot(sp, aes(sp, method)) + geom_point() + xlab("Spearman's correlation")
ggsave("figs/ensemble-spearman.pdf", width = 5, height = 6)

pdf("figs/gbm-partial.pdf", width = 8, height = 5)
par(cex = 0.6)
par(mfrow = c(2, 3));for(i in 1:5) plot(m_gbm, i.var = i, ylim = c(0.8, 1.4))
dev.off()

pdf("figs/lm-coef.pdf", width = 8, height = 8)
par(cex = 0.8, mar = c(1, 15, 3, 1))
arm::coefplot(m_lm)
dev.off()
