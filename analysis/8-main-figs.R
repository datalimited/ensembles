# Make the main figures
library(dplyr)
library(ggplot2)
source("0-ensemble-functions.R")

# bring in the simulated data and shape it for plotting:
d_sim <- readRDS("generated-data/cv_sim_long.rds")

clean_names <- dplyr::data_frame(
  method = c("CMSY", "COMSIR", "Costello", "SSCOM",
    "gbm_ensemble", "rf_ensemble", "lm_ensemble", "mean_ensemble"),
  clean_method = c("CMSY", "COM-SIR", "mPRM", "SSCOM",
    "GBM Ensemble", "RF Ensemble", "LM Ensemble", "Mean Ensemble"),
  order = c(1, 2, 4, 3, 8, 7, 6, 5),
  label_fudge_x = c(
    0, -0.03,0,-0.04,
    0, 0,0,-0.05),
  label_fudge_y = c(
    0, -0.02,0.03,0.04,
    0, 0.02,0,-0.02))

d_sim <- suppressWarnings(inner_join(d_sim, clean_names))

d_sim_perf <- d_sim %>% group_by(test_iter, type, clean_method, order, label_fudge_x, label_fudge_y) %>%
  summarise(
    mare = median(abs(re)),
    mre = median(re),
    corr = cor(bbmsy_true_trans, bbmsy_est_trans, method = "spearman",
      use = "pairwise.complete.obs"))

d_sim_perf$clean_method <- as.factor(d_sim_perf$clean_method)
d_sim_perf$clean_method <- reorder(d_sim_perf$clean_method, d_sim_perf$order)

d_sim_perf_slope <- filter(d_sim_perf, type == "slope")
d_sim_perf <- filter(d_sim_perf, type == "mean")

d_sim_perf_long <- reshape2::melt(d_sim_perf,
  id.vars = c("test_iter", "clean_method", "type", "order", "label_fudge_x", "label_fudge_y"))
d_sim_perf_summ <- d_sim_perf_long %>%
  group_by(clean_method, order, variable, label_fudge_x, label_fudge_y) %>%
  summarise(
    m = median(value),
    l = quantile(value, 0.25),
    u = quantile(value, 0.75)) %>%
  as.data.frame
d_sim_perf_summ$text_col <- ifelse(grepl("Ensemble", d_sim_perf_summ$clean_method),
  "grey20", "grey50")
d_sim_perf_wide <- reshape2::dcast(d_sim_perf_summ, clean_method + order ~ variable,
  value.var = "m")

# Format basic simulation ensemble without covariates for plotting:
d_sim_basic <- readRDS("generated-data/cv_sim_mean_basic_long.rds")
d_sim_basic$bbmsy_est[d_sim_basic$bbmsy_est > 10] <- NA
d_sim_basic <- na.omit(d_sim_basic)
d_sim_basic <- suppressWarnings(inner_join(d_sim_basic, clean_names))

# Format the RAM ensemble data for plotting:
d_ram <- readRDS("generated-data/ram-ensemble-predicted.rds")
clean_names_ram <- clean_names
d_ram <- suppressWarnings(inner_join(d_ram, clean_names_ram))
d_ram$bbmsy_est <- as.numeric(as.character(d_ram$bbmsy_est))
d_ram <- filter(d_ram, bbmsy_true < 4, bbmsy_est < 4)

re_ram <- d_ram %>% mutate(
  sq_er = (bbmsy_est - bbmsy_true)^2,
  re    = (bbmsy_est - bbmsy_true) / bbmsy_true)

re_ram_sum <- re_ram %>% group_by(clean_method) %>%
  summarise(mare = median(abs(re)),
    msqe = mean(sq_er),
    mre = median(re),
    corr = cor(bbmsy_true, bbmsy_est, method = "spearman")) %>%
  as.data.frame()

# now make these data match the simulation data in format:
re_ram_sum_long <- reshape2::melt(re_ram_sum, id.vars = "clean_method") %>%
  inner_join(clean_names) %>%
  rename(m = value) %>%
  mutate(l = m, u = m) %>%
  mutate(text_col = ifelse(grepl("Ensemble", clean_method), "grey20", "grey50"))
re_ram_sum_long$clean_method <- reorder(re_ram_sum_long$clean_method, re_ram_sum_long$order)
re_ram_sum_long$label_fudge_x <- 0
re_ram_sum_long$label_fudge_y <- 0



get_roc <- function(true, est) {
  library("pROC") # required or pROC:: will generate errors
  y <- pROC::roc(response = true, predictor = est, 5)
  data.frame(sens = y$sensitivities, spec = y$specificities,
    auc = as.numeric(y$auc))
}

# get ROC and AUC for models applied to simulation data with b/bmsy as the response:
d_sim_mean <- d_sim %>% filter(type == "mean") %>%
  mutate(above1 = ifelse(bbmsy_true > 1, 1, 0))
rocs_sim_mean <- d_sim_mean %>% group_by(clean_method, order) %>%
  do({get_roc(true = .$above1, est = .$bbmsy_est)}) %>%
  as.data.frame
rocs_sim_mean$Ensemble <- ifelse(grepl("Ensemble", rocs_sim_mean$clean_method),
  "(A) Ensemble", "(B) Individual")
auc_sim_mean <- rocs_sim_mean %>%
  group_by(clean_method) %>%
  summarise(auc = auc[1])
rocs_sim_mean <- rocs_sim_mean %>% as.data.frame() %>%
  mutate(clean_method = paste(order, clean_method, sep = "-"))

# and do the same for the RAM data:
d_ram_mean <- d_ram %>%
  mutate(above1 = ifelse(bbmsy_true > 1, 1, 0))
rocs_ram_mean <- d_ram_mean %>% group_by(clean_method, order) %>%
  do({get_roc(true = .$above1, est = .$bbmsy_est)}) %>%
  as.data.frame
rocs_ram_mean$Ensemble <- ifelse(grepl("Ensemble", rocs_ram_mean$clean_method),
  "(A) Ensemble", "(B) Individual")
auc_ram_mean <- rocs_ram_mean %>%
  group_by(clean_method) %>%
  summarise(auc = auc[1])
rocs_ram_mean <- rocs_ram_mean %>% as.data.frame() %>%
  mutate(clean_method = paste(order, clean_method, sep = "-"))

make_roc_curves <- function(dat) {
  p <- ggplot(dat, aes(spec, sens, colour = clean_method, group = clean_method)) +
    geom_line() + coord_equal() +
    geom_abline(intercept = 1, slope = 1, lty = 2, col = "darkgrey") +
    xlim(1, 0) +
    xlab("Specificity") + ylab("Sensitivity") + theme_bw() +
    theme(plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) + facet_wrap(~Ensemble) +
    scale_colour_brewer(palette="Spectral", guide = guide_legend(title = "Model"))
    p
}

## Figure making:

## ROC curves:
# make_roc_curves(rocs_sim_all)
p <- make_roc_curves(rocs_sim_mean)
ggsave("../figs/roc-sim.pdf", width = 8, height = 5)

p <- make_roc_curves(rocs_ram_mean)
ggsave("../figs/roc-ram.pdf", width = 8, height = 5)

## Hexagon figures:
d_mean_plot <- filter(d_sim, type == "mean")
d_mean_plot$bbmsy_est[d_mean_plot$bbmsy_est > 10] <- NA
d_mean_plot <- na.omit(d_mean_plot)

d_mean_sim_and_ram <- bind_rows(
  mutate(d_mean_plot, data = "Simulated"),
  filter(d_ram, grepl("Ensemble", clean_method)) %>%
    mutate(order = order + 4, data = "RAM Legacy"))

ram_third_row_lims <- filter(d_ram, grepl("Ensemble", clean_method)) %>%
  summarise(lower = min(c(bbmsy_est, bbmsy_true)), upper = max(c(bbmsy_est, bbmsy_true))) %>%
  as.numeric

# Basic hexagon plot for simulation mean bbmsy:
pdf("../figs/hex-mean-sim-cv.pdf", width = 7, height = 3.9)
plot_hex_fig(d_mean_plot, xbins = 100L)
dev.off()

# Same as above but add the RAM ensembles as a third row:
pdf("../figs/hex-mean-sim-ram-cv.pdf", width = 6.0, height = 3.9)
plot_hex_fig(d_mean_sim_and_ram, xbins = 100L, xbins3 = 30L, lims_hex3 = ram_third_row_lims,
  count_transform3 = 10)
dev.off()

d_slope_plot <- filter(d_sim, type == "slope")
d_slope_plot$bbmsy_est[d_slope_plot$bbmsy_est > 10] <- NA
d_slope_plot <- na.omit(d_slope_plot)
pdf("../figs/hex-slope-sim.pdf", width = 7, height = 3.9)
plot_hex_fig(d_slope_plot, xbins = 80L,
  lims_hex = range(c(d_slope_plot$bbmsy_est, d_slope_plot$bbmsy_true)),
  xlim_plot = c(-0.5, 0.5), ylim_plot = c(-0.5, 0.5), axis_ticks = c(-0.4, 0, 0.4),
  xlab = expression(B/B[BMSY]~slope), ylab = expression(widehat(B/B[BMSY])~slope))
dev.off()

pdf("../figs/hex-mean-ram-cv.pdf", width = 8, height = 4)
plot_hex_fig(d_ram, add_hex = TRUE, alpha = 80,
  lims_hex = range(c(d_ram$bbmsy_est, d_ram$bbmsy_true)), xbins = 25L)
dev.off()

pdf("../figs/hex-mean-sim-basic-cv.pdf", width = 8, height = 4)
plot_hex_fig(d_sim_basic, xbins = 100L)
dev.off()

## Performance panel figures:
pdf("../figs/performance.pdf", width = 6, height = 3.1)
par(mfrow = c(1, 2), mgp = c(1.5, 0.4, 0), las = 1, tck = -0.012,
  oma = c(2.7, 3.5, .5, .5), cex = 0.8, mar = c(0, 0, 0, 0),
  xaxs = "i", yaxs = "i", col.axis = "grey50", col.lab = "grey50")
par(xpd = FALSE)
performance_panel(d_sim_perf_summ, xlim = c(0.2, 0.6), ylim = c(0.1, 0.7))
performance_panel(re_ram_sum_long, xlim =  c(0.48, 0.65), ylim = c(0.15, 0.45))
mtext("Within population inaccuracy (MAPE)", side = 1, line = 1.4, col = "grey40", cex = 0.8, outer = TRUE)
mtext("Across population correlation", side = 2, line = 2.2, col = "grey40", las = 0, cex = 0.8, outer = TRUE)
#text(0.54, 0.5, "Bias (MPE)", col = "grey30", pos = 4)
dev.off()


performance_ggplot <- function(dat) {
  p <- dat %>% ggplot(aes(mare, corr)) + geom_point(aes(colour = mre), size = 5) +
  geom_text(aes(label = clean_method), size = 2.7, vjust = 2) +
  scale_colour_gradient2(guide = guide_legend(title = "Bias\n(MRE)")) +
  theme_bw() + xlab("Inaccuracy (MARE)") + ylab("Rank-order correlation") +
  xlim(range(dat$mare) + c(-0.03, 0.01)) +
  ylim(range(dat$corr) + c(-0.05, 0.00)) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank())
  p
}
p1 <- performance_ggplot(d_sim_perf_wide)
p2 <- performance_ggplot(re_ram_sum)
pdf("../figs/performance-gg.pdf", width = 10, height = 4)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


