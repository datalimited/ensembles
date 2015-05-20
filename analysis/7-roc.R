# Check ROC curves for over/under B/BMSY = 1

get_roc <- function(true, est) {
  library("pROC") # required or pROC:: will generate errors
  y <- pROC::roc(response = true, predictor = est, 5)
  data.frame(sens = y$sensitivities, spec = y$specificities,
    auc = as.numeric(y$auc))
}

cv_sim_binary <- readRDS("generated-data/cv_sim_binary.rds") %>% as.data.frame
cv_ram_binary <- readRDS("generated-data/cv_ram_binary.rds") %>% as.data.frame

# If you want to be picky and take the individual models
# at their word and use a threshold of 1.
# Otherwise, they get treated as continuous variables
# with a moving threshold by the ROC function.

# cv_sim_binary <- cv_sim_binary %>% mutate(
#   CMSY = ifelse(CMSY > 1, 1, 0),
#   COM.SIR = ifelse(COM.SIR > 1, 1, 0),
#   Costello = ifelse(Costello > 1, 1, 0),
#   SSCOM = ifelse(SSCOM > 1, 1, 0),
#   mean_ensemble = ifelse(mean_ensemble > 1, 1, 0))

# switch to long format for ROCing:
cv_sim_binary_long <- select(cv_sim_binary, -bbmsy_true_mean) %>%
 reshape2::melt(id.vars = c("stockid", "iter", "test_iter", "LH", "above_bbmsy1_true"),
  variable.name = "method", value.name = "bbmsy_est")

cv_ram_binary_long <- select(cv_sim_binary, -bbmsy_true_mean) %>%
 reshape2::melt(id.vars = c("stockid", "test_iter", "habitat", "above_bbmsy1_true"),
  variable.name = "method", value.name = "bbmsy_est")

rocs_sim <- cv_sim_binary_long %>% group_by(method) %>%
  do({get_roc(true = .$above_bbmsy1_true, est = .$bbmsy_est)})

rocs_ram <- cv_ram_binary_long %>% group_by(method) %>%
  do({get_roc(true = .$above_bbmsy1_true, est = .$bbmsy_est)})

p <- ggplot(rocs_sim, aes(spec, sens, colour = method)) + geom_line() + coord_equal() +
  geom_abline(intercept = 1, slope = 1, lty = 2, col = "darkgrey") +
  xlim(1, 0) +
  xlab("Specificity") + ylab("Sensitivity") + theme_bw() +
  theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave("../figs/roc-sim.pdf", width = 7, height = 5)

p <- ggplot(rocs_ram, aes(spec, sens, colour = method)) + geom_line() + coord_equal() +
  geom_abline(intercept = 1, slope = 1, lty = 2, col = "darkgrey") +
  xlim(1, 0) +
  xlab("Specificity") + ylab("Sensitivity") + theme_bw() +
  theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
ggsave("../figs/roc-ram.pdf", width = 7, height = 5)
