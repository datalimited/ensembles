# Check ROC curves for over/under B/BMSY = 1

get_roc <- function(true, est) {
  library("pROC") # required or pROC:: will generate errors
  y <- pROC::roc(response = true, predictor = est, 5)
  data.frame(sens = y$sensitivities, spec = y$specificities,
    auc = as.numeric(y$auc))
}

# cv_sim_long <- readRDS("generated-data/cv_sim_long.rds") %>% as.data.frame
# cv_sim_long <- cv_sim_long %>% mutate(
#   above1_est = ifelse(bbmsy_est > 1, 1, 0),
#   above1_true = ifelse(bbmsy_true > 1, 1, 0))
#
# rocs <- filter(cv_sim_long, type == "mean") %>% group_by(method_id) %>%
#   do({get_roc(true = .$above1_true, est = .$above1_est)})
# ggplot(rocs, aes(spec, sens, colour = method_id)) + geom_line()

cv_sim_binary <- readRDS("generated-data/cv_sim_binary.rds") %>% as.data.frame

# cv_sim_binary <- cv_sim_binary %>% mutate(
#   CMSY = ifelse(CMSY > 1, 1, 0),
#   COM.SIR = ifelse(COM.SIR > 1, 1, 0),
#   Costello = ifelse(Costello > 1, 1, 0),
#   SSCOM = ifelse(SSCOM > 1, 1, 0),
#   mean_ensemble = ifelse(mean_ensemble > 1, 1, 0))

# switch to long format for ROCing:
cv_sim_binary_long <- select(cv_sim_binary, -bbmsy_true_mean) %>%
 reshape2::melt(id.vars = c("stock_id", "iter", "test_iter", "LH", "above_bbmsy1_true"),
  variable.name = "method_id", value.name = "bbmsy_est")

rocs <- cv_sim_binary_long %>% group_by(method_id) %>%
  do({get_roc(true = .$above_bbmsy1_true, est = .$bbmsy_est)})

p <- ggplot(rocs, aes(spec, sens, colour = method_id)) + geom_line() + coord_equal() +
  geom_abline(intercept = 1, slope = 1, lty = 2, col = "darkgrey") +
  xlim(1, 0) +
  xlab("Specificity") + ylab("Sensitivity") + theme_bw() +
  theme(plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
print(p)
