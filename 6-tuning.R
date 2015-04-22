# try various models and tunning them using the caret package

# library(randomForest)

# make sure dplyr isn't loaded since caret relies on plyr

library(caret)

ram_fits <- readRDS("generated-data/ram_fits.rds")

x <- dplyr::select(ram_fits, year_before_end, COMSIR, Costello, SSCOM, CMSY_new_prior,
  resilience, sigma, slope, habitat)

library(doMC)
registerDoMC(4)
set.seed(1) # for cross-validation cuts
m <- caret::train(x = x, y = ram_fits$b2bmsy_true, metric = "RMSE",
  trControl = trainControl(method = "repeatedCV", number = 2, repeats = 2),
  method = "gbm",
  verbose = FALSE,
  tuneGrid = expand.grid(
    interaction.depth = c(2, 4, 6, 8, 10),
    n.trees = c(1000, 2500, 6000, 8000, 10000),
    shrinkage = c(0.1, 0.05, 0.01)))
ggplot(m, metric = "Rsquared")
ggplot(m, metric = "RMSE")
saveRDS(m, file = "generated-data/caret-gbm.rds")
# so, we'll pick interaction.depth = 10, n.trees = 2500, shrinkage = 0.05

# cforest
# gam
# lm
# rf
# # Random Forest (method = 'rf')
# # For classification and regression using package randomForest with tuning parameters:
# # Number of Randomly Selected Predictors (mtry, numeric)
# gbm
