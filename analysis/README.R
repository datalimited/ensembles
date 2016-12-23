install.packages(c("mblm", "randomForest", "dplyr", "plyr", "ggplot2",
  "reshape2", "gbm", "devtools", "doParallel", "hexbin", "pROC",
  "RColorBrewer", "Rcpp", "R2jags", "arm", "bbmle"))
devtools::install_github("datalimited/datalimited")

# start in the root "ensembles" folder, then:
setwd("analysis")
source("0-ensemble-functions.R")
source("1-read-old-ram-fits.R")
source("2-merge-ram-fits.R")
source("3-load-sim-data.R")
source("4-add-spectral.R")
source("5-simulation-ensembles.R")
source("6-sim-to-ram.R")
source("7-plot-sims.R")
source("8-main-figs.R")
source("8-motivate-fig.R")
source("9-output-values.R")

