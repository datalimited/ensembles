<!-- README.md is generated from README.Rmd. Please edit that file -->
Improving estimates of population status and trend with superensemble models
============================================================================

[![DOI](https://zenodo.org/badge/33008246.svg)](https://zenodo.org/badge/latestdoi/33008246)

This repository holds the analysis for: Anderson, S. C., A. B. Cooper, O. P. Jensen, C. Minto, J. T. Thorson, J. C. Walsh, M. Dickey-Collas, K. M. Kleisner, C. Longo, G. C. Osio, D. Ovando, I. Mosqueira, A. A. Rosenberg, and E. R. Selig. 2017. Improving estimates of population status and trend with superensemble models. Fish and Fisheries. In press.

To recreate the analysis, source the `.R` files in the `analysis` folder in order. Or run the following in R:

``` r
packages <- c("mblm", "randomForest", "plyr", "dplyr", "ggplot2",
  "reshape2", "gbm", "devtools", "doParallel", "hexbin", "pROC",
  "RColorBrewer", "Rcpp", "R2jags", "arm", "bbmle")
```

``` r
install.packages(packages)
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
```

Generated data ends up in `analysis/generated-data` and figures in `figs`. Intermediate versions of some raw data have been cached in the `analysis/generated-data` folder so that the scripts will run from a freshly cloned version of this repository. (Some of the raw data files are too large for Git.)

Package versions used:
----------------------

| package      | version     | date       |
|:-------------|:------------|:-----------|
| arm          | 1.9-3       | 2016-11-27 |
| bbmle        | 1.0.18      | 2016-02-11 |
| datalimited  | 0.1.0       | 2016-12-23 |
| devtools     | 1.12.0.9000 | 2016-12-18 |
| doParallel   | 1.0.10      | 2015-10-14 |
| dplyr        | 0.5.0       | 2016-06-24 |
| gbm          | 2.1.1       | 2015-03-11 |
| ggplot2      | 2.2.0       | 2016-11-11 |
| hexbin       | 1.27.1      | 2015-08-19 |
| mblm         | 0.12        | 2013-12-30 |
| plyr         | 1.8.4       | 2016-06-08 |
| pROC         | 1.8         | 2015-05-05 |
| R2jags       | 0.5-7       | 2015-08-23 |
| randomForest | 4.6-12      | 2015-10-07 |
| RColorBrewer | 1.1-2       | 2014-12-07 |
| Rcpp         | 0.12.8      | 2016-11-17 |
| reshape2     | 1.4.2       | 2016-10-22 |
