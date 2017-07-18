# compare the datalimited::prm function to the original fits
# start with setwd("analysis") from the ensembles root folder

library("datalimited")
library("dplyr")
library("ggplot2")

d <- dplyr::inner_join(datalimited::ram_ts, datalimited::spp_categories,
  by = "scientificname")
ram_prm_dat <- plyr::ddply(d, "stockid", function(x) {
  datalimited::format_prm(year = x$year, catch = x$catch,
    bbmsy = x$bbmsy_ram,
    species_cat = x$spp_category[1L])
})

m <- datalimited::fit_prm(ram_prm_dat)
d <- ram_prm_dat
d <- cbind(select(d, stockid), datalimited::predict_prm(d, ci = TRUE))
d <- dplyr::rename(d, package = bbmsy_q50)

ram <- readRDS("generated-data/ram_fits.rds") # original fits
qq <- subset(ram, method == "Costello")
q <- dplyr::inner_join(select(qq, stockid, tsyear, b_bmsy_est),
  select(d, stockid, year, package) %>% rename(tsyear = year))
qlong <- reshape2::melt(q, id.vars = c("stockid", "tsyear"))

set.seed(99)
ids <- sample(unique(qlong$stockid), size = 24)
ggplot(dplyr::filter(qlong, stockid %in% ids),
  aes(tsyear, value)) + facet_wrap(~stockid) +
  geom_line(aes(colour = variable))

ggplot(q, aes(b_bmsy_est, package)) + geom_point(alpha = 0.1) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1, colour = "red")

## CIs:
# ggplot(dplyr::filter(d, stockid %in% unique(d$stockid)[1:10]),
#   aes(year, package)) +
#   geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), fill = "red", alpha = 0.3) +
#   geom_line() +
#   ylab(expression(B/B[MSY])) +
#   facet_wrap(~stockid)
