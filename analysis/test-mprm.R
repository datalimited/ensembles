# compare the datalimited::prm function to the original fits
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
d <- cbind(d, datalimited::predict_prm(d, ci = TRUE))
d <- dplyr::rename(d, package = fit)

ram <- readRDS("generated-data/ram_fits.rds")
qq <- subset(ram, method == "Costello")
q <- dplyr::inner_join(select(qq, stockid, tsyear, b_bmsy_est),
  select(d, stockid, year, package) %>% rename(tsyear = year))
qlong <- reshape2::melt(q, id.vars = c("stockid", "tsyear"))

ggplot(dplyr::filter(qlong, stockid %in% unique(qlong$stockid)[1:25]),
  aes(tsyear, value)) + facet_wrap(~stockid) +
  geom_line(aes(colour = variable))

# CIs:
ggplot(dplyr::filter(d, stockid %in% unique(d$stockid)[1:10]),
  aes(year, package)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#FF000040") +
  geom_line() +
  ylab(expression(B/B[MSY])) +
  facet_wrap(~stockid)

