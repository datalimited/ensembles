# Example time series plot to motivate the study

library(dplyr)

# try: SBWHITACIR
ram <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE) %>%
  filter(stockid == "SBWHITACIR") %>%
  select(tsyear, Bbmsy_toUse) %>%
  rename(yr = tsyear, b_bmsy = Bbmsy_toUse)

load("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/CMSY/cmsy_rlegacy_results_table_v0.RData")
cmsy <- cmsy.rlegacy.df0 %>%
  filter(stock_id == "SBWHITACIR") %>%
  select(yr, b_bmsyiq25, b_bmsy, b_bmsyiq75) %>%
  mutate(method = "CMSY")

load("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/COMSIR/comsir_rlegacy_results_table_v0.RData")
comsir <- comsir.rlegacy.df0 %>%
  filter(stock_id == "SBWHITACIR") %>%
  select(yr, b_bmsyiq25, b_bmsy, b_bmsyiq75) %>%
  mutate(method = "COM-SIR")

load("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/SSCOM/sscom_rlegacy_results_table_v0.RData")
sscom <- sscom.rlegacy.df0 %>%
  filter(stock_id == "SBWHITACIR") %>%
  select(yr, b_bmsy_iq25, b_bmsy, b_bmsy_iq75) %>%
  rename(b_bmsyiq25 = b_bmsy_iq25, b_bmsyiq75 = b_bmsy_iq75) %>%
  mutate(method = "SSCOM")

costello <- read.csv("~/Dropbox/FisheriesWorkingGroupPhaseII/RAM_Legacy_fits/Costello/Costello_rlegacy_results.csv")
cdat <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE)
costello <- merge(costello, unique(cdat[,c("stockid","stocklong")]), by="stocklong", all.x=TRUE)
costello <- costello %>% filter(stockid == "SBWHITACIR") %>%
  select(year, BvBmsy, BvBmsy_LogSD) %>%
  mutate(year, b_bmsyiq25 = exp(log(BvBmsy) - 0.5 * BvBmsy_LogSD),
    b_bmsy = BvBmsy,
    b_bmsyiq75 = exp(log(BvBmsy) + 0.5 * BvBmsy_LogSD)) %>%
  select(-BvBmsy, -BvBmsy_LogSD) %>%
  rename(yr = year) %>%
  mutate(method = "mPRM")

dl_dat <- bind_rows(cmsy, comsir, sscom, costello) %>%
  arrange(method, yr)

#library(ggplot2)
#ggplot(dl_dat, aes(yr, b_bmsy, colour = method)) + geom_line() + xlim(1990, 2011) +
#geom_ribbon(aes(ymax = b_bmsyiq75, ymin = b_bmsyiq25, fill = method), alpha = 0.2) +
#geom_line(data = ram, aes(yr, b_bmsy), colour = "black", fill = "black")

pal <- paste0(RColorBrewer::brewer.pal(4, "Dark2")) # colour-blind proof
pal <- c("#728905", "#2075C7", "#465A61", "#D01B24") # solarized
pal <- c("#BD1550", "#E97F02", "#F8CA00", "#8A9B0F") # http://www.colourlovers.com/palette/848743/(◕_”_◕)
pal <- c("#00A0B0", "#6A4A3C", "#CC333F", "#EDC951") #"EDC951" #http://www.colourlovers.com/palette/1473/Ocean_Five
pal <- paste0(RColorBrewer::brewer.pal(4, "Set2"))
library(RColorBrewer)
pals <- c("BuGn", "BuPu", "GnBu", "OrRd")
pals <- c("Reds", "Greens", "Blues", "Purples")
pal <- sapply(pals, function(x) brewer.pal(9, x)[7]) %>% as.character

col_df <- data_frame(col = pal,
  method = c("CMSY", "COM-SIR", "SSCOM", "mPRM"))
dl_dat$col <- NULL # for re-running
dl_dat <- inner_join(dl_dat, col_df)

plot_method <- function(dat) {
  polygon(c(dat$yr, rev(dat$yr)), c(dat$b_bmsyiq25, rev(dat$b_bmsyiq75)), border = NA,
    col = paste0(dat$col, 20))
  lines(dat$yr, dat$b_bmsy, col = dat$col, lwd = 2.5)
}

xlim <- c(1990, 2010)
ylim <- range(c(dl_dat$b_bmsyiq25, dl_dat$b_bmsyiq75))

pdf("../figs/motivate.pdf", width = 4.6, height = 2.9)
par(mar = c(3, 3.4, .5, 4.5), cex = 0.8, oma = c(0, 0, 0, 0), tck = -0.015,
  mgp = c(2, 0.5, 0), col.axis = "grey40", las = 1)
#par(family="serif")
plot(1, 1, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE)
abline(h = 1, lty = 2, col = "grey40", lwd = 2)
plyr::d_ply(dl_dat, "method", plot_method)
lines(ram$yr, ram$b_bmsy, lwd = 3.5, col = "grey30")
axis(2, col = "grey40", cex.axis = 0.9)
par(mgp = par()[["mgp"]] + c(0, -0.25, 0))
axis(1, col = "grey40", cex.axis = 0.9)
box(col = "grey40")
mtext("Year", side = 1, line = 1.7, col = "grey20", cex = 0.8)
mtext(expression(widehat(B/B[MSY])), side = 2, line = 1.7, col = "grey20", cex = 0.8,
  las = 0)

lab <- filter(dl_dat, yr == max(dl_dat$yr))
par(xpd = NA)
text(lab$yr+0, lab$b_bmsy, lab$method, pos = 4, col = lab$col, cex = 0.9)
text(unique(lab$yr)+0, 0.05 + ram[nrow(ram), "b_bmsy"][[1]], "Assessed", pos = 4, cex = 0.9)
text(xlim[1]-0.5, ylim[2]-0.1,
  paste0("Example stock:\n", unique(filter(cdat, stockid == "SBWHITACIR")$stocklong)),
  cex = 0.9, pos = 4, col = "grey20")

dev.off()
