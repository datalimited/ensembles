library("dplyr")

lh <- readRDS("raw-data/OldRAMData-for-lifehistory-categories.RDS")
lh <- lh %>% filter(database == "Ram") %>%
  select(spname, habitat, pricecategory, envtemp, tl,
  agematmax, resilience, LifeHist2)
lh <- lh[!duplicated(lh), ]

ram <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE)

ram_fits <- readRDS("generated-data/ram-orig-fits.rds")

cmsy_2015 <- readRDS("generated-data/ram-cmsy-2015.rds")

ram_fits <- ram_fits %>% rename(stockid = stock)

ram_fits <- filter(ram_fits, !method %in% c("CMSY_old_prior", "CMSY_new_prior"))

cmsy_2015 <- cmsy_2015 %>% select(stockid, year, bbmsy_q50) %>%
  mutate(method = "CMSY_new_prior") %>%
  rename(b2bmsy = bbmsy_q50)

ram_fits <- rbind(ram_fits, cmsy_2015)

ram_fits <- rename(ram_fits, tsyear = year)

d <- ram_fits %>% inner_join(select(ram, tsyear, stockid, stocklong,
  Bbmsy_toUse, res, CtoUse))

ram_sp <- ram[,c("stockid", "scientificname")]
ram_sp <- ram_sp[!duplicated(ram_sp), ]

# some duplicates on join - no idea why, checking later
d <- d[!duplicated(d), ]

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

ram_fits <- d2_long
ram_fits <- ram_fits %>% rename(b2bmsy_true = Bbmsy_toUse)

# Merge in exploitation:
ram <- src_sqlite("ram-data/ramlegacy.sqlite3")
ramts <- tbl(ram, "timeseries_values_views") %>% as.data.frame %>%
  select(stockid, year, Utouse, Ctouse) %>%
  rename(tsyear = year, effort = Utouse, catch = Ctouse)

ram_fits <- inner_join(ram_fits, ramts)
ram_fits <- inner_join(ram_fits, ram_sp)
ram_fits <- inner_join(ram_fits, rename(lh, scientificname = spname))

ram_fits <- mutate(ram_fits,
  habitat = as.factor(habitat),
  pricecategory = as.factor(pricecategory),
  envtemp = as.factor(envtemp),
  resilience = as.factor(resilience))

# Categorize the exploitation: increasing, decreasing

categorize_effort <- function(effort) {
  d <- data.frame(year = seq_along(effort), log(effort))
  m <- lm(effort ~ year, data = d)
  slope <- summary(m)$coefficient["year","Estimate"]
  p_slope <- summary(m)$coefficient["year","Pr(>|t|)"]
  sigma <- summary(m)$sigma
  if (p_slope < 0.1 & slope >= 0)
    s <- "increasing"
  if (p_slope < 0.1 & slope < 0)
    s <- "decreasing"
  if (p_slope >= 0.1)
    s <- "unknown"
  # sigma = 0.1 is about the median in the RAM dataset
  sigma <- ifelse(sigma > 0.1, "high", "low")
  data.frame(sigma = sigma, slope = s)
}

qq <- plyr::ddply(ram_fits, c("method", "stockid"), function(x) {
  categorize_effort(x$effort)
})
ram_fits$sigma <- NULL
ram_fits$slope <- NULL
ram_fits <- inner_join(ram_fits, qq)

# add max catch:
ram_fits <- ram_fits %>%
  group_by(stockid) %>%
  mutate(log_max_catch = log(max(catch)),
    log_total_catch = log(sum(catch))) %>%
  as.data.frame

# strip out all but the last 10 years, which we'll focus on:
# we'll also drop all with less than 15 years
# that removes ballpark 25 stocks
# first, make sure in order:
ram_fits <- arrange(ram_fits, method, stockid, tsyear)
ram_fits <- plyr::ddply(ram_fits, c("method", "stockid"), function(x) {
  if (nrow(x) >= 15) {
    out <- x[(nrow(x)-9):nrow(x),]
    out$year_before_end <- seq_len(10)
    out
  }
})
ram_fits <- select(ram_fits, -agematmax, -LifeHist2)
ram_fits[ram_fits$tl < 0, "tl"] <- NA # were some -999 values

# Now widen the data for model fitting etc:
zz <- reshape2::dcast(ram_fits, tsyear + stockid ~ method, value.var = "b2bmsy")
zz <- select(zz, -Minto)

ram_fits <- select(ram_fits, -method, -b2bmsy)
ram_fits <- ram_fits[!duplicated(ram_fits), ]
ram_fits <- inner_join(ram_fits, zz)

saveRDS(ram_fits, "generated-data/ram_fits.rds")
