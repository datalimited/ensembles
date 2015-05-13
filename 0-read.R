library("dplyr")
# Uncomment to re-cache the RAM data:
# dir.create("ram-data")
# setwd("ram-data")
# ramlegacy::make_ramlegacy() # devtools::install_github("seananderson/ramlegacy")
# setwd("../")
# ram <- src_sqlite("ram-data/ramlegacy.sqlite3")

cdat_old <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv",
  stringsAsFactors = FALSE) %>%
  select(stockid, res, tsyear, CtoUse, Bbmsy_toUse, scientificname) %>%
  filter(!is.na(CtoUse))
priors_old <- read.csv("RL_pred_Bbmsy.csv", stringsAsFactors = FALSE) %>%
  select(stockid, log.SD, log.Mean, post.median, stocklong) %>%
  rename(
    log_sd = log.SD,
    log_mean = log.Mean,
    post_median = post.median) %>%
  inner_join(cdat_old, by = "stockid") %>%
  mutate(res = tolower(res))
priors_old <- priors_old[!duplicated(priors_old), ]

clean_ramnames <- function(x) {
  names(x) <- tolower(names(x))
  names(x) <- gsub("\\.", "_", names(x))
  names(x) <- gsub("touse", "_touse", names(x))
  x
}

ramts <- priors_old
ramts <- clean_ramnames(ramts)
ramts <- rename(ramts, year = tsyear, bbmsy_ram = bbmsy__touse, catch = c_touse)

saveRDS(ramts, file = "generated-data/ramts.rds")
