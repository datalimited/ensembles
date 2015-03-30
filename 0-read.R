library("dplyr")
# Uncomment to re-cache the RAM data:
dir.create("ram-data")
setwd("ram-data")
ramlegacy::make_ramlegacy() # devtools::install_github("seananderson/ramlegacy")
setwd("../")
ram <- src_sqlite("ram-data/ramlegacy.sqlite3")

cdat_old <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors = FALSE) %>%
  select(stockid, res)
priors_old <- read.csv("RL_pred_Bbmsy.csv", stringsAsFactors = FALSE) %>%
  select(stockid, log.SD, log.Mean, post.median) %>%
  rename(
    log_sd = log.SD,
    log_mean = log.Mean,
    post_median = post.median) %>%
  inner_join(cdat_old, by = "stockid") %>%
  mutate(res = tolower(res))
priors_old <- priors_old[!duplicated(priors_old), ]

ramts <- tbl(ram, "timeseries_values_views") %>% as.data.frame %>%
  left_join(priors_old, by = "stockid") %>%
  select(-assessid)

clean_ramnames <- function(x) {
  names(x) <- tolower(names(x))
  names(x) <- gsub("\\.", "_", names(x))
  names(x) <- gsub("touse", "_touse", names(x))
  x
}

ramts <- clean_ramnames(ramts) %>%
  select(-b_bmsy, -ssb_ssbmsy, -f_fmsy, - tc, -tl, -ssb, -tb) %>%
  arrange(stockid, year)

# Don't have priors or weren't in original data:
dropped <- ramts %>% filter(is.na(post_median)) %>%
  select(stockid) %>%
  filter(!duplicated(stockid))

if(length(dropped[[1]] > 0))
  warning(paste(length(dropped[[1]]), "stocks were dropped:",
    "they were not in the original prior data."))

saveRDS(ramts, file = "generated-data/ramts.rds")
