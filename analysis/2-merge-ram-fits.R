# This script does a bunch of messy merging of RAM data It takes the original
# fits, merges in the revised CMSY fits, renames
# some columns, and shapes the data appropriately for the ensemble models.

library("dplyr")

ram <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE)
ram_fits <- readRDS("generated-data/ram-orig-fits.rds")
load("raw-data/cmsy_rlegacy_results_table_v0.RData")
cmsy_refit <- cmsy.rlegacy.df0
cmsy_refit <- cmsy_refit %>%
  rename(stockid = stock_id, year = yr, b2bmsy = b_bmsy) %>%
  select(stockid, year, b2bmsy) %>%
  mutate(method = "CMSY")

ram_fits <- ram_fits %>% rename(stockid = stock)
ram_fits <- filter(ram_fits, !method %in% c("CMSY_old_prior", "CMSY_new_prior"))
ram_fits <- dplyr::bind_rows(ram_fits, cmsy_refit)
ram_fits <- rename(ram_fits, tsyear = year)

d <- ram_fits %>% inner_join(select(ram, tsyear, stockid, stocklong,
  Bbmsy_toUse, res, CtoUse))

ram_sp <- ram[,c("stockid", "scientificname")]
ram_sp <- ram_sp[!duplicated(ram_sp), ]

# some duplicates on join:
d <- d[!duplicated(d), ]

# need to widen and lengthen to strip out years where anything is NA
assessed_bbmsy <- select(d, stockid, tsyear, Bbmsy_toUse)
assessed_bbmsy <- assessed_bbmsy[!duplicated(assessed_bbmsy), ]
# how is convergence affecting lost stocks?
# all_na <- function(x) ifelse(sum(is.na(x)) == length(x), TRUE, FALSE)
# check_nas <- d2 %>% group_by(stockid)%>% summarise(COMSIR = all_na(COMSIR),
#   CMSY = all_na(CMSY), Costello = all_na(Costello), SSCOM = all_na(SSCOM))
# colSums(check_nas[,-1])
# COMSIR     CMSY Costello    SSCOM
# 80           59        0       44
d2 <- reshape2::dcast(d, stockid + tsyear + CtoUse + stocklong ~ method,
  value.var = "b2bmsy") %>%
  inner_join(assessed_bbmsy)
d2 <- na.omit(d2) # losing 121 stock_ids

d2_long <- reshape2::melt(select(d2, -Bbmsy_toUse), id.vars = c("stockid", "stocklong", "tsyear", "CtoUse"))
d2_long <- inner_join(d2_long, d2[,c("stockid", "tsyear", "Bbmsy_toUse")]) %>%
  rename(b2bmsy = value, method = variable)

ram_fits <- d2_long
ram_fits <- ram_fits %>% rename(b2bmsy_true = Bbmsy_toUse, catch = CtoUse)
ram_fits <- inner_join(ram_fits, ram_sp)

ram_fits <- arrange(ram_fits, method, stockid, tsyear)

# to match the simulation data frame:
ram_fits <- rename(ram_fits,
  b_bmsy_true = b2bmsy_true,
  b_bmsy_est = b2bmsy)

saveRDS(ram_fits, "generated-data/ram_fits.rds")
