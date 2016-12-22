library("dplyr")

# Bring in all results exept CMSY:
d1 <- readRDS("raw-data/batch1result.rds")

# bringing in corrected fits below
d1 <- dplyr::filter(d1, method_id != "CMSY")

d1$stock_id <- as.character(d1$stock_id)
d1$ID <- as.numeric(as.character(d1$ID))
d1$b_bmsyUpper <- as.numeric(as.character(d1$b_bmsyUpper))
d1$b_bmsyLower <- as.numeric(as.character(d1$b_bmsyLower))
d1$b_bmsy_iq25 <- as.numeric(as.character(d1$b_bmsy_iq25))
d1$b_bmsy_iq75 <- as.numeric(as.character(d1$b_bmsy_iq75))
d1$seed <- as.numeric(as.character(d1$seed))
d1$n_iterations <- as.numeric(as.character(d1$n_iterations))
d1$effective_sample_size <- as.numeric(as.character(d1$effective_sample_size))

load("raw-data/cmsy_fits_22052015.RData")
d2 <- cmsynewfits
rm(cmsynewfits)

# match stock_id values in the new dataset to match the old one:
d2$stock_id <- as.character(d2$stock_id)
d2$stock_id <- gsub(
  "([A-Z]+_[ID0-9]+_[A-Z]+_[A-Z0-9.]+_[A-Z0-9.]+_[UR0-9]+_[TS0-9]+)_[sR0-9]+_[sC0-9]+_[0-9]+",
  "\\1", d2$stock_id)
d2 <- dplyr::rename(d2, b_bmsy_est = b_bmsy)

d <- dplyr::bind_rows(d1, d2) %>% as.data.frame()

saveRDS(d, "generated-data/dsim.rds")
