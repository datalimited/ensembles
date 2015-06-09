library("dplyr")
library("ggplot2")

dsim <- readRDS("generated-data/dsim.rds")
dsim$method_id <- sub("COM.SIR", "COMSIR", dsim$method_id) # to match RAM fits
dsim$method_id <- sub("Costello", "mPRM", dsim$method_id)
dsim <- rename(dsim, stockid = stock_id, method = method_id) # to match RAM fits
dsim <- dsim %>%
  mutate(cv_id = paste0(stockid, "_lh_", LH)) %>%
  mutate(stockid =
      paste0(stockid, "_sigmaR_", sigmaR, "_lh_", LH)) %>%
  arrange(stockid, iter, year) # critical since not all in order

simple <- theme_bw() +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank()) +
    theme(strip.text.x = element_text(size = 5))

p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY") %>%
  ggplot(aes(year, catch, colour = as.factor(sigmaC), group = sigmaC)) +
    geom_line() +
    facet_wrap(~stockid, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-catch-by-sigmaC-iter1.pdf", width = 15, height = 13)

p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY") %>%
  ggplot(aes(year, b_bmsy_est, colour = as.factor(sigmaC), group = sigmaC)) +
    geom_line() +
    facet_wrap(~stockid, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-CMSY-by-sigmaC-iter1.pdf", width = 15, height = 13)

p <- dsim %>% filter(iter == 1, LH == "DE", method == "SSCOM") %>%
  ggplot(aes(year, b_bmsy_est, colour = as.factor(sigmaC), group = sigmaC)) +
    geom_line() +
    facet_wrap(~stockid, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-SSCOM-by-sigmaC-iter1.pdf", width = 15, height = 13)

p <- dsim %>% filter(iter == 1, LH == "DE", method == "COMSIR") %>%
  ggplot(aes(year, b_bmsy_est, colour = as.factor(sigmaC), group = sigmaC)) +
    geom_line() +
    facet_wrap(~stockid, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-COMSIR-by-sigmaC-iter1.pdf", width = 15, height = 13)

dsim <- dsim %>%
  mutate(stockid2 = paste(LH, ID, ED, SEL, UR, TS, sigmaR, sigmaC, sep = "_"))

p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY") %>%
  ggplot(aes(year, catch, colour = as.factor(AR))) + geom_line() +
    facet_wrap(~stockid2, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-catch-by-AR-iter1.pdf", width = 15, height = 13)

dsim <- dsim %>%
  mutate(stockid3 = paste(LH, ID, ED, SEL, UR, sigmaR, sigmaC, AR, sep = "_"))

p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY") %>%
  ggplot(aes(year, catch, colour = as.factor(TS))) + geom_line() +
    facet_wrap(~stockid3, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-catch-by-TS-iter1.pdf", width = 15, height = 13)

p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY") %>%
  ggplot(aes(year, b_bmsy_est, colour = as.factor(TS))) + geom_line() +
    facet_wrap(~stockid3, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-CMSY-by-TS-iter1.pdf", width = 15, height = 13)

dsim <- dsim %>%
  mutate(stockid4 = paste(LH, ID, ED, SEL, UR, TS, sigmaC, AR, sep = "_"))

p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY") %>%
  ggplot(aes(year, catch, colour = as.factor(sigmaR))) + geom_line() +
    facet_wrap(~stockid4, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-catch-by-sigmaR-iter1.pdf", width = 15, height = 13)

p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY") %>%
  ggplot(aes(year, catch, colour = as.factor(sigmaR))) + geom_line() +
    facet_wrap(~stockid4, scales = "free_y", ncol = 9) + simple
ggsave("../figs/sim-catch-by-sigmaR-iter1.pdf", width = 15, height = 13)

# p <- dsim %>% filter(iter == 1, LH == "DE", method == "CMSY", TS == 60) %>%
#   ggplot(aes(year, catch)) + geom_line() +
#     facet_wrap(~stockid, scales = "free_y", ncol = 9) + simple
#ggsave("../figs/sim-catch-by-sigmaR-iter1.pdf", width = 15, height = 13)
