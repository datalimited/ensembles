cdat <- read.csv("raw-data/RAM_bmsy_Ctousev4.csv", stringsAsFactors=FALSE)
names(cdat)[names(cdat)=="tsyear"]<-"year"
names(cdat)[names(cdat)=="stockid"]<-"stock"
names(cdat)[names(cdat)=="CtoUse"]<-"catch"
names(cdat)[names(cdat)=="Bbmsy_toUse"]<-"Assess.b2bmsy"
cdat<-cdat[order(cdat$stock, cdat$year),]
all.stocks<-unique(cdat$stock)
all.stocks<-all.stocks[order(all.stocks)]

## old prior
setwd("../ensembles-old/RAM_Legacy_fits/All_fits/")
load("../CMSY/cmsy_rlegacy_results_table_v0.RData")
cmsy.oldprior.df<-cmsy.rlegacy.df0
## new prior
load("../CMSY/cmsy_rlegacy_results_table_v1.RData")
cmsy.newprior.df<-cmsy.rlegacy.df1

##--------
## COMSIR
##--------
load("../COMSIR/comsir_rlegacy_results_table_v0.RData")
comsir.df<-comsir.rlegacy.df0

##--------
## SSCOM
##--------
load("../SSCOM/sscom_rlegacy_results_table_v0.RData")
sscom.df<-sscom.rlegacy.df0

##--------
## Costello
##--------
costello.df0<-read.csv("../Costello/Costello_rlegacy_results.csv")
costello.df<-merge(costello.df0, unique(cdat[,c("stock","stocklong")]), by="stocklong", all.x=TRUE)

## make a master data frame
## change some names
names(cmsy.oldprior.df)[names(cmsy.oldprior.df)=="yr"]<-"year"
names(cmsy.oldprior.df)[names(cmsy.oldprior.df)=="stock_id"]<-"stock"
names(cmsy.oldprior.df)[names(cmsy.oldprior.df)=="b_bmsy"]<-"b2bmsy"
cmsy.oldprior.df$method<-"CMSY_old_prior"
##
names(cmsy.newprior.df)[names(cmsy.newprior.df)=="yr"]<-"year"
names(cmsy.newprior.df)[names(cmsy.newprior.df)=="stock_id"]<-"stock"
names(cmsy.newprior.df)[names(cmsy.newprior.df)=="b_bmsy"]<-"b2bmsy"
cmsy.newprior.df$method<-"CMSY_new_prior"
##
names(costello.df)[names(costello.df)=="yr"]<-"year"
names(costello.df)[names(costello.df)=="BvBmsy"]<-"b2bmsy"
costello.df$method<-"Costello"
##
names(comsir.df)[names(comsir.df)=="stock_id"]<-"stock"
names(comsir.df)[names(comsir.df)=="b_bmsy"]<-"b2bmsy"
names(comsir.df)[names(comsir.df)=="yr"]<-"year"
comsir.df$method<-"COMSIR"
##
names(sscom.df)[names(sscom.df)=="yr"]<-"year"
names(sscom.df)[names(sscom.df)=="stock_id"]<-"stock"
names(sscom.df)[names(sscom.df)=="b_bmsy"]<-"b2bmsy"
sscom.df$method<-"SSCOM"
keep.names<-c("stock","year","b2bmsy","method")
all.fits.df<-rbind(cmsy.oldprior.df[,keep.names],
  cmsy.newprior.df[,keep.names],
  costello.df[,keep.names],
  comsir.df[,keep.names],
  sscom.df[,keep.names]
)
all.fits.df$stock <- as.character(all.fits.df$stock)
setwd("../../../ensembles/")
saveRDS(all.fits.df, "generated-data/ram-orig-fits.rds")
