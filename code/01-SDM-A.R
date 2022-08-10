## ABIOTIC MODEL

library(raster)
library(sp) 
library(rgdal) 
library(rgeos) 
library(gtools)

library(glmnet) 
library(maxnet)
library(ENMeval)
library(enmSdm) 
library(dismo) 
library(ENMGadgets) 

library(maptools) 
library(wesanderson)

library(omnibus)
library(statisfactory)
library(legendary)

setwd("C:/HAEA-BioticSDM")
rasterOptions(progress="text")
memory.limit()

# load filtered occurrence points
locs <- read.csv("./data/HAEA_locs_5km.csv")
head(locs)
nrow(locs)

# load countries shapefile
shp <- readOGR(dsn="C:/HAEA-BioticSDM/data", layer="countries")
summary(shp)
#plot(shp, col="lightgrey", border = "grey70")

# load range extent shapefile
ext <- readOGR(dsn="C:/HAEA-BioticSDM/data", layer="Biomes_extent")
summary(ext)
#plot(ext, add=T, col=alpha("darkkhaki",0.5), border="darkkhaki")
#points(locs$lon, locs$lat, col="steelblue", pch=20, cex=0.25)

# load environmental raster layers
list.env <- mixedsort(list.files(paste(getwd(), "/data/", sep = ""), 
                                 full.names = TRUE, 
                                 pattern = ".tif"))
list.env 

# stack raster files
env <- stack(list.env) 
env
names(env)
plot(env)

# assign CRS
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs ")
proj4string(env) <- crs.geo 
env

# Test multi-collinearity
library(fuzzySim)

locs <- locs[, -1]
dim(locs)
head(locs)

# extract env data at species locations
bio <- extract(x=env[[1:3]], 
               y=locs[,1:2],
               method="simple")
head(bio)

# combine into dataframe
df <- cbind(locs, bio)
head(df)
nrow(df)
ncol(df)

# remove NAs
df <- na.omit(df)

# check VIF scores
multicol(vars=df[, 3:5], reorder=TRUE)

# test using VIF and Spearman correlation
VIF <- corSelect(data=df,
                 sp.cols=NULL,
                 var.cols=3:5,
                 cor.thresh=0.7,
                 select="VIF", 
                 method="spearman", 
                 use="everything")
VIF

# Mantel test
library(ecospat)

ecospat.mantel.correlogram(dfvar = df,
                           colxy = 1:2,
                           n = 100,
                           colvar = 3:5,
                           max = 1000,
                           nclass = 50,
                           nperm = 999)

# Evaluate best supported model using AICc

# data frame of lon/lat
occ <- cbind.data.frame(locs$lon, locs$lat) 
head(occ)

# run models
evalc2 <- ENMevaluate(occ=occ, 
                      env=env, 
                      RMvalues=seq(1, 5, 0.5), 
                      fc=c("L", "LQ"),
                      algorithm="maxnet",
                      method='checkerboard2',
                      aggregation.factor=c(5,5),
                      n.bg=10000, 
                      rasterPreds=TRUE,
                      clamp=TRUE, 
                      parallel = TRUE, 
                      numCores = 2)

evalc2
aicmods <- which(evalc2@results$AICc == min(na.omit(evalc2@results$AICc)))
evalc2@results[aicmods,]
evalc2@results

# Fit penalized logistic regression model

# withold a 20% sample for testing
fold <- kfold(x=occ, k=5)
test <- occ[fold == 1, ]
train <- occ[fold != 1, ]

# extract environment values for occurrences
TrainEnv <- extract(x=env, y=train)
head(TrainEnv)

# generate random background/absence points 
bg <- randomPoints(mask=env, 
                   n=10000,
                   p=occ,
                   ext=ext,
                   #extf=1.1,
                   #excludep=TRUE,
                   #prob=TRUE,
                   lonlatCorrection=TRUE)
head(bg)
nrow(bg)

# extract values from model environment for bg points
absvals <- extract(x=env, y=bg)
head(absvals)

# presence/absence column for background points
presabs <- c(rep(1, nrow(TrainEnv)), rep(0, nrow(absvals)))
head(presabs)

# combine presence/absence data with model environment values
sdmdata <- data.frame(cbind(presabs, rbind(TrainEnv, absvals)))
head(sdmdata)

# check for NA values
row.has.na <- apply(sdmdata, 1, function(x){any(is.na(x))})
sum(row.has.na)

# remove NA values
sdmdata <- sdmdata[!row.has.na,]
nrow(sdmdata)
head(sdmdata)

# subset of that dataset without the presence and absence values
env.data <- sdmdata[ ,-1]
head(env.data)

# run maxnet function to fit the SDM using glmnet
mn <- maxnet(p = sdmdata$presabs, 
             data = env.data,
             f = maxnet.formula(p = sdmdata$presabs, 
                                data = env.data, 
                                classes = "lq"),
             regmult = 1)

summary(mn)
round(mn$betas,3)

# save model as RData
#saveRDS(mn, file = './results/SDM_A.RData')

# Plot response curves
plot(x=mn,
     common.scale=TRUE,
     type="cloglog",
     ylab="Predicted value",
     mod=mn)

par(mfrow=c(1,1))

# predict
pred.mn <- maxnet.predictRaster(mod=mn,
                                env=env, 
                                type="cloglog", 
                                clamp=TRUE)
pred.mn

writeRaster(pred.mn, 
            "./results/SDM_A.tif", 
            overwrite=TRUE)

# colour palette
pal <- wes_palette(100, name = "Zissou1", type = "continuous")

plot(shp, col="lightgrey", border = "grey70")

plot(x=pred.mn, 
     add=T,
     col=pal,
     xlab="Longitude", 
     ylab="Latitude", las=1,
     zlim=c(0.0,1.0), 
     legend.only=F) 

#plot(ext, add=TRUE, border = "grey50", lwd=0.5)

points(locs$lon, locs$lat, pch=21, col="black",cex=0.5)


# Model evaluation

# pROC
library(ENMGadgets)

# data frame of lat/lon
#occ <- cbind.data.frame(locs$lon, locs$lat) 

# add an index that makes k (ten) random groups of observations
#fold <- kfold(occ, k=5) 

# hold out one fifth as test data 
#test <- occ[fold == 1, ] 
nrow(test)

write.csv(test, "./pROC/testdata_A.csv")
writeRaster(pred.mn, "./pROC/SDM_A.asc", overwrite=TRUE)

proc <- PartialROC(PresenceFile="./pROC/testdata_A.csv", 
                   PredictionFile="./pROC/SDM_A.asc", 
                   OmissionVal=0.1, 
                   RandomPercent=50, 
                   NoOfIteration=1000,
                   OutputFile="./pROC/TestRoc_SDM_A.csv") 

head(proc)
mean(proc$AUC_ratio)  
sd(proc$AUC_ratio)    
range(proc$AUC_ratio) 

# Continuous Boyce Index (CBI)
predPres <- extract(x = pred.mn, 
                    y = cbind(locs$lon, locs$lat),
                    na.rm=TRUE)

predPres <- na.omit(predPres)
head(predPres)
min(predPres)

bg <- as.data.frame(bg)

predBg <- extract(x = pred.mn, y = cbind(bg$x, bg$y))
predBg <- na.omit(predBg)
head(predBg)

# get proportion of presences or background sites in each bin
presDistrib <- hist(predPres, 
                    plot=FALSE, 
                    breaks=seq(0, 1, length.out=11))$counts

bgDistrib <- hist(predBg, 
                  plot=FALSE, 
                  breaks=seq(0, 1, length.out=11))$counts

# convert to proportion of sites
presDistrib <- presDistrib / sum(presDistrib)
bgDistrib <- bgDistrib / sum(bgDistrib)

# P/E plot
pe <- presDistrib / bgDistrib

plot(pe,
     xlab='Prediction class (1 to 10)',
     ylab='Predicted-to-Expected Ratio', 
     main='P/E Plot',
     col='red', 
     pch=16)

cbi <- contBoyce(pres=predPres, 
                 bg=predBg,
                 numBins = 101, 
                 binWidth = 0.1,
                 autoWindow = TRUE, 
                 method = "spearman", 
                 dropZeros = TRUE,
                 graph = FALSE,
                 na.rm=TRUE)
cbi


# CBI Random k-fold cross-validation

# calculate k-folds for presences and background sites
head(df)
kPres <- kfold(x = df, k=5)

bg2 <- cbind(bg, absvals)
head(bg2)
nrow(bg2)

kBg <- kfold(x = bg2, k=5)

head(kPres)
head(kBg)

# map
plot(shp)
points(df$lon, df$lat)

points(df$lon[kPres==1],
       df$lat[kPres==1],
       bg='red',
       pch=21)

legend('topright',
       legend=c('Training presence', 'Test presence'),
       pch=c(1, 16),
       col=c('black', 'red'),
       bg='white',
       cex=0.8)


# for storing CBI
cbiRandom <- rep(NA, 5)
names(env)

predictors <- c("ClimaticMoistureIndex", 
                "MinTempWarmestMonth",
                "TRI")

# cycle through each k-fold
for (i in 1:5) {
  
  omnibus::say('K-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0 for
  # presence/background... using only points not in this k-fold
  envData <- rbind(df[kPres!=i, predictors],bg2[kBg!=i, predictors])
  
  head(envData)
  
  presBg <- c(rep(1, sum(kPres!=i)), rep(0, sum(kBg!=i)))
  head(presBg)
  
  trainData <- cbind(presBg, envData)
  head(trainData)
  
  # check for NA values
  row.has.na <- apply(trainData, 1, function(x){any(is.na(x))})
  sum(row.has.na)
  
  # remove NA values
  trainData.filtered <- trainData[!row.has.na,]
  trainData <- trainData.filtered
  
  # tuned model
  model <- trainMaxNet(data=trainData,
                       regMult=1,
                       classes='lq',
                       verbose=FALSE
  )
  
  # predict to presences and background sites
  predPres <- raster::predict(model, 
                              newdata=df[kPres==i, ], 
                              type='cloglog')
  
  predBg <- raster::predict(model, 
                            newdata=bg2[kBg==i, ], 
                            type='cloglog')
  
  # evaluate and remember result
  thisEval <- evaluate(p=as.vector(predPres), a=as.vector(predBg))
  
  thisCbi <- contBoyce(pres=predPres, 
                       bg=predBg,
                       numBins = 101, 
                       binWidth = 0.1,
                       autoWindow = TRUE, 
                       method = "spearman", 
                       dropZeros = TRUE,
                       graph = FALSE,
                       na.rm = TRUE)
  
  omnibus::say(': CBI = ', round(thisCbi, 3))
  
  cbiRandom[i] <- thisCbi
  
}

omnibus::say('Mean CBI:', round(mean(cbiRandom), 3))
