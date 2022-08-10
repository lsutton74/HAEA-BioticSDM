## BIOTIC MODEL

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

# load countries shapefile
shp <- readOGR(dsn="C:/HAEA-BioticSDM/data", layer="countries")
summary(shp)
#plot(shp, col="lightgrey", border = "grey70")

# load range extent shapefile
ext <- readOGR(dsn="C:/HAEA-BioticSDM/data", layer="Biomes_extent")
summary(ext)
#plot(ext, add=T, col=alpha("darkkhaki",0.5), border="darkkhaki")

# load prey distribution models
capuchin <- raster("./results/capuchin.tif")
howler <- raster("./results/howler_monkey.tif")
slo3 <- raster("./results/three-toed_sloth.tif")
porcupine <- raster("./results/tree_porcupine.tif")
slo2 <- raster("./results/two-toed_sloth.tif")

# Species Richness raster

# Standardize raster values
c <- scale(capuchin)
h <- scale(howler)
s3 <- scale(slo3)
p <- scale(porcupine)
s2 <- scale(slo2)

# stack
prey <- stack(h, c, p, s2, s3)
prey
plot(prey)

# check collinearity
usdm::vif(x=prey)
usdm::vifstep(x=prey, th=10)

# sum of prey models 
spr <- sum(prey)
spr

writeRaster(spr, "./results/SpeciesRichness.tif")

# colour palette
pal <- wes_palette(100, name = "Zissou1", type = "continuous")

plot(spr, 
     xlab="Longitude", 
     ylab="Latitude", las=1,
     zlim=c(-8, 14), 
     legend.only=FALSE,
     col = pal)

# load scaled species richness raster
spr <- raster("./results/SpeciesRichness.tif")

# stack all biotic rasters
bio <- stack(capuchin, howler, slo3, porcupine, slo2, spr)
bio
plot(bio)

# assign CRS
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs ")
proj4string(bio) <- crs.geo 
bio

# Evaluate best supported model using AICc

# load filtered Harpy Eagle occurrence points
locs <- read.csv("./data/HAEA_locs_5km.csv")
head(locs)
nrow(locs)

# data frame of lat/lon 
occ <- cbind.data.frame(locs$lon, locs$lat) 
head(occ)

evalc2 <- ENMevaluate(occ=occ, 
                      env=bio, 
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

# create training environment
TrainEnv <- extract(x=bio, y=train)
head(TrainEnv)

# background data
bg <- randomPoints(mask=bio, 
                   n=10000,
                   p=occ,
                   ext=ext,
                   #extf=1.1,
                   #excludep=TRUE,
                   #prob=TRUE,
                   lonlatCorrection=TRUE)
head(bg)
nrow(bg)

# extract values for model environment for bg points
absvals <- extract(x=bio, y=bg)
head(absvals)

# presence/absence column 
presabs <- c(rep(1, nrow(TrainEnv)), rep(0, nrow(absvals)))
head(presabs)

# dataset that combines the presence and absence data
sdmdata <- data.frame(cbind(presabs, rbind(TrainEnv, absvals)))
head(sdmdata)

# check for NA values
row.has.na <- apply(sdmdata, 1, function(x){any(is.na(x))})
sum(row.has.na)

# remove NA values
sdmdata <- sdmdata[!row.has.na,]
head(sdmdata)

# subset of that dataset without the presence and absence values
env.data <- sdmdata[,-1]
head(env.data)

# run maxnet function to fit the SDM using glmnet
mn <- maxnet(p=sdmdata$presabs, 
             data=env.data,
             f=maxnet.formula(p=sdmdata$presabs, 
                              data=env.data, 
                              classes="lq"),
             regmult=1)

summary(mn)
round(mn$betas,3)

# save model as RData
#saveRDS(mn, file = './results/SDM_B.RData')

# Plot response curves
plot(x=mn,
     common.scale=TRUE,
     type="cloglog",
     ylab="Predicted value",
     mod=mn)

par(mfrow=c(1,1))

# predict
pred.mn <- maxnet.predictRaster(mod=mn,
                                env=bio, 
                                type="cloglog", 
                                clamp=TRUE)
pred.mn

writeRaster(pred.mn, 
            "./results/SDM_B.tif", 
            overwrite=TRUE)

plot(shp, col="lightgrey", border = "grey70")

plot(x=pred.mn, 
     add=T,
     col=pal,
     xlab="Longitude", 
     ylab="Latitude", las=1,
     zlim=c(0.0,1.0), 
     legend.only=F) 

#plot(ext, add=TRUE, border = "grey30", lwd=0.5)

points(locs$lon, 
       locs$lat, 
       pch=21, 
       col="black",
       cex=0.5) 


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

write.csv(test, "./pROC/testdata_B.csv")
writeRaster(pred.mn, "./pROC/SDM_B.asc", overwrite=TRUE)

proc <- PartialROC(PresenceFile="./pROC/testdata_B.csv", 
                   PredictionFile="./pROC/SDM_B.asc", 
                   OmissionVal=0.1, 
                   RandomPercent=50, 
                   NoOfIteration=1000,
                   OutputFile="./pROC/TestRoc_SDM_B.csv") 

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

predBg <- extract(x = pred.mn, 
                  y = cbind(bg$x, bg$y))

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

# CBI Random cross-validation

locs <- locs[, -1]
dim(locs)
head(locs)

df <- extract(x=bio, y=occ)
head(df)

df <- cbind(locs, df)
head(df)

# calculate k-folds for presences and background sites
kPres <- kfold(x = df, k=5)

bg1 <- cbind(bg, absvals)
head(bg1)
nrow(bg1)

kBg <- kfold(x = bg1, k=5)

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
names(bio)

predictors <- c("SpeciesRichness",
                "two.toed_sloth",        
                "three.toed_sloth",
                "howler_monkey",
                "capuchin",
                "tree_porcupine")


# cycle through each k-fold
for (i in 1:5) {
  
  omnibus::say('K-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0 for
  # presence/background... using only points not in this k-fold
  envData <- rbind(
    df[kPres!=i, predictors],
    bg1[kBg!=i, predictors]
  )
  
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
  model <- trainMaxNet(
    data=trainData,
    regMult=1,
    classes='lq',
    verbose=FALSE
  )
  
  # predict to presences and background sites
  predPres <- raster::predict(model, 
                              newdata=df[kPres==i, ], 
                              type='cloglog')
  
  predBg <- raster::predict(model, 
                            newdata=bg1[kBg==i, ], 
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

omnibus::say('Mean CBI:', round(mean(cbiRandom), 3)) # Mean CBI:0.709
