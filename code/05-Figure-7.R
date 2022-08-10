## Resource Stability - Figure 7

library(raster)
library(sp) 
library(rgdal) 
library(rgeos)
library(maptools) 

library(sdStaf)
library(rangeBuilder)
library(scales)

library(wesanderson)
library(colorRamps)

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

# colour palette
pal <- wes_palette(25, name = "Zissou1", type = "continuous")

# load Harpy Eagle SDMs
a <- raster("./results/SDM_A.tif")
b <- raster("./results/SDM_B.tif")
ab <- raster("./results/SDM_AB.tif")
spr <- raster("./results/SDM_ASR.tif")

# stack
st <- stack(b, ab, spr)
st
plot(st, col=pal)

# weighted average based on CBI scores from paper
w <- c(0.709, # B
       0.766, # AB
       0.788) # ASR

# CBI-based ensemble
wm <- weighted.mean(x = st, 
                    w = w,
                    na.rm = TRUE)
wm
plot(wm, col=pal)

# load prey distribution models
c <- raster("./results/capuchin.tif")
h <- raster("./results/howler_monkey.tif")
p <- raster("./results/tree_porcupine.tif")
s2 <- raster("./results/two-toed_sloth.tif")
s3 <- raster("./results/three-toed_sloth.tif")

# stack
prey <- stack(c, h, p, s2, s3)
prey
plot(prey, col=pal)

prey <- resample(prey, wm, method = "bilinear")

# resource stability
rs <- stability(current = wm,
                project = prey,
                thr.value = 0.5, 
                continue = TRUE)
print(rs)
str(rs)

plot(shp, 
     axes=T,
     xlab="Longitude", 
     ylab="Latitude", #las=1,
     col="lightgrey", 
     border = "grey70",
     main = "")

plot(rs@map, 
     add=T,
     col=alpha(matlab.like(25), 0.9),
     zlim=c(-2.0,2.0), 
     legend=F)

addRasterLegend(r=rs@map, 
                direction="vertical",
                side=4,
                location=c(-95,-93, -30,0),
                nTicks=3,
                shortFrac = 0.02, 
                longFrac = 0.3, 
                axisOffset = 0, 
                border = TRUE,
                ramp=matlab.like(25),
                isInteger = 'auto', 
                #ncolors = 100, 
                breaks = NULL,
                minmax = NULL, 
                locs = NULL, 
                cex.axis = 0.75, 
                labelDist = 0.7, 
                digits=0)

# Niche Overlap
library(ENMeval)

s1 <- stack(a, b, ab, spr)

no1 <- calc.niche.overlap(predictive.maps = s1, stat = "D")
no1

s2 <- stack(a, prey)

no2 <- calc.niche.overlap(predictive.maps = s2, stat = "D")
no2
