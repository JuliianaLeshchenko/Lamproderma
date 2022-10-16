#GBIF phenology visualisation
data <- read.csv("yourdataset.csv")
library(ggplot2)
fruiting_plot <- 
  ggplot(data, aes(month, species, fill = species)) +
  geom_tile() +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme(
    axis.title.x = element_blank(),  # no need for x-axis title
    axis.text.y  = element_text(hjust = 0)) # left-justify species names
fruiting_plot
#legend.position = "none" # Remove all legends from plot
fruiting_plot +
  # Make the plot circular
  coord_polar() +
  # Remove non-sensical y-axis details now that plot is circular
  theme(axis.text.y  = element_blank(), axis.title.y = element_blank()) #This is much more realistic.
library(dismo)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(rJava)
library(gbm)
library(raster)
library(sf)
library(ggplot2)
library(tidyr)
library(usdm)
#import data from worldclim (previsly croped in QGIS):
bio1 <- raster("bioclim1.tif")
bio2 <- raster("bioclim2.tif")
bio3 <- raster("bioclim3.tif")
bio4 <- raster("bioclim4.tif")
bio5 <- raster("bioclim5.tif")
bio6 <- raster("bioclim6.tif")
bio7 <- raster("bioclim7.tif")
bio8 <- raster("bioclim8.tif")
bio9 <- raster("bioclim9.tif")
bio10 <- raster("bioclim10.tif")
bio11 <- raster("bioclim11.tif")
bio12 <- raster("bioclim12.tif")
bio13 <- raster("bioclim13.tif")
bio14 <- raster("bioclim14.tif")
bio15 <- raster("bioclim15.tif")
bio16 <- raster("bioclim16.tif")
bio17 <- raster("bioclim17.tif")
bio18 <- raster("bioclim18.tif")
bio19 <- raster("bioclim19.tif")
plot(bio1)
preds <- stack(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19)
names (preds)
class(preds)
rm(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19)
memory.limit()
memory.limit(size=20000)
library(usdm)
v <- vifstep (preds)
v
significantPred <- exclude(preds, v)
significantPred
shape <- read_sf(dsn = ".", layer = "allpoints")
shape.df <- as(shape, "data.frame")
class(shape)
data <- read.csv(file = 'GreatDataset.csv')
coordinates(data) <- c("Lon", "Lat")
library(sdm)
d <- sdmData(fromula=occurrence~.,train = data, predictors = significantPred)
d
m <- sdm(occurrence~.,d, methods=c("Bioclim","Maxent","RF"), replication = "boot",
         test.percent=20, n=5)
m
p1 <- predict(m,significantPred,filename="pr1.img", overwrite=TRUE)
