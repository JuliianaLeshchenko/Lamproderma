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
plot(significantPred[[1]])
points (occurrence,cex=0.5,pch=16)
library(mapview)
proj4string(occurrence) <- projection(raster())
mapview(occurrence)
shape <- readOGR(dsn = ".", layer = "allpoints")
points (shape, col="blue")
sp <- data.frame(shape@coords)
occ <- rep(1,nrow(sp))
occurrence <-cbind(sp, occ)
names(occurrence)[1] <- "lon"
names(occurrence)[2] <- "lat"
rm(occ, sp)
library(sdm)
coordinates(occurrence) <- c("lon", "lat")
d <- sdmData(fromula=occ~.,train = occurrence, predictors = significantPred, bg =
               list(method="gRandom", n=3000))
m <- sdm(occ~.,d, methods="Bioclim", replication = "boot",
         test.percent=20, n=5)
m
p1 <- predict(m,significantPred,filename="pr1.img", overwrite=TRUE)

###MAXENT
rm(list = ls())
# Loadind tools and data
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
bio1 <- raster("wc2.1_10m_bio_1.tif")
bio2 <- raster("wc2.1_10m_bio_2.tif")
bio3 <- raster("wc2.1_10m_bio_3.tif")
bio4 <- raster("wc2.1_10m_bio_4.tif")
bio5 <- raster("wc2.1_10m_bio_5.tif")
bio6 <- raster("wc2.1_10m_bio_6.tif")
bio7 <- raster("wc2.1_10m_bio_7.tif")
bio8 <- raster("wc2.1_10m_bio_8.tif")
bio9 <- raster("wc2.1_10m_bio_9.tif")
bio10 <- raster("wc2.1_10m_bio_10.tif")
bio11 <- raster("wc2.1_10m_bio_11.tif")
bio12 <- raster("wc2.1_10m_bio_12.tif")
bio13 <- raster("wc2.1_10m_bio_13.tif")
bio14 <- raster("wc2.1_10m_bio_14.tif")
bio15 <- raster("wc2.1_10m_bio_15.tif")
bio16 <- raster("wc2.1_10m_bio_16.tif")
bio17 <- raster("wc2.1_10m_bio_17.tif")
bio18 <- raster("wc2.1_10m_bio_18.tif")
bio19 <- raster("wc2.1_10m_bio_19.tif")
elev <-raster('wc2.1_10m_elev.tif')
predictors <- stack(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,bio11,bio12,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19,elev)
names (predictors)
oc <- readOGR(dsn = ".", layer = "greatdataset")
EU <- readOGR(dsn = ".", layer = "BiogeoRegions2016")
summary(EU)
EUt <- spTransform(EU, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# Cleaning data
str(oc)
#View(oc@data)
oc <- remove.duplicates(oc)
plot(predictors)
plot(predictors, 1)
points (oc, col="blue")
#dev.off()
m2_maxent <- maxent(predictors, oc@coords)
maxentpred <- predict(predictors, m2_maxent)
plot(maxentpred)
points(oc, col="blue")
# Simulating background data
presvals <- raster::extract(predictors, oc)
backgr <- spsample (oc, 1000, type = "random")
plot(predictors, 1)
points (oc, col="blue")
points (backgr, col="red")
backgvals <- raster::extract(predictors, backgr)
pr <- c(rep(1, nrow(presvals)), rep(0, nrow(backgvals)))
sdmdata <- data.frame(cbind(pr, rbind(presvals, backgvals)))
View(sdmdata)
sdmdata = na.omit(sdmdata)
rm(pr, files)
group <- kfold(oc, 5) # randomly divide the data into 5 equal parts
pres_train <- oc[group != 1, ] # 80% (parts 1-5) as training data
pres_test <- oc[group == 1, ] # 20% (part 1) as testing data
group1 <- kfold(backgr, 5)
backg_train <- backgr[group1 != 1, ]
backg_test <- backgr[group1 == 1, ]
plot(predictors, 1)
points(pres_train, pch= "+", cex=1.5, col="red")
points(backg_train, pch="-", cex=1.5, col="blue")
points(pres_test, pch="+", cex=1.5, col="maroon")
points(backg_test, pch="-", cex=1.5, col="black")
m2_maxent <- maxent(predictors, pres_train)
#visualy evaluate the significance of the predictors
plot(m2_maxent)
predictors_new = predictors[[c()]] # specify layers' numbers
e_maxent <- evaluate(pres_test, backg_test, m2_maxent, predictors)
e_maxent
plot(e_maxent, "ROC")
maxentpred <- predict(predictors, m2_maxent)
plot (maxentpred)
tr <- threshold(e_maxent, "spec_sens")
plot(maxentpred > tr, main="presence/absence")
# Visualisation
maxentpred.df =  as.data.frame(maxentpred, xy=TRUE) %>% drop_na()
head

# Plotting
EU_plot = st_as_sf(EUt) # convert to simple features object
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=maxentpred.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "Maxent",
       caption = "Data: Bioclim and elevation 10 min, GBIF, literature extracted points, L. ovoideum, L. aeneum, L. preudomaculatum, L. splendens, L. arcyrioides, L. sauteri")+
  theme_minimal()

