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
# SDM preparation
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
# or bio<-raster::getData('worldclim', var='bio',res=2.5) and crop in R
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

#code run 2022-12-16
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
plot(bio1)
preds <- stack(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19) # making a raster object
preds # see the specification of the raster layers (e.g., cell size, extent, etc.)
class(preds)
rm(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19)
memory.limit()
memory.limit(size=20000)
library(usdm)
v <- vifcor (preds, th=0.9)
v
significantPred <- exclude(preds, v)
significantPred
plot(significantPred[[1]])
points (occurrence,cex=0.5,pch=16)
library(mapview)
proj4string(occurrence) <- projection(raster())
mapview(occurrence)
shape <- readOGR(dsn = ".", layer = "greatdataset")
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
getmethodNames()
m <- sdm(occ~.,d, methods=c("Bioclim", 'rf', 'maxent'), replication = "boot",
         test.percent=20, n=100)
m
gui(m)
p1 <- predict(m,significantPred)
plot(p1)
#bioclim
en_2 <- ensemble(m, significantPred, filename = "en_2.img",
               setting=list(id=c(1:100),method="weighted", stat="tss", opt=2))

plot(en_2)
points(occurrence,col="red")
enseb.df = as.data.frame(en_2, xy=TRUE)
head(enseb.df)
EU <- readOGR(dsn = ".", layer = "BiogeoRegions2016")
EUt <- spTransform(EU, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
EU_plot = st_as_sf(EUt) # convert to simple features object
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=enseb.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "Bioclim",
       caption = "Data: Bioclim and elevation 10 min, GBIF, literature extracted points, L. ovoideum, L. aeneum, L. preudomaculatum, L. splendens, L. arcyrioides, L. sauteri")+
  theme_minimal()
ggsave("Bioclim_10_min_12_preds.png", width=14, height=10, unit="cm", dpi=300)
#rf
rf <- ensemble(m, significantPred, filename = "rf1.img", overwrite = TRUE,
                 setting=list(id=c(101:200),method="weighted", stat="tss", opt=2))

plot(rf)
points(occurrence,col="red")
ensebRF.df = as.data.frame(rf, xy=TRUE)
head(enseb.df)
EU <- readOGR(dsn = ".", layer = "BiogeoRegions2016")
EUt <- spTransform(EU, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
EU_plot = st_as_sf(EUt) # convert to simple features object
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=ensebRF.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "Random Forest",
       caption = "Data: Bioclim and elevation 10 min, GBIF, literature extracted points, L. ovoideum, L. aeneum, L. preudomaculatum, L. splendens, L. arcyrioides, L. sauteri")+
  theme_minimal()
ggsave("RF_10_min_12_preds.png", width=14, height=10, unit="cm", dpi=300)
#maxent
maxent <- ensemble(m, significantPred, filename = "maxent.img", overwrite = TRUE,
               setting=list(id=c(201:300),method="weighted", stat="tss", opt=2))
plot(maxent)
ensebmax.df = as.data.frame(maxent, xy=TRUE)
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=ensebmax.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "MaxEnt",
       caption = "Data: Bioclim and elevation 10 min, GBIF, literature extracted points, L. ovoideum, L. aeneum, L. preudomaculatum, L. splendens, L. arcyrioides, L. sauteri")+
  theme_minimal()
ggsave("MaxEnt_10_min_12_preds.png", width=14, height=10, unit="cm", dpi=300)
###future prediction
biof <- raster::getData('CMIP5', var='bio', res=10, rcp=85,year=70,model='AC')
names(biof)<-names(preds)
pf<-predict(m,biof)
enf<-calc(pf,mean)

#separate datasets
shape_cryo <- readOGR(dsn = ".", layer = "cryophilous")
sp1 <- data.frame(shape_cryo@coords)
occ1 <- rep(1,nrow(sp1))
occurrence1 <-cbind(sp1, occ1)
names(occurrence1)[1] <- "lon"
names(occurrence1)[2] <- "lat"
coordinates(occurrence1) <- c("lon", "lat")
library(sdm)
d_cryo <- sdmData(fromula=occ1~.,train = occurrence1, predictors = significantPred, bg =
               list(method="gRandom", n=3000))
d_cryo
m_cryo <- sdm(occ1~.,d_cryo, methods=c("Bioclim", 'rf', 'maxent'), replication = "boot",
         test.percent=20, n=100)
gui(m_cryo)
p_cryo <- predict(m_cryo,significantPred)
plot(p_cryo)
#bioclim_cryo
en_cryo <- ensemble(m_cryo, significantPred, filename = "en_cryo.img",
                 setting=list(id=c(1:100),method="weighted", stat="tss", opt=2))
plot(en_cryo)
ensebmax_cryo.df = as.data.frame(en_cryo, xy=TRUE)
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=ensebmax_cryo.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "Bioclim",
       caption = "Data: Bioclim and elevation 10 min, GBIF, literature extracted points, L. arcyrioides, L. sauteri")+
  theme_minimal()
ggsave("RF_cryo_10_min_12_preds.png", width=14, height=10, unit="cm", dpi=300)
#rf_cryo
rf_cryo <- ensemble(m_cryo, significantPred, filename = "rf_cryo.img", overwrite = TRUE,
               setting=list(id=c(101:200),method="weighted", stat="tss", opt=2))
ggsave("RF_cryo_10_min_12_preds.png", width=14, height=10, unit="cm", dpi=300)
plot(rf_cryo)
ensebRF_cryo.df = as.data.frame(rf_cryo, xy=TRUE)
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=ensebRF_cryo.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "RF",
       caption = "Data: Bioclim and elevation 10 min, GBIF, literature extracted points, L. arcyrioides, L. sauteri")+
  theme_minimal()
ggsave("RF_cryo_10_min_12_preds.png", width=14, height=10, unit="cm", dpi=300)
#maxent_cryo
maxent_cryo <- ensemble(m_cryo, significantPred, filename = "rf_cryo.img", overwrite = TRUE,
                        setting=list(id=c(101:300),method="weighted", stat="tss", opt=2))
plot(maxent_cryo)
ensebMAxent_cryo.df = as.data.frame(maxent_cryo, xy=TRUE)
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=ensebMAxent_cryo.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "MaxEnt",
       caption = "Data: Bioclim and elevation 10 min, GBIF, literature extracted points, L. arcyrioides, L. sauteri")+
  theme_minimal()


shape_facult <- readOGR(dsn = ".", layer = "facult")
sp2 <- data.frame(shape_facult@coords)
occ2 <- rep(1,nrow(sp2))
