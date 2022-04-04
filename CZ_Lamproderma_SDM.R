rm(list=ls())
library(sdm)
installAll() #Run this only 1 time after installing sdm package
library(sp)
library(raster)
library(dismo)
library(rgeos)
library(rgdal)
library(maptools)
library(rJava)
library(gbm)

CzechBorders <- shapefile("CZE_adm0.shp") #loading the territoty 
#loading predictors from file (10 arcmin)
bio1 <- raster("wc2.0_bio_10m_01.tif")
bio2 <- raster("wc2.0_bio_10m_02.tif")
bio3 <- raster("wc2.0_bio_10m_03.tif")
bio4 <- raster("wc2.0_bio_10m_04.tif")
bio5 <- raster("wc2.0_bio_10m_05.tif")
bio6 <- raster("wc2.0_bio_10m_06.tif")
bio7 <- raster("wc2.0_bio_10m_07.tif")
bio8 <- raster("wc2.0_bio_10m_08.tif")
bio9 <- raster("wc2.0_bio_10m_09.tif")
bio10 <- raster("wc2.0_bio_10m_10.tif")
bio11 <- raster("wc2.0_bio_10m_11.tif")
bio12 <- raster("wc2.0_bio_10m_12.tif")
bio13 <- raster("wc2.0_bio_10m_13.tif")
bio14 <- raster("wc2.0_bio_10m_14.tif")
bio15 <- raster("wc2.0_bio_10m_15.tif")
bio16 <- raster("wc2.0_bio_10m_16.tif")
bio17 <- raster("wc2.0_bio_10m_17.tif")
bio18 <- raster("wc2.0_bio_10m_18.tif")
bio19 <- raster("wc2.0_bio_10m_19.tif")
files <- list.files(pattern="tif")
predictors <- stack(files)
names (predictors)
#territory of interest
CzechBorders <- shapefile("CZE_adm0.shp")
save(CzechBorders, file = "Ukraine.Rdata")
plot(CzechBorders)
#our data (created shp file in QGIS from csv)
points<- shapefile("czhlenky.shp")
#crop predictord by our territory
predictors <- crop(predictors, CzechBorders)
plot(predictors)
plot(predictors, 1) #show plot of the first predictor
points (points, col="blue")
# Maxent model
system.file("java", package="dismo") # dismo requires maxent.java for running
#or put a file here
#(paste(system.file(package="dismo"), "/java/maxent.jar", sep=""))
m2_maxent <- maxent(predictors, points@coords)
maxentpred <- predict(predictors, m2_maxent)
plot(maxentpred)
points(points, col="blue")
backgvals <- raster::extract(predictors, backgr)
pr <- c(rep(1, nrow(presvals)), rep(0, nrow(backgvals)))
sdmdata <- data.frame(cbind(pr, rbind(presvals, backgvals)))
View(sdmdata)
sdmdata = na.omit(sdmdata)
rm(pr, files)
# Dividing into train and test subsamples for both data and background
group <- kfold(points, 5) # randomly divide the data into 5 equal parts
pres_train <- points[group != 1, ] # 80% (parts 1-5) as training data
pres_test <- points[group == 1, ] # 20% (part 1) as testing data
group1 <- kfold(backgr, 5)
backg_train <- backgr[group1 != 1, ]
backg_test <- backgr[group1 == 1, ]
plot(predictors, 1)
points(pres_train, pch= "+", cex=1.5, col="red")
points(backg_train, pch="-", cex=1.5, col="blue")
points(pres_test, pch="+", cex=1.5, col="maroon")
points(backg_test, pch="-", cex=1.5, col="black")
dev.off()
#Maxent model
m2_maxent <- maxent(predictors, pres_train)
response(m2_maxent)
par(mfrow=c(1,1))
plot(m2_maxent)
e_maxent <- evaluate(pres_test, backg_test, m2_maxent, predictors)
e_maxent
plot(e_maxent, "ROC")
maxentpred <- predict(predictors, m2_maxent)
plot (maxentpred)
tr <- threshold(e_maxent, "spec_sens")
plot(maxentpred > tr, main="presence/absence")
library(ggplot2)  
maxentpred.df =  as.data.frame(maxentpred, xy=TRUE)
head(maxentpred.df)
library(sf)
CZ = st_as_sf(CzechBorders)
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=maxentpred.df)+
  geom_sf(fill='transparent', colour = "white", data=CZ)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "Maxent",
       caption = "Data: Bioclim 10 min, Lamproderma sp.")+
  theme_minimal()
library(corrplot)
