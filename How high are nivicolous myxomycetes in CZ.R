#import preloaded data from worldclim:
bio1 <- raster("wc2.1_30s_bio_1.tif")
bio2 <- raster("wc2.1_30s_bio_2.tif")
bio3 <- raster("wc2.1_30s_bio_3.tif")
bio4 <- raster("wc2.1_30s_bio_4.tif")
bio5 <- raster("wc2.1_30s_bio_5.tif")
bio6 <- raster("wc2.1_30s_bio_6.tif")
bio7 <- raster("wc2.1_30s_bio_7.tif")
bio8 <- raster("wc2.1_30s_bio_8.tif")
bio9 <- raster("wc2.1_30s_bio_9.tif")
bio10 <- raster("wc2.1_30s_bio_10.tif")
bio11 <- raster("wc2.1_30s_bio_11.tif")
bio12 <- raster("wc2.1_30s_bio_12.tif")
bio13 <- raster("wc2.1_30s_bio_13.tif")
bio14 <- raster("wc2.1_30s_bio_14.tif")
bio15 <- raster("wc2.1_30s_bio_15.tif")
bio16 <- raster("wc2.1_30s_bio_16.tif")
bio17 <- raster("wc2.1_30s_bio_17.tif")
bio18 <- raster("wc2.1_30s_bio_18.tif")
bio19 <- raster("wc2.1_30s_bio_19.tif")
predictors <- stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)
names (predictors1)
plot(CzechBorders)
class(predictors)
CzechBorders <- shapefile("CZE_adm0.shp") #loading the territoty 

#crop predictord by our territory
predictors1 <- crop(predictors, CzechBorders)
predictors1 <-stack(predictors1)
class(predictors1)
plot(predictors)
plot(predictors, 1) #show plot of the first predictor

#our data (created shp file in QGIS from csv)
points<- shapefile("czhlenky.shp")
points (points, col="blue")
#this is all evidence of the presence of occurrence
Lamproderma <- data.frame(points@coords)
occ <- rep(1,nrow(Lamproderma))
occurrence <-cbind(Lamproderma, occ)
occurrence
names(occurrence)[1] <- "lon"
names(occurrence)[2] <- "lat"
#transform df into Spatial Points df
coordinates(occurrence) <- c("lon", "lat")
#VIF calculation
library(usdm)
v <- vifstep (predictors)
v
pred <- exclude(predictors, v)
pred
#model
premodel <- sdmData(fromula=occ~.,train = occurrence, predictors = pred, bg =
               list(method="gRandom", n=100))
model <- sdm(occ~.,premodel, methods=c("Bioclim","Maxent","RF"), replication = "boot",
         test.percent=20, n=5)
model
p1 <- predict(model,pred,filename="pr1.img", overwrite=TRUE)
plot(p1)
#view by index
plot(p1[[13]])
Bioclim <- ensemble(model, pred, filename = "bioclim.img",
               setting=list(id=c(1,2,3,4,5),method="weighted", stat="tss", opt=2))
plot(Bioclim)
points(occurrence,col="blue")
# Visualisation
maxentpred.df =  as.data.frame(Bioclim, xy=TRUE)
head(maxentpred.df)
ggplot()+
  geom_raster(aes(x=x,y=y,fill=bioclim),data=maxentpred.df)+
  geom_sf(fill='transparent', colour = "white", data=CZ)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "BIOCLIM")+
  theme_minimal()
ggsave("Bioclim_30_sec_7_preds.png", width=14, height=10, unit="cm", dpi=300)


RF <- ensemble(model, pred, filename = "RF.img",
               setting=list(id=c(11,12,13,14,15),method="weighted", stat="tss", opt=2))
plot(RF)
points(occurrence,col="blue")
library(ggplot2)
enseb1.df = as.data.frame(RF, xy=TRUE)
head(enseb1.df)
library(sf)
CZ = st_as_sf(CzechBorders) # convert to simple features object
ggplot()+
  geom_raster(aes(x=x,y=y, fill=RF),data=enseb1.df)+
  geom_sf(fill='transparent', colour = "white", data=CZ)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "RF")+
  theme_minimal()
ggsave("RF_30_sec_7_preds.jpg", width=14, height=10, unit="cm", dpi=300)
#Naxent visual
maxent <- ensemble(model, pred, filename = "maxent.img",
               setting=list(id=c(11,12,13,14,15),method="weighted", stat="tss", opt=2))
plot(maxent)
points(occurrence,col="blue")
library(ggplot2)
enseb2.df = as.data.frame(maxent, xy=TRUE)
head(enseb2.df)
library(sf)
CZ = st_as_sf(CzechBorders) # convert to simple features object
ggplot()+
  geom_raster(aes(x=x,y=y, fill=maxent),data=enseb2.df)+
  geom_sf(fill='transparent', colour = "white", data=CZ)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "Maxent")+
  theme_minimal()

ggsave("Maxent_30_sec_7_preds.jpg", width=14, height=10, unit="cm", dpi=300)
#BART
devtools::install_github('cjcarlson/embarcadero')
library(embarcadero)
library(sp)
library(raster)
library(dismo)
library(rgeos)
library(rgdal)
library(maptools)
library(rJava)
library(gbm)
install.packages('reshape')
library(reshape)
presvals <- raster::extract(predictors1, points)
backgr <- spsample (points, 300, type = "random")
backgvals <- raster::extract(predictors1, backgr)
pr <- c(rep(1, nrow(presvals)), rep(0, nrow(backgvals)))
sdmdata <- data.frame(cbind(pr, rbind(presvals, backgvals)))
View(sdmdata)
sdmdata = na.omit(sdmdata)
rm(pr,backgvals, presvals)
varimp.diag(sdmdata[,-1],
            sdmdata[,1], iter=50)

step.model <- variable.step(x.data=sdmdata[,-1],
                            y.data=sdmdata[,1])
step.model
predictors_new <- predictors1 [[c(3,7,13)]]
presvals1 <- raster::extract(predictors_new, points)
backgr <- spsample (points, 300, type = "random")
backgvals1 <- raster::extract(predictors_new, backgr)
pr1 <- c(rep(1, nrow(presvals1)), rep(0, nrow(backgvals1)))
sdmdata1 <- data.frame(cbind(pr1, rbind(presvals1, backgvals1)))
bart_m2 <- bart(y.train=sdmdata[,1],
                x.train = sdmdata[,-1],
                keeptrees = TRUE)
bart_pred2 <- predict(bart_m2, predictors_new)
plot(bart_pred2)
summary(bart_m2)

bart_pred2.df =  as.data.frame(bart_pred2, xy=TRUE)
head(bart_pred2.df)
##final visual BART
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),bart_pred2.df)+
  geom_sf(fill='transparent', colour = "white", data=CZ)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "BART")+
  theme_minimal()
ggsave("BART_30_sec_3_preds.jpg", width=14, height=10, unit="cm", dpi=300)
