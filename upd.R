rm(list=ls())

lamproderma <- read.csv(file = 'lamproderma.csv', header = TRUE)
summary(lamproderma)
attach(lamproderma)
dateOnly <- as.Date(lamproderma$eventDate)
df <- lamproderma 
detach(lamproderma)
attach(df)
#cleaning the dataset
df <- data.frame(lamproderma, dateOnly)
ncols <- max(stringr::str_count(df$eventDate, "T")) + 1
ncols
colmn <- paste("col", 1:ncols)
df <- cbind(df, stringr::str_split_fixed(df$eventDate, "T", ncols))
head(df)
df$eventDate <- NULL
df$"2" <- NULL
attach(df)
library(ggplot2)
library(dismo)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(rJava)
library(gbm)
library(raster)
library(sf)
library(tidyr)
library(usdm)

base <- ggplot(df, aes(dateOnly, species)) +
  geom_line()
base + scale_x_date(date_labels = "%b")

library(maptools)
data(wrld_simpl)
p<-plot(wrld_simpl, xlim=c(-20,160), ylim=c(30,80), axes=TRUE, col="light yellow")
# restore the box around the map
box()
# add the points
points(df$decimalLongitude, df$decimalLatitude, col='orange', pch=20, cex=0.75)


library(rgdal)
library(rgeos)

URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
fil <- basename(URL)
if (!file.exists(fil)) download.file(URL, fil)
fils <- unzip(fil)
oceans <- readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
                  stringsAsFactors=FALSE, verbose=FALSE)

test.points <- data.frame(lon = c(-5,-3,-4,19,19,19),
                          lat = c(45,45,45,50,51,52))

coordinates(test.points) <- ~lon+lat
proj4string(test.points) <- CRS(proj4string(oceans))

over(test.points, oceans)

library(ggplot2)
ggplot() +
  geom_polygon(data = rwa2, aes(x = long, y = lat, group= group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = df, aes(x = Lon, y = Lat, z = z, fill = z), alpha = 0.8) +
  ggtitle("State Data") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_distiller(type = "div", palette = "Spectral")+
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 10)) +
  coord_map()