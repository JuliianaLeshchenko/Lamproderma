data <- read.csv("GBIF_phenology.csv")
# Sample data in a vector
datetime_vector <- c(data$eventDate)
# Separate date and time using 'T' as a separator
date_vector <- sapply(strsplit(datetime_vector, "T"), `[`, 1)
time_vector <- sapply(strsplit(datetime_vector, "T"), `[`, 2)
# Create a data frame with separated date and time
data$date <- date_vector
# Display the resulting data frame
print(data)
# Delete the 'datetime_column' using the subset() function
data <- subset(data, select = -eventDate)
library(ggplot2)

# Create a new column for species subgroup
data$species_group <- ifelse(data$species %in% c("Lamproderma ovoideum", "Lamproderma aeneum"), "Strictly Nivicolous",
                             ifelse(data$species %in% c("Lamproderma pseudomaculatum", "Lamproderma splendens"), "Facultatively Nivicolous",
                                    ifelse(data$species %in% c("Lamproderma arcyrioides", "Lamproderma sauteri"), "Cryophilous", NA)))

# Create the plot
fruiting_plot <- ggplot(data, aes(month, species, fill = species_group)) +
  geom_tile() +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme(
    axis.title.x = element_blank(),  # no need for x-axis title
    axis.text.y = element_text(hjust = 0))  # left-justify species names

# Make the plot circular
fruiting_plot <- fruiting_plot + coord_polar()

# Remove non-sensical y-axis details now that plot is circular
fruiting_plot <- fruiting_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank())

# Set custom colors for the rings (pink, blue, violet)
fruiting_plot <- fruiting_plot + scale_fill_manual(values = c("pink", "blue", "violet"))

# Arrange species in a single circle
fruiting_plot <- fruiting_plot + facet_grid(rows = vars(species_group))

# Print the plot
print(fruiting_plot)
ggsave("fruiting_plots.png", width=14, height=10, unit="cm", dpi=300)
# SDM preparation
bio <- raster::getData('worldclim', var='bio', res=10)
bio
names(bio)
library(dismo)
require(usdm)
library(sp)
library(rgeos)
library(rgdal)
library(sf)
EU <- readOGR(dsn = ".", layer = "BiogeoRegions2016")
EU <- spTransform(EU, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
shapeF <- readOGR(dsn = ".", layer = "facult")

preds<- crop(bio, EU)
v <- vifstep (preds)
significantPred <- exclude(preds, v)
significantPred
plot(significantPred[[1]])
points (shape1, col="blue")
sp <- data.frame(shapeF@coords)
occ <- rep(1,nrow(sp))
occurrence <-cbind(sp, occ)
names(occurrence)[1] <- "lon"
names(occurrence)[2] <- "lat"
library(sdm)
installAll()
coordinates(occurrence) <- c("lon", "lat")
d <- sdmData(fromula=occ~.,train = occurrence, predictors = significantPred, bg =
               list(method="gRandom", n=3000))
m1 <- sdm(occ~.,d, methods=c("Bioclim", 'rf', 'maxent'), replication = "boot",
         test.percent=20, n=100)
m1
gui(m1)
biof <- raster::getData(name = 'CMIP5', var = 'bio', res = 10,
                rcp = 45, model = 'IP', year = 70,
                path = '/Diploma/futuredata/')
biof_croped<- crop(biof, EU)
names(biof_croped)<- names(preds)
en_str <- ensemble(m, preds,setting=list(id=c(201:300),method="weighted", stat="tss", opt=2))

plot(en2)
cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))
plot(en1, col=cl(200))
plot(en2, col=cl(200))
en_str.df = as.data.frame(en_str, xy=TRUE)
library(ggplot2)
EU <- readOGR(dsn = ".", layer = "BiogeoRegions2016")
EU <- spTransform(EU, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
library(sf)
EU_plot = st_as_sf(EU) # convert to simple features object
library(ggplot2)
ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=en_str.df)+
  geom_sf(fill='transparent', colour = "white", data=EU_plot)+
  scale_fill_viridis_c('')+
  labs(x="", y = "",
       title = "maxent str")+
  theme_minimal()
ggsave("maxent_str.png", width=14, height=10, unit="cm", dpi=300)

gui(m)

biof1 <- raster::getData(name = 'CMIP5', var = 'bio', res = 10,
                        rcp = 85, model = 'IP', year = 70,
                        path = '/Diploma/futuredata/')
biof1_croped<- crop(biof1, EU)
names(biof1_croped)<-names(preds)
en_facul5 <- ensemble(m, biof_croped,setting=list(id=c(201:300),method="weighted", stat="tss", opt=2))
en_facul5.df = as.data.frame(en_facul5, xy=TRUE)


shape2 <- readOGR(dsn = ".", layer = "strictly")
points (shape2, col="blue")
sp <- data.frame(shape2@coords)
occ <- rep(1,nrow(sp))
occurrence <-cbind(sp, occ)
names(occurrence)[1] <- "lon"
names(occurrence)[2] <- "lat"
library(sdm)

coordinates(occurrence) <- c("lon", "lat")
d <- sdmData(fromula=occ~.,train = occurrence, predictors = significantPred, bg =
               list(method="gRandom", n=3000))
m <- sdm(occ~.,d, methods=c("Bioclim", 'rf', 'maxent'), replication = "boot",
         test.percent=20, n=100)

#QR code for posters
install.packages('qrcode')
library(qrcode)
code <- qr_code('https://github.com/JuliianaLeshchenko/Lamproderma')
plot(code)
