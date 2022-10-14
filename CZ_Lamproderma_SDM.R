#GBIF phenology visualisation
data <- read.csv("GBIFEurope.csv")
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
