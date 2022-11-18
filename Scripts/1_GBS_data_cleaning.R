##########################
#### user: Lisbeth Hordley
#### date: July 2022
#### info: Cleaning & exploring GBS data for analysis

# load packages
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(rnrfa)
library(store)

# install.packages("remotes")
# remotes::install_github("edgararuiz/connections")
# install.packages("processx")
# install.packages("rlang")
# update.packages("rlang")
# devtools::install_github("gcfrench/store")

# read in data
gbs_data <- read.csv("../Data/GBS_raw_2016_2021.csv", header=TRUE)

# extract info from date
gbs_data$date <- dmy(gbs_data$date)
gbs_data$year <- year(gbs_data$date)
gbs_data$month <- month(gbs_data$date)
gbs_data$day <- day(gbs_data$date)

# only keep years 2016 - 2021
gbs_data <- gbs_data[gbs_data$year>=2016 & gbs_data$year<=2021,]
length(unique(gbs_data$year))
length(unique(gbs_data$month))
length(unique(gbs_data$grid_reference)) # 4705

# check grid references
unique(nchar(gbs_data$grid_reference)) # varying lengths: 7,8,9,10,11,12,14 
# odd numbers are wrong - these are Irish GRs which only have one letter
gbs_data$grid_reference <- ifelse(nchar(gbs_data$grid_reference) %% 2==0, gbs_data$grid_reference, NA)  # if number of characters is odd - change to NA
gbs_data <- gbs_data[!is.na(gbs_data$grid_reference),]

# now only have even grid references - all UK now
# 14 characters is odd - maximum is 12 (10 figures = 1x1m square) - remove these
gbs_data <- gbs_data[!nchar(as.character(gbs_data$grid_reference)) >= 14, ]

# convert these to lat/lon to see whether they are correct grid references
gridrefs <- data.frame(unique(gbs_data$grid_reference)) ## 4665 gridrefs
colnames(gridrefs) <- "grid_reference"
lon_lat <- as.data.frame(osg_parse(grid_refs=gridrefs$grid_reference, coord_system = "WGS84"))
# doesn't work first time round - some gridrefs must be incorrect
gridrefs <- data.frame(gridrefs[!grepl('^[L]', gridrefs$grid_reference),]) # incorrect - location is Grange-over-sands
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[I]', gridrefs$grid_reference),]) # incorrect [two of these I think the 'I' is a mistake - works without it - but in Ireland so don't want them anyway] other one should be TR
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[W]', gridrefs$grid_reference),]) # jersey
colnames(gridrefs) <- "grid_reference"
# try lat/lon again
lon_lat <- as.data.frame(osg_parse(grid_refs=gridrefs$grid_reference, coord_system = "WGS84"))
# still not working - issue with lots of zeros e.g. SV000001
gridrefs <- data.frame(gridrefs[!(gridrefs$grid_reference=="SV0020100008" | gridrefs$grid_reference=="SV000001" |
                         gridrefs$grid_reference=="SV0002300004" | gridrefs$grid_reference=="SV00000000" |
                         gridrefs$grid_reference=="SV000000" | gridrefs$grid_reference=="SV0000000000"),])
colnames(gridrefs) <- "grid_reference"
# try lat/lon again
lon_lat <- as.data.frame(osg_parse(grid_refs=gridrefs$grid_reference, coord_system = "WGS84"))
# works! # 4654 grid refs (11 removed)
# merge back into main dataframe
gridrefs <- cbind(gridrefs, lon_lat)

# some lat/lon points are way off
gridrefs <-  data.frame(gridrefs[-grep('^[O]', gridrefs$grid_reference),]) # incorrect
colnames(gridrefs) <- c("grid_reference", "lon","lat")
gridrefs <-  data.frame(gridrefs[-grep('TO', gridrefs$grid_reference),]) # incorrect
colnames(gridrefs) <- c("grid_reference", "lon","lat")
gridrefs <-  data.frame(gridrefs[-grep('TN', gridrefs$grid_reference),]) # incorrect
colnames(gridrefs) <- c("grid_reference", "lon","lat")
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "TM958269",]) # incorrect
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "NP518653",]) # incorrect
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "SY309304",]) # incorrect
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "SB793485",]) # incorrect- in the middle of the sea
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "NW6316015144",]) # N ireland
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "NW0218121336",]) # N ireland
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "NV4000505690",]) # N ireland
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "NW5087827122",]) # N ireland
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "SZ62955889",]) # In the sea south of IoW
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "SY66476589",]) # In the sea south of Weymouth
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "SY789553",]) # In the sea south of Weymouth
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "NR088648",]) # In the sea west of Islay
# all done now
unique(nchar(gridrefs$grid_reference))

# Need to remove Isle of Man as we're looking at GB only
# get easting and northing from grid references
east_north <- as.data.frame(osg_parse(grid_refs=gridrefs$grid_reference))
gridrefs <- cbind(gridrefs, east_north)

ggplot(gridrefs, aes(easting, northing))+
  geom_point()+coord_fixed()

## Need to remove Isle of Man records (GB only) (already removed Northern Ireland by gridref above and no Channel Island records)
# Filter out Isle of Man records ----
sqIM <- gridrefs[gridrefs$northing < 505000 & gridrefs$northing > 450000 &
                   gridrefs$easting > 200000 & gridrefs$easting < 250000,]
ggplot(gridrefs, aes(easting, northing))+
  geom_point()+coord_fixed()+
  geom_point(aes(easting, northing), data = sqIM, color="red") # looks good

ggplot(gridrefs[!gridrefs$grid_reference %in% sqIM$grid_reference,], aes(easting, northing))+
  geom_point()+coord_fixed()

## This removes the isle of man hectads
gridrefs <- gridrefs[!gridrefs$grid_reference %in% sqIM$grid_reference,] ## 4627 sites (6 sites are IoM)

# Now merge in gridrefs so gbs_data has the correct gridreferences 
gbs_data2 <- merge(gbs_data, gridrefs, by="grid_reference", all.y=TRUE)
length(unique(gbs_data2$grid_reference)) # 4627

# add in species common names
names <- read.csv("../Data/Butterfly_names.csv", header=TRUE)
gbs_data2 <- merge(gbs_data2, names, by.x="species", by.y="scientific_name", all.x=TRUE)
gbs_data2$id <- NULL
gbs_data2$location <- NULL
gbs_data2 <- gbs_data2[,c(2,14:17,3,13,12,11,18,1,4,5:10)]
# save file
write.csv(gbs_data2, file="../Data/GBS_2016_2021_cleaned.csv", row.names=FALSE)
  
# ## extract points from SX99 hectad - Exeter - to match with OSMM data in QGIS
# exeter <- grid_references[grid_references$hectad=="SX99",]
# exeter$hectad <- NULL
# # change all grid references to 100x100m square (smallest resolution)
# exeter2 <- exeter %>%
#   rowwise() %>%
#   mutate(hectare = hectare(grid_reference))
# # calculate XY coordinates from middle of 100x100m square
# exeter2 <- exeter2 %>%
#   rowwise() %>%
#   mutate(northing_centre = northing(grid_reference, centre=TRUE),
#          easting_centre = easting(grid_reference, centre = TRUE))
# 
# write.csv(exeter2, file="../Data/Exeter_coordinates_test.csv", row.names=FALSE)
