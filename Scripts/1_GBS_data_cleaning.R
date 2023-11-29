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

# read in data
gbs_data <- read.csv("Data/Raw GBS data/GBS_raw_2016_2021.csv", header=TRUE)
names <- read.csv("Data/Raw GBS data/Butterfly_names.csv", header=TRUE)
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
# 14 characters is odd - maximum is 12 (= 1x1m square) - remove these
# regardless these 4 only have one day of recording each so would be removed anyway
gbs_data <- gbs_data[!nchar(as.character(gbs_data$grid_reference)) >= 14, ]

# convert these to easting and northing to see whether they are correct grid references
gridrefs <- data.frame(unique(gbs_data$grid_reference)) ## 4663 gridrefs
colnames(gridrefs) <- "grid_reference"
lat_lon <- as.data.frame(osg_parse(grid_refs=gridrefs$grid_reference, coord_system = "WGS84"))
# doesn't work first time round - some gridrefs must be incorrect
gridrefs <- data.frame(gridrefs[!grepl('^[L]', gridrefs$grid_reference),]) # incorrect - location is Grange-over-sands
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[I]', gridrefs$grid_reference),]) # incorrect [two of these I think the 'I' is a mistake - works without it - but in Ireland so don't want them anyway] other one should be TR
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[W]', gridrefs$grid_reference),]) # jersey
colnames(gridrefs) <- "grid_reference"
# try lat/lon again
lat_lon <- as.data.frame(osg_parse(grid_refs=gridrefs$grid_reference, coord_system = "WGS84"))
# works
gridrefs <- cbind(gridrefs, lat_lon)
# try plotting on a map
UK <- map_data(map = "world", region = "UK") # changed map to "world"
ggplot() + 
  geom_polygon(data = UK, aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  geom_point(data = gridrefs, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), shape=20, size=2)+
  coord_map()
# quite a few definitely not right
# remove grid refs with lots of zeros e.g. SV000001
gridrefs <- data.frame(gridrefs[!(gridrefs$grid_reference=="SV0020100008" | gridrefs$grid_reference=="SV000001" |
                         gridrefs$grid_reference=="SV0002300004" | gridrefs$grid_reference=="SV00000000" |
                         gridrefs$grid_reference=="SV000000" | gridrefs$grid_reference=="SV0000000000"),])
# plot on a map again
ggplot() + 
  geom_polygon(data = UK, aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  geom_point(data = gridrefs, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), shape=20, size=2)+
  coord_map()
# still lots wrong
# remove more incorrect ones
gridrefs <-  data.frame(gridrefs[-grep('^[O]', gridrefs$grid_reference),]) # incorrect
gridrefs <-  data.frame(gridrefs[-grep('TO', gridrefs$grid_reference),]) # incorrect
gridrefs <-  data.frame(gridrefs[-grep('TN', gridrefs$grid_reference),]) # incorrect
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
# plot map again
ggplot() + 
  geom_polygon(data = UK, aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  geom_point(data = gridrefs, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), shape=20, size=2)+
  coord_map() # looks right apart from ones on Isle of Man that need removed

# 4633 grid references
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

# plot map again
ggplot() + 
  geom_polygon(data = UK, aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  geom_point(data = gridrefs, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), shape=20, size=2)+
  coord_map() # looks good

## Now need to make sure that each lat/lon and easting/northing value are in the middle of each grid square
# This is so that when we calculate the landscape buffer we more accurately capture the landscape 
# especially for sites that are 100x100m resolution with the 100m landscape buffer
# This will depend on whether the square is 1x1m (12 figure), 10x10m (10 figure), or 100x100m (8 figure)
# First remove lat/lon - this is easiest to change from easting and northing
gridrefs <- subset(gridrefs, select = -c(lat, lon))

# 8 figure grid references need 50 adding on to each easting and northing
# 10 figure grid references need 5 adding on to each easting and northing
# 12 figure grid references need 0.5 adding on to each easting and northing

# Easting 
gridrefs <- gridrefs %>% mutate(easting_centre=case_when(
  nchar(grid_reference)==8 ~ easting + 50,
  nchar(grid_reference)==10 ~ easting + 5,
  nchar(grid_reference)==12 ~ easting + 0.5,
))
# Northing
gridrefs <- gridrefs %>% mutate(northing_centre=case_when(
  nchar(grid_reference)==8 ~ northing + 50,
  nchar(grid_reference)==10 ~ northing + 5,
  nchar(grid_reference)==12 ~ northing + 0.5,
))

# Then recalculate lat/lon based on those easting and northing values
wgs84 = "+init=epsg:4326"
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
+ellps=airy +datum=OSGB36 +units=m +no_defs'

ConvertCoordinates <- function(easting,northing) {
  out = cbind(easting,northing)
  mask = !is.na(easting)
  sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS(bng)),sp::CRS(wgs84))
  out[mask,]=sp@coords
  out
}

x <- ConvertCoordinates(gridrefs$easting_centre, gridrefs$northing_centre)
colnames(x) <- c("lon_centre","lat_centre")
gridrefs <- cbind(gridrefs, x) 

# remove old easting and northing
gridrefs <- subset(gridrefs, select = -c(easting, northing))

# Now merge in gridrefs so gbs_data has the correct gridreferences 
gbs_data2 <- merge(gbs_data, gridrefs, by="grid_reference", all.y=TRUE)
length(unique(gbs_data2$grid_reference)) # 4627

# add in species common names
names <- read.csv("Data/Raw GBS data/Butterfly_names.csv", header=TRUE)
gbs_data2 <- merge(gbs_data2, names, by.x="species", by.y="scientific_name", all.x=TRUE)
gbs_data2$id <- NULL
gbs_data2$location <- NULL
gbs_data2 <- gbs_data2[,c(2,14:17,3,13,12,11,18,1,4:10)]
# save file
write.csv(gbs_data2, file="Data/Raw GBS data/GBS_2016_2021_cleaned.csv", row.names=FALSE)
  