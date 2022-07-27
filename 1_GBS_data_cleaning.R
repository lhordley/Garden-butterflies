##########################
#### user: Lisbeth Hordley
#### date: July 2022
#### info: Cleaning & exploring GBS data for analysis

rm(list = ls())

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
# gbs_data <- merge(gbs_data, gridrefs, by="grid_reference", all.y=TRUE)
# length(unique(gbs_data$grid_reference))

#
# make maps
worldmap = map_data('world')
all_recs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = gridrefs, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), colour="black", size=2, shape=20) + 
  theme_void() +
  theme(title = element_text(size = 12))
all_recs
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

## still a few more to check but nearly there 












# Change gridrefs to 10km for plotting
gridrefs <- data.frame(unique(gbs_data$grid_reference))
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[L]', gridrefs$grid_reference),]) # incorrect - location is Grange-over-sands
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[I]', gridrefs$grid_reference),]) # incorrect [two of these I think the 'I' is a mistake - works without it - in Ireland so maybe don't want them anyway] other one should be TR
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[W]', gridrefs$grid_reference),]) # jersey
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[E]', gridrefs$grid_reference),]) # incorrect - location is Wokingham, should be SU
colnames(gridrefs) <- "grid_reference"
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "ST336843138284",]) # maybe too long? 10 figures is 1x1m - what is 12 figures?!
colnames(gridrefs) <- "grid_reference"
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "TQ543636197552",]) # maybe too long?
colnames(gridrefs) <- "grid_reference"
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "SN236989255749",]) # maybe too long?
colnames(gridrefs) <- "grid_reference"
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "SE422720428371",]) # maybe too long?
colnames(gridrefs) <- "grid_reference"

grid_references <- gridrefs %>%
  select(grid_reference) %>%
  mutate(grid_reference = as_gridref(grid_reference))
#grid_references <- unique(grid_references)
grid_references <- grid_references %>%
  rowwise() %>%
  mutate(hectad = hectad(grid_reference))

# make a map - first just at the hectad level, then coloured by number of gardens within each hectad with at least one record between 2016 and 2021

uk_grid_refs <- data.frame(uk_ireland_tenkm_grid_squares)
uk_grid_refs <- data.frame(uk_grid_refs[uk_grid_refs$geographical !="sea",])
uk_grid_refs <- data.frame(uk_grid_refs[-grep('^[Ireland]', uk_grid_refs$country),])
uk_grid_refs <- na.omit(uk_grid_refs)

true_hecs <- data.frame(uk_grid_refs$ten_km)
colnames(true_hecs) <- "ten_km"

grid_references2 <- merge(grid_references, true_hecs, by.x="hectad", by.y="ten_km")
# now convert these to lat/lon
hectads <- data.frame(unique(grid_references2[,1])) ## 1377 10km gridrefs
colnames(hectads) <- "hectad"
lon_lat <- as.data.frame(osg_parse(grid_refs=hectads$hectad, coord_system = "WGS84"))

hectads <- cbind(hectads, lon_lat)
grid_references2 <- merge(grid_references2, hectads, by="hectad", all.x=TRUE)

# create summary - number of gardens surveyed within each 10km square (i.e. number of rows for each 10km square)
gbs_summary <- grid_references2 %>% group_by(hectad, lon, lat) %>% summarise(n_gardens=n())

# make maps
worldmap = map_data('world')
hectad_records <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = gbs_summary, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), colour="black", size=2, shape=0) + 
  theme_void() +
  theme(title = element_text(size = 12))
hectad_records
  
worldmap = map_data('world')
hectad_records2 <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = gbs_summary, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour=n_gardens), shape=15, size=2) + 
  scale_color_viridis_c(name="Number of gardens") + 
  theme_void() +
  theme(title = element_text(size = 12))
hectad_records2  

  
## extract points from SX99 hectad - Exeter - to match with OSMM data in QGIS
exeter <- grid_references[grid_references$hectad=="SX99",]
exeter$hectad <- NULL
# change all grid references to 100x100m square (smallest resolution)
exeter2 <- exeter %>%
  rowwise() %>%
  mutate(hectare = hectare(grid_reference))
# calculate XY coordinates from middle of 100x100m square
exeter2 <- exeter2 %>%
  rowwise() %>%
  mutate(northing_centre = northing(grid_reference, centre=TRUE),
         easting_centre = easting(grid_reference, centre = TRUE))

write.csv(exeter2, file="../Data/Exeter_coordinates_test.csv", row.names=FALSE)




#################################################################################################################

# Look at the data in more detail

# First remove incorrect gridrefs
gbs_data <-  data.frame(gbs_data[-grep('^[L]', gbs_data$grid_reference),]) # incorrect - location is Grange-over-sands
gbs_data <-  data.frame(gbs_data[-grep('^[I]', gbs_data$grid_reference),]) # incorrect [two of these I think the 'I' is a mistake - works without it - in Ireland so maybe don't want them anyway] other one should be TR
gbs_data <-  data.frame(gbs_data[-grep('^[W]', gbs_data$grid_reference),]) # jersey
gbs_data <-  data.frame(gbs_data[-grep('^[E]', gbs_data$grid_reference),]) # incorrect - location is Wokingham, should be SU

# Some GRs too long - 12 figure. Remove or change to 10 figure?
# GRs don't match closely to location given - remove for now

# remove if length == 14
gbs_data <- gbs_data[!nchar(as.character(gbs_data$grid_reference)) >= 14, ]
# now only have 10 and 6 figure GR - not sure why we don't have 8 figure?

# also remove all Irish records - remove if second character is a number 
gbs_data$grid_reference <- ifelse(nchar(gbs_data$grid_reference) %% 2==0, gbs_data$grid_reference, NA)  # if number of characters is odd - change to NA (Irish GRs missing extra letter so all are odd)
gbs_data <- gbs_data[!is.na(gbs_data$grid_reference),]

# Now change all GRs to same scale - 100m (6 figure)
# Don't do this - some gardens have the same 100m grid reference but are actually different records - maybe neighbours
# gridrefs <- gbs_data %>%
#   select(grid_reference) %>%
#   mutate(grid_reference = as_gridref(grid_reference)) # change class of grid references
# 
# gridrefs <- gridrefs %>% group_by(grid_reference) %>%
#   rowwise() %>%
#   mutate(GR_100m = hectare(grid_reference)) # calculate 100m grid reference 
# gridrefs <- unique(gridrefs)
# gbs_data <- merge(gbs_data, gridrefs, by="grid_reference", all.x=TRUE) # merge back in

# calculate species richness for each location in each year
#gbs_data$quantity <- NULL
gbs_sp_rich <- gbs_data %>% group_by(grid_reference, year, garden_direction, garden_hard_surface, garden_ivy,
                                     garden_long_grass, garden_long_grass_area, garden_size, day, month) %>% select(-quantity) %>%
                                     summarise(n_spp=n())
# This is the number of species seen on a specific date
# calculate number of species recorded in each month (regardless of year, then including year) - how much affect does month have?
gbs_sp_rich <- gbs_sp_rich[!gbs_sp_rich$year==2022,]
gbs_sp_rich <- gbs_sp_rich[!gbs_sp_rich$year<2016,]

month_summ <- gbs_sp_rich %>% group_by(month) %>% summarise(n_spp=sum(n_spp))
plot(n_spp ~ month, data=month_summ)
# big peaks in july and august

# how does this vary over years?
gbs_data <- gbs_data[!gbs_data$year==2022,]
gbs_data <- gbs_data[!gbs_data$year<2016,]
month_year_summ <- gbs_data %>% group_by(date, day, month, year) %>% select(-quantity) %>% summarise(n_spp=n_distinct(species))
# we have records for each month in each year 2016 - 2021
nspp <- ggplot(month_year_summ, aes(x=date, y=n_spp)) +
  geom_point() +
  theme_bw() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(y="Number of species", x="Date") +
  theme(axis.text = element_text(size=12))
nspp
ggsave(nspp, file="../Graphs/Summaries/N_species_2016_2021.png")
# possible annomoly on the 1/1/2016 - 10 species found which is unlikely (only from one site)

