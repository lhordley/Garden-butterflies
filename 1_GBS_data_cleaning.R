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

# read in data
gbs_data <- read.csv("../Data/GBS_raw_2016_2021.csv", header=TRUE)

# extract info from date
gbs_data$date <- dmy(gbs_data$date)
gbs_data$year <- year(gbs_data$date)
gbs_data$month <- month(gbs_data$date)
gbs_data$day <- day(gbs_data$date)


library(hectad)

install.packages("remotes")
remotes::install_github("edgararuiz/connections")
install.packages("processx")
install.packages("rlang")
update.packages("rlang")
devtools::install_github("gcfrench/store")

library(store)
gridrefs <- data.frame(unique(gbs_data$grid_reference))
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[L]', gridrefs$grid_reference),]) # incorrect
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[I]', gridrefs$grid_reference),]) # incorrect
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[W]', gridrefs$grid_reference),]) # jersey
colnames(gridrefs) <- "grid_reference"
gridrefs <-  data.frame(gridrefs[-grep('^[E]', gridrefs$grid_reference),]) # incorrect
colnames(gridrefs) <- "grid_reference"
gridrefs <- data.frame(gridrefs[gridrefs$grid_reference != "ST336843138284",]) # maybe too long?
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
                 y = as.numeric(lat)), colour="black", size=3, shape=20) + 
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
                 y = as.numeric(lat), colour=n_gardens), size=3, shape=20) + 
  theme_void() +
  theme(title = element_text(size = 12))
hectad_records2  

  
