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
ggsave(all_recs, file="../Graphs/Summaries/All_GBS_recs_2016_2021.png")


# Convert to 10km grid references and plot heatmap 
grid_references <- gridrefs %>%
  select(grid_reference) %>%
  mutate(grid_reference = as_gridref(grid_reference))
grid_references <- grid_references %>%
  rowwise() %>%
  mutate(hectad = hectad(grid_reference))
# convert hectad to lat/lon
lon_lat <- as.data.frame(osg_parse(grid_refs=grid_references$hectad, coord_system = "WGS84"))
grid_references <- cbind(grid_references, lon_lat)

# create summary - number of gardens surveyed within each 10km square (i.e. number of rows for each 10km square)
gbs_summary <- grid_references %>% group_by(hectad, lon, lat) %>% summarise(n_gardens=n())

worldmap = map_data('world')
hectad_recs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = gbs_summary, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour=n_gardens), shape=20, size=1.2) + 
  scale_color_viridis_c(name="Number of gardens") + 
  theme_void() +
  theme(title = element_text(size = 12))
hectad_recs  
ggsave(hectad_recs, file="../Graphs/Summaries/Hectad_GBS_recs_2016_2021.png")

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




#################################################################################################################

# Look at the data in more detail
gbs_data <- read.csv("../Data/GBS_2016_2021_cleaned.csv", header=TRUE)

# calculate species richness for each location in each year
#gbs_data$quantity <- NULL
gbs_sp_rich <- gbs_data2 %>% group_by(grid_reference, year, garden_direction, garden_hard_surface, garden_ivy,
                                     garden_long_grass, garden_long_grass_area, garden_size, day, month) %>% select(-quantity) %>%
                                     summarise(n_spp=n())
# This is the number of species seen on a specific date

month_summ <- gbs_sp_rich %>% group_by(month) %>% summarise(n_spp=sum(n_spp))
plot(n_spp ~ month, data=month_summ)
# big peaks in july and august

# how does this vary over years?
month_year_summ <- gbs_data2 %>% group_by(date, day, month, year) %>% select(-quantity) %>% summarise(n_spp=n_distinct(species))
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



# Filter the data to ensure gardens have been sampled throughout the year

# Each garden must have:
  # At least one count between mid-Apil and mid-May
  # At least two counts between July and August 
  # At least one count between mid-September and mid-October
# ALL COUNTS MUST BE AT LEAST 10 DAYS APART FOLLOWING THE WCBS GUIDELINES # 

gbs_data <- read.csv("../Data/GBS_2016_2021_cleaned.csv", header=TRUE)

# subset by columns needed for filter
gbs2 <- gbs_data[,c("grid_reference","date","year"),]
gbs2 <- unique(gbs2)

# convert date to correct format
gbs2$date <- as.Date(gbs2$date)
gbs2$date2<-format(gbs2$date, format="%m-%d")

###########
# Filter 1: gardens must be sampled once in spring, twice in summer, and once in autumn
###########

# create function to categorise a date as spring1, spring 2, summer1 etc. based on the dates above
getSeason <- function(d) {
  spring1 <- as.Date("04-15", format="%m-%d") # mid-April
  spring2 <- as.Date("05-15", format="%m-%d") # mid-May
  summer1 <- as.Date("07-01", format="%m-%d") # beginning of July
  summer2 <- as.Date("08-31", format="%m-%d") # end of August
  autumn1 <- as.Date("09-15", format="%m-%d") # mid-September
  autumn2 <- as.Date("10-15", format="%m-%d") # mid-October
  
  spring1<-format(spring1, format="%m-%d")
  spring2<-format(spring2, format="%m-%d")
  summer1<-format(summer1, format="%m-%d")
  summer2<-format(summer2, format="%m-%d")
  autumn1<-format(autumn1, format="%m-%d")
  autumn2<-format(autumn2, format="%m-%d")
  
  ifelse (d >= spring1 & d <= spring2, "Spring",
          ifelse (d >= summer1 & d <= summer2, "Summer",
                  ifelse (d >= autumn1 & d <= autumn2, "Autumn", NA))) # rows must be between these dates
}

gbs2$season <- getSeason(gbs2$date2) # use function on GBS dates
gbs2 <- na.omit(gbs2) # remove NAs where rows do not meet above conditions (just for now - added in again later once minimum conditions are met)

# now filter to ensure that there is at least 1 visit in spring, 2 in summer and 1 in autumn for each site, each year
test <-  gbs2 %>%
  group_by(grid_reference, season, year) %>%
  filter((season == "Spring" & n() >= 1) | (season == "Summer" & n() >= 2) | (season == "Autumn" & n() >= 1))
# next filter to ensure that each site has a spring, summer and autumn visit each year
# should just be able to use & instead of | above, but doesn't work
test2 <-  test %>%
  group_by(grid_reference, year) %>%
  filter(all(c("Spring", "Summer", "Autumn") %in% season))

# ND2474066857 - shouldn't be in
# NB423324 - should be in

length(unique(test2$grid_reference)) # 827 sites 
check <- test2 %>% group_by(grid_reference, year) %>% summarise(n=n()) # just to show each site has a minimum of 4 visits

###########
# Filter 2: both summer records must be at least 10 days apart
###########

# within each season/site/year, calculate the difference in days between the date and first date
test3 <- test2 %>% group_by(grid_reference, year, season) %>%  mutate(date = as.Date(date)) %>%
  arrange(date) %>% mutate(Diff = as.integer(date - first(date))) 

# for summer season only, filter those sites so that the difference between dates is at least 10 days - this removes some sites
test4 <- test3 %>% group_by(grid_reference, year) %>% filter(season=="Summer") %>% filter(max(Diff) >= 10)

test2 <- test2[!test2$season=="Summer",] # remove summer records from test2 - use those from test4 as they have the correct filter applied
test4$Diff <- NULL
test2 <- rbind(test2, test4) # put the filtered summer records (test4) and spring and autumn records (test2) together
# again make sure each site and year has each season 
# this is because some sites won't have summer records now because they didn't meet the 10 day filter
test2 <-  test2 %>%
  group_by(grid_reference, year) %>%
  filter(all(c("Spring", "Summer", "Autumn") %in% season))
length(unique(test2$grid_reference)) # 823 FINAL 
check2 <- test2 %>% group_by(grid_reference, year) %>% summarise(n=n()) # just to show each site has a minimum of 4 visits

lat_lon <- unique(gbs_data[,c("grid_reference","lat","lon"),]) # lat/lon from all GBS sites
test2_sites <- data.frame(unique(test2$grid_reference)) # grid references from filtered 823 sites
colnames(test2_sites) <- "grid_reference"
lat_lon_filtered <- merge(lat_lon, test2_sites, by="grid_reference", all.y=TRUE) # filter so the filtered sites have lat/lon to plot


worldmap = map_data('world')
filtered_sites <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = lat_lon_filtered, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), shape=20, size=2) + 
  theme_void() +
  theme(title = element_text(size = 12))
filtered_sites
ggsave(filtered_sites, file="../Graphs/Summaries/Filtered_GBS_sites_2016_2021.png")
# not bad coverage 

# filter back into gbs data
test2 <- unique(test2[,c("grid_reference", "year")]) # only need grid reference and year - date doesn't matter as we've already done the filter on date
gbs_filter_final <- merge(gbs_data, test2, by=c("grid_reference", "year"), all.y=TRUE) # merge back into GBS data
# this filters the right sites and years in - but keeps ALL dates that those sites have been visited in those years that meet the criteria
length(unique(gbs_filter_final$grid_reference)) # 823

check <- gbs_filter_final %>% group_by(grid_reference, year) %>% summarise(n_dates=n_distinct(date)) # just to show each site has a minimum of 4 visits
gbs_filter_final$day_month <- NULL
gbs_filter_final <- gbs_filter_final[,c(1,3:9,2,10:18)]
write.csv(gbs_filter_final, file="../Data/GBS_2016_2021_cleaned_filtered.csv", row.names=FALSE)


## Some more plots using filtered data
rm(list = ls())

gbs_data <- read.csv("../Data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
# Number of records for each species - each year and total

sp_summ_year <- gbs_data %>% group_by(species, year) %>% summarise(n_records=n())
sp_summ <- gbs_data %>% group_by(species) %>% summarise(n_records=n())
# quite a few species with 1 or 2 records over the whole time period
# most species with less than e.g. 4 records per year won't be included 

# how many times has a site recorded butterflies in each year?
gbs_data$date <- as.Date(gbs_data$date)
check <- gbs_data %>% group_by(grid_reference, year) %>% summarise(n=n_distinct(date)) # min 5, max 148
# sites that meet the minimum criteria of 4 sites happen to also have an extra visit in there so min is 5

# Maybe need an additional filter that species must be recorded at least 4 times a year, each year to be included? 
# Maybe also a filter on the number of sites a species is found at? 

sp_summ_site <- gbs_data %>% group_by(species) %>% summarise(n_sites=n_distinct(grid_reference))
sp_summ_site_year <- gbs_data %>% group_by(species,year) %>% summarise(n_sites=n_distinct(grid_reference))


mean_abund <- gbs_data %>% group_by(date, species) %>% summarise(abund=mean(quantity))
mean_abund$date <- as.Date(mean_abund$date, format = "%Y/%m/%d")

# extract day-month and year separately to plot
mean_abund$year <- format(as.Date(mean_abund$date),'%Y')
mean_abund$monthday <- format(as.Date(mean_abund$date), "%m-%d")
mean_abund$month <- format(as.Date(mean_abund$date), "%m")

mean_abund <- mean_abund %>%
  mutate(monthday2 = paste0("2020-", monthday))
str(mean_abund)
mean_abund$monthday2 <- as.Date(mean_abund$monthday2, format = "%Y-%m-%d")

# save all these plots for each species

species <- unique(mean_abund$species)

for(i in species) {
    print(i)
# Printing ggplot within for-loop
temp_plot <- ggplot(data=mean_abund[mean_abund$species==i,], mapping=aes(x=monthday2, y=abund, colour=year))+
  geom_point()+
  facet_grid(facets = year ~ ., margins = FALSE)+
  scale_x_date(date_breaks="1 month", date_labels="%m")+
  labs(x="Month", y="Mean abundance")+
  ggtitle(i) +
  theme_bw()
ggsave(temp_plot, file=paste0("../Graphs/Species raw data/Mean_daily_abundance_2016_2021_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}
# calculate the maximum weekly abundance value for each site/species
gbs_data$week <- lubridate::week(ymd(gbs_data$date))

max_weekly <- gbs_data %>% group_by(grid_reference, species, week, year) %>% summarise(max_abund=max(quantity))

# take the mean weekly abundance value for each species across all sites
mean_max_weekly <- max_weekly %>% group_by(species, week, year) %>% summarise(mean_abund=mean(max_abund))
mean_max_weekly$year <- as.factor(mean_max_weekly$year)

# plot this for all species
for(i in species) {
  print(i)
  # Printing ggplot within for-loop
  temp_plot <- ggplot(data=mean_max_weekly[mean_max_weekly$species==i,], mapping=aes(x=week, y=mean_abund, colour=year))+
    geom_point()+
    facet_grid(facets = year ~ ., margins = FALSE)+
    labs(x="Week", y="Mean abundance")+
    ggtitle(i) +
    theme_bw()
  ggsave(temp_plot, file=paste0("../Graphs/Species raw data/Max_mean_weekly_abundance_2016_2021_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}



