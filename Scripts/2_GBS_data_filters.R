##########################
#### user: Lisbeth Hordley
#### date: July 2022
#### info: Filtering GBS data for analysis

# load packages
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(rnrfa)
library(store)

# read in data
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)

####################### Filters for main analysis ###################################

# Filter the data to ensure gardens have been sampled throughout the year

# Each garden must have:
# At least one count between mid-Apil and mid-May
# At least two counts between July and August 
# At least one count between mid-September and mid-October
# ALL COUNTS MUST BE AT LEAST 10 DAYS APART FOLLOWING WCBS GUIDELINES # 

# subset by columns needed for filter
gbs2 <- gbs_data[,c("grid_reference","date","year"),]
gbs2 <- unique(gbs2)

# convert date to correct format
gbs2$date <- as.Date(gbs2$date)
gbs2$date2<-format(gbs2$date, format="%m-%d")

# How many individual recording days does each garden have? 
rec_days <- gbs2 %>% group_by(grid_reference, year) %>% summarise(n_days=n())

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
  filter((season == "Spring" & n() >= 2) | (season == "Summer" & n() >= 4) | (season == "Autumn" & n() >= 2))
# next filter to ensure that each site has a spring, summer and autumn visit each year
# should just be able to use & instead of | above, but doesn't work
test2 <-  test %>%
  group_by(grid_reference, year) %>%
  filter(all(c("Spring", "Summer", "Autumn") %in% season))

# Use these gridrefs to check this worked:
# ND2474066857 - shouldn't be in
# NB423324 - should be in

length(unique(test2$grid_reference)) # 827 sites (617 with new filters)
check <- test2 %>% group_by(grid_reference, year) %>% summarise(n=n()) # just to show each site has a minimum of 4 visits (8 with new filters)

###########
# Filter 2: records within a season must be at least 10 days apart
###########

# within each season/site/year, calculate the difference in days between the date and first date
test3 <- test2 %>% group_by(grid_reference, year, season) %>%  mutate(date = as.Date(date)) %>%
  arrange(date) %>% mutate(Diff = as.integer(date - first(date))) 

# for summer season only, filter those sites so that the difference between dates is at least 10 days - this removes some sites
test4 <- test3 %>% group_by(grid_reference, year, season) %>% filter(max(Diff) >= 10)

# test2 <- test2[!test2$season=="Summer",] # remove summer records from test2 - use those from test4 as they have the correct filter applied
# test4$Diff <- NULL
# test2 <- rbind(test2, test4) # put the filtered summer records (test4) and spring and autumn records (test2) together

# again make sure each site and year has each season 
# this is because some sites won't have summer records now because they didn't meet the 10 day filter
test5 <-  test4 %>%
  group_by(grid_reference, year) %>%
  filter(all(c("Spring", "Summer", "Autumn") %in% season))
length(unique(test5$grid_reference)) # 823 FINAL (479 with new filters)
check2 <- test5 %>% group_by(grid_reference, year) %>% summarise(n=n()) # just to show each site has a minimum of 4 visits (8 visits with new filter)

# lat_lon <- unique(gbs_data[,c("grid_reference","lat_centre","lon_centre"),]) # lat/lon from all GBS sites
# test2_sites <- data.frame(unique(test2$grid_reference)) # grid references from filtered 823 sites
# colnames(test2_sites) <- "grid_reference"
# lat_lon_filtered <- merge(lat_lon, test2_sites, by="grid_reference", all.y=TRUE) # filter so the filtered sites have lat/lon to plot


worldmap = map_data('world')
filtered_sites <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = lat_lon_filtered, 
             aes(x = as.numeric(lon_centre), 
                 y = as.numeric(lat_centre)), shape=20, size=2) + 
  theme_void() +
  theme(title = element_text(size = 12))
filtered_sites
ggsave(filtered_sites, file="../Graphs/Filtered_GBS_sites_2016_2021.png")
# not bad coverage 

# filter back into gbs data
test5 <- unique(test5[,c("grid_reference", "year")]) # only need grid reference and year - date doesn't matter as we've already done the filter on date
gbs_filter_final <- merge(gbs_data, test5, by=c("grid_reference", "year"), all.y=TRUE) # merge back into GBS data
# this filters the right sites and years in - but keeps ALL dates that those sites have been visited in those years that meet the criteria
length(unique(gbs_filter_final$grid_reference)) # 823 (479 with new filters)

check <- gbs_filter_final %>% group_by(grid_reference, year) %>% summarise(n_dates=n_distinct(date)) # just to show each site has a minimum of 4 visits
gbs_filter_final$day_month <- NULL
gbs_filter_final <- gbs_filter_final[,c(1,3:9,2,10:18)]


###########
# Filter 3: remove rare species
###########

# Remove species that are recorded in less than 1% of ALL sites (use cleaned data, but not filtered)
gbs_raw <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)

# number of sites each species has been recorded at across years
# and the percentage compared to total GBS sites (ALL sites, not just the subset we're analysing)
species_sites <- gbs_raw %>% group_by(common_name, species) %>% summarise(n_sites=n_distinct(grid_reference))
species_sites$total_sites <- length(unique(gbs_raw$grid_reference)) # 4627 sites
species_sites$percentage <- (species_sites$n_sites/species_sites$total_sites)*100
write.csv(species_sites, file="Data/Raw GBS data/GBS_species_sites.csv", row.names=FALSE)
# this will remove 26 species
species_sites <- species_sites[species_sites$percentage>=1,] # from 57 to 31 species
new_species <- species_sites$common_name
gbs_filter_final <- gbs_filter_final[gbs_filter_final$common_name %in% new_species,]
length(unique(gbs_filter_final$common_name)) # 31
length(unique(gbs_filter_final$grid_reference)) # 823

write.csv(gbs_filter_final, file="Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", row.names=FALSE)



###### Also create data frame for sites with records in the top 50% 

gbs_filter_final <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv")
gbs2 <- gbs_filter_final[,c("grid_reference","date","year"),]
gbs2 <- unique(gbs2)

# convert date to correct format
gbs2$date <- as.Date(gbs2$date)
gbs2$date2<-format(gbs2$date, format="%m-%d")

# How many individual recording days does each garden have? 
rec_days <- gbs2 %>% group_by(grid_reference, year) %>% summarise(n_days=n()) # 5 - 181 recording days (10 - 181 with new filters)
rec_days_yr <- rec_days %>% group_by(grid_reference) %>% summarise(n_yrs=n())
hist(rec_days_yr$n_yrs)

rec_days <- rec_days[order(rec_days$n_days, decreasing = TRUE),]  
rec_days_top <- head(rec_days,0.50*nrow(rec_days)) # 41 - 181 recording days
length(unique(rec_days_top$grid_reference)) # 397 gardens
rec_days_top$n_days <- NULL
gbs_filter_final2 <- merge(gbs_filter_final, rec_days_top, by=c("grid_reference", "year"), all.y=TRUE)
length(unique(gbs_filter_final2$grid_reference)) # 397 gardens
write.csv(gbs_filter_final2, file="Data/Raw GBS data/GBS_2016_2021_cleaned_filtered_top_rec.csv", row.names=FALSE)

# note: n_days are not necessarily within sampling period, could be in Jan or Nov 
# variation in number of sites with x years of sampling after filtering is v similar
# regardless of original filters, 50% cut off or 25% cut off






####################### Autumn ivy data filters ###################################




#################################################################################################################################
## Step 1: filter sites using season of flowering ivy (Sept-Nov)

# read in data
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)

# subset by columns needed for filter
gbs2 <- gbs_data[,c("grid_reference","date","year"),]
gbs2 <- unique(gbs2)

# convert date to correct format
gbs2$date <- as.Date(gbs2$date)
gbs2$date2<-format(gbs2$date, format="%m-%d")

###########
# Filter 1: gardens must be sampled at least twice between sept and nov
###########

# create function to categorise a date as spring1, spring 2, summer1 etc. based on the dates above
getSeason <- function(d) {
  sept <- as.Date("09-01", format="%m-%d") # 1st sept
  nov <- as.Date("11-30", format="%m-%d") # 30th november
  
  sept<-format(sept, format="%m-%d")
  nov<-format(nov, format="%m-%d")
  
  ifelse (d >= sept & d <= nov, "Ivy", NA) # rows must be between these dates
}

gbs2$season <- getSeason(gbs2$date2) # use function on GBS dates
gbs2 <- na.omit(gbs2) # remove NAs where rows do not meet above conditions (just for now - added in again later once minimum conditions are met)

# now filter to ensure that there is at least 2 visits during 'ivy' for each site, each year
test <-  gbs2 %>%
  group_by(grid_reference, season, year) %>%
  filter((season == "Ivy" & n() >= 2))

length(unique(test$grid_reference)) # 1256 sites for 2 visits of 876 for 4 visits
check <- test %>% group_by(grid_reference, year) %>% summarise(n=n()) # just to show each site has a minimum of 4 visits

###########
# Filter 2: both ivy records must be at least 30 days apart
###########

# within each season/site/year, calculate the difference in days between the date and first date
test2 <- test %>% group_by(grid_reference, year, season) %>%  mutate(date = as.Date(date)) %>%
  arrange(date) %>% mutate(Diff = as.integer(date - first(date))) 

# for summer season only, filter those sites so that the difference between dates is at least 10 days - this removes some sites
test4 <- test2 %>% group_by(grid_reference, year) %>% filter(season=="Ivy") %>% filter(max(Diff) >= 30)
# NO42262196 in 2017 shouldn't be in - records are 02/09 and 06/09

length(unique(test4$grid_reference)) # 753 sites if records are at least 30 days apart
check2 <- test4 %>% group_by(grid_reference, year) %>% summarise(n=n()) # just to show each site has a minimum of 2 visits


lat_lon <- unique(gbs_data[,c("grid_reference","lat_centre","lon_centre"),]) # lat/lon from all GBS sites
test2_sites <- data.frame(unique(test4$grid_reference)) 
colnames(test2_sites) <- "grid_reference"
lat_lon_filtered <- merge(lat_lon, test2_sites, by="grid_reference", all.y=TRUE) # filter so the filtered sites have lat/lon to plot


worldmap = map_data('world')
filtered_sites <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = lat_lon_filtered, 
             aes(x = as.numeric(lon_centre), 
                 y = as.numeric(lat_centre)), shape=20, size=2) + 
  theme_void() +
  theme(title = element_text(size = 12))
filtered_sites
# not bad coverage 

# filter back into gbs data
test4 <- unique(test4[,c("grid_reference", "date")]) # only need grid reference and year - date doesn't matter as we've already done the filter on date
test4$date <- as.character(test4$date)
gbs_filter_final <- merge(gbs_data, test4, by=c("grid_reference", "date"), all.y=TRUE) # merge back into GBS data
# this filters the right sites and years in - but keeps ALL dates that those sites have been visited in those years that meet the criteria
length(unique(gbs_filter_final$grid_reference)) # 753

check <- gbs_filter_final %>% group_by(grid_reference, year) %>% summarise(n_dates=n_distinct(date)) # just to show each site has a minimum of 4 visits

gbs_filter_final <- gbs_filter_final[,c(1,3:9,2,10:18)]

###########
# Filter 3: remove rare species
###########

# Remove species that are recorded in less than 1% of ALL sites (use cleaned data, but not filtered)
gbs_raw <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)

# number of sites each species has been recorded at across years
# and the percentage compared to total GBS sites (ALL sites, not just the subset we're analysing)
species_sites <- gbs_raw %>% group_by(common_name, species) %>% summarise(n_sites=n_distinct(grid_reference))
species_sites$total_sites <- length(unique(gbs_raw$grid_reference)) # 4627 sites
species_sites$percentage <- (species_sites$n_sites/species_sites$total_sites)*100

species_sites <- species_sites[species_sites$percentage>=1,] # from 57 to 31 species
new_species <- species_sites$common_name
gbs_filter_final <- gbs_filter_final[gbs_filter_final$common_name %in% new_species,]
length(unique(gbs_filter_final$common_name)) # 26
length(unique(gbs_filter_final$grid_reference)) # 753

write.csv(gbs_filter_final, file="Data/Raw GBS data/GBS_2016_2021_cleaned_filtered_ivy.csv", row.names=FALSE)
# no observed/visits in the model like a recording effort covariate # 

#### Create dataframe with high sampling effort sites (top 50% of recording days within filter period)
gbs_filter_final <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv")
gbs2 <- gbs_filter_final[,c("grid_reference","date","year"),]
gbs2 <- unique(gbs2)

# convert date to correct format
gbs2$date <- as.Date(gbs2$date)
gbs2$date2<-format(gbs2$date, format="%m-%d")

# How many individual recording days does each garden have? 
rec_days <- gbs2 %>% group_by(grid_reference, year) %>% summarise(n_days=n()) # 2 - 60 recording days 

rec_days <- rec_days[order(rec_days$n_days, decreasing = TRUE),]  
write.csv(rec_days, file="Data/Recording_effort_day_ivy.csv", row.names=FALSE)

rec_days_top <- head(rec_days,0.50*nrow(rec_days)) # 9 - 60 recording days
length(unique(rec_days_top$grid_reference)) # 389 gardens
write.csv(rec_days_top, file="Data/Recording_effort_day_top_rec.csv", row.names=FALSE)

