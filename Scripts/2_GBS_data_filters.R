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
gbs_data <- read.csv("Data/GBS_2016_2021_cleaned.csv", header=TRUE)


# Filter the data to ensure gardens have been sampled throughout the year

# Each garden must have:
# At least one count between mid-Apil and mid-May
# At least two counts between July and August 
# At least one count between mid-September and mid-October
# ALL COUNTS MUST BE AT LEAST 10 DAYS APART FOLLOWING THE WCBS GUIDELINES # 

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

# Use these gridrefs to check this worked:
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


###########
# Filter 3: remove rare species
###########

# Remove species that are recorded in less than 1% of ALL sites (use cleaned data, but not filtered)
gbs_raw <- read.csv("../Data/GBS_2016_2021_cleaned.csv", header=TRUE)

# number of sites each species has been recorded at across years
species_sites <- gbs_raw %>% group_by(common_name, species) %>% summarise(n_sites=n_distinct(grid_reference))
species_sites$total_sites <- length(unique(gbs_raw$grid_reference))
species_sites$percentage <- (species_sites$n_sites/species_sites$total_sites)*100
write.csv(species_sites, file="../Data/GBS_species_sites.csv", row.names=FALSE)
# this will remove 26 species
species_sites <- species_sites[species_sites$percentage>=1,] # from 57 to 31 species
new_species <- species_sites$common_name
gbs_filter_final <- gbs_filter_final[gbs_filter_final$common_name %in% new_species,]
length(unique(gbs_filter_final$common_name)) # 31
length(unique(gbs_filter_final$grid_reference)) # 823

write.csv(gbs_filter_final, file="../Data/GBS_2016_2021_cleaned_filtered.csv", row.names=FALSE)

