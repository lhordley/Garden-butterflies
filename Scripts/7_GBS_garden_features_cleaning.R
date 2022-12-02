##########################
#### user: Lisbeth Hordley
#### date: November 2022
#### info: Initial GBS abundance/richness and landscape analysis
options(scipen = 100)

# libraries 
library(ggplot2)
library(dplyr)

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)

## Abundance
# put garden features and landscape parameters into same dataframe
garden_features <- garden_features[,c("grid_reference", "date", "year", "garden_direction", "garden_hard_surface",
                                      "garden_ivy", "garden_long_grass", "garden_long_grass_area",
                                      "garden_size")]
# this takes the unique of the above rows MINUS date - this should = 823 but its 891
garden_features2 <- garden_features[!duplicated(garden_features[,c("grid_reference", "garden_direction", "garden_hard_surface",
                                                                   "garden_ivy", "garden_long_grass", "garden_long_grass_area",
                                                                   "garden_size")]),]
# issues with gardens within the same year have different features
# not so much of an issue with e.g. area of long grass as could change through the year - take the mean?
# but it's an issue where gardens have different sizes within one year

duplicated_sites <- garden_features2 %>%
  group_by(grid_reference) %>% 
  filter(n() > 1) %>%
  ungroup() 
length(unique(duplicated_sites$grid_reference)) # 57 out of 823 sites with duplicated values

# NH520239 x 2: area long grass and garden size (within 2016) - EXCLUDE
# NJ561196 x 2: area long grass (within 2017) - KEEP
# NO6119107748 x 2: area long grass and garden size (within 2020) - EXCLUDE
# NZ226279 x 2: garden size (within 2016) - EXCLUDE
# NZ6202820871 x 2: across years (garden ivy) - KEEP
# SD645045 x 2: across years (garden ivy) - KEEP
# SD886021 x 3: within years (2016) long grass and area of long grass, 2017 same - KEEP
# SE2520818916 x 2: across years garden ivy - KEEP
# SE281362 x 2: across years - KEEP
# SE541080 x 4: area long grass - KEEP
# SJ811843 x 3: within year garden size - EXCLUDE
# SJ9089234892 x 2: across years - KEEP
# SK1612373407 x 2: across years - KEEP
# SK245305 x 2: across years - KEEP
# SK4497114268 x 2: area long grass - KEEP
# SK4945309892 x 3: within year garden size - EXCLUDE
# SK626008 x 2: KEEP but remove row with all zeros
# SP0974394960 x 2: within year long grass KEEP
# SP817007 x 2: across years KEEP
# SP837578 x 2: within year long grass KEEP
# SS964468 x 2: across years KEEP
# ST1405907628 x 5: within years garden size - EXCLUDE
# ST3769432785 x 2: within year ivy KEEP
# ST485012 x 2: across years KEEP
# ST501486 x 2: across years KEEP
# ST63765701 x 2: within years long grass KEEP
# ST678188 x 2: across years KEEP
# ST6960849215 x 2: across years KEEP
# ST739641 x 2: across years KEEP
# SU131295 x 2: within year garden size EXCLUDE
# SU3695501940 x 2: across years KEEP
# SU4067314308 x 3: across and between long grass KEEP
# SU4388654840 x 2: within years garden size EXCLUDE
# SU743714 x 3: within years long grass KEEP (remove the row of zeros) 
# SU855492 x 2: within years garden size EXCLUDE
# SU952400 x 2: within years long grass KEEP
# SX031531 x 2: garden size EXCLUDE
# SX793451 x 2: garden size EXCLUDE
# SX869512 x 2: across years KEEP
# SZ029796 x 2: within year long grass KEEP
# TA076300 x 2: across years KEEP
# TF1043024048 x 2: within year long grass KEEP
# TF950216 x 2: within year long grass KEEP
# TL163055 x 2: across years KEEP
# TL2787944002 x 3: across years KEEP
# TL417972 x 2: across years KEEP (but check zero in garden direction)
# TL888400 x 2: KEEP (remove duplicated rows with zeros)
# TM089243 x 2: across years KEEP
# TM097812 x 2: across years KEEP
# TM127918 x 2: across years KEEP
# TM303561 x 2: across years KEEP
# TQ1470431132 x 2: within year ivy KEEP
# TQ2087051674 x 2: across years KEEP
# TQ290051 x 2: across years KEEP
# TQ3044318486 x 2: across years KEEP
# TQ7028009690 x 2: across years KEEP
# TR1380467079 x 2: within years ivy KEEP

# REMOVE 11 SITES - these are ones that have different garden directions or sizes within a year - implausible 
sites_exclude <- c("NH520239", "NO6119107748", "NZ226279", "SJ811843", "SK4945309892", "ST1405907628", "SU131295",
                   "SU4388654840", "SU855492", "SX031531", "SX793451")

garden_features2 <- garden_features2[!garden_features2$grid_reference %in% sites_exclude,]
length(unique(garden_features2$grid_reference)) # 812

## What are the categories for each column?

# Garden direction: categories, 1-8 (north, north-east, east, south-east, south, south-west, west, north-west)
# Garden hard surface: How much of your garden is covered by hard surface? Categories 1-5
# Garden ivy: Do you have flowering ivy in your garden? 1 = yes, 2 = no
# Garden long grass: 1 = yes, 2 = no
# Area of long grass: area in square meters
# Garden size: categories, 1-4 (50 sqm or less, 50-100sqm, 100-450sqm, more than 450sqm)

# Checks to make on the remaining data:

# 1. Which gardens have zeros in the columns? This is only possible for area of long grass (only if long grass == 2)
colSums(garden_features2==0)
# 24 in garden direction
# 8 in hard surface
# 9 in long grass
# 430 long grass area (these are fine)
# 8 in garden size
# These gardens can be removed if and when I need to e.g. if using garden direction in a model,
# remove the gardens that have zeros in garden direction
# Bigger issue is gardens that have all zeros - are these errors or do they have data for other dates?
zeros <- garden_features2[rowSums(garden_features2[, 4:7,9]) == 0,]
# 5 sites that have zeros for all columns (apart from area of long grass)
# Do these sites have data on other dates?
# SK626008 - yes same year (2016) just remove one row
# SO9533919864 - no, remove site
# SP7779565063 - no, remove site
# SU743714 - yes same year (2016) just remove one row
# TL888400 - yes but for different year - just remove for that year
garden_features3 <- garden_features2[!rowSums(garden_features2[, 4:7,9]) == 0,] # removes 2 sites and 5 rows
# make sure when I merge again, merge with year too

# 2. Do Y/N for long grass and area all match up? 
check <- garden_features3[,c("grid_reference", "year", "garden_long_grass", "garden_long_grass_area")]
# These grid refs have said no long grass, but have put a positive value for area
# NH562225
# SN656896 
# SO791457
# SP2574077407
# SU4037990084

# These gridrefs have negative values of area
# TQ121959 
# SJ9213172288 
# SU174653 
# SS825128 
# SX817557 

# These grdirefs have said yes to long grass, but have put a zero value in for area
# 12 gardens - these should get NAs in both columns because the long grass data is implausible, but they have usable data for 
# other features

# These gridrefs have zeros in long grass, so will be convered to NA below (they have data for other features)
# SP986078
# ST515753
# SZ799969
# TL037084



# 3. Does area of long grass and garden size match up?

# Garden size 1 category is correct (all long grass <50)
# Garden size 2 category (50-100):
# TQ121959 negative long grass area
# NJ9209908413 long grass area is 300
# Garden size 3 category (100-450):
# SS825128 negative long grass area
# SJ9213172288 negative long grass area
# TM253564 long grass area is 600
# TG075067 long grass area is 3000
# Garden size 4 category (>450): 
# SX817557 negative long grass area
# SU174653 negative long grass area


# Negative values in long grass area and zeros in all other columns get NAs
# These can be removed easily if needed
garden_features_final <- garden_features3 %>% 
  mutate(across(c(garden_direction, garden_hard_surface, garden_ivy, garden_long_grass, garden_size), ~ na_if(., 0)))
garden_features_final <- garden_features_final %>% mutate(garden_long_grass_area = replace(garden_long_grass_area, 
                                                                                           which(garden_long_grass_area<0), NA))
length(unique(garden_features_final$grid_reference)) # 810

# Put NAs where garden long grass == 2 and area is a positive value
garden_features_final$garden_long_grass <- ifelse(garden_features_final$garden_long_grass==2 & garden_features_final$garden_long_grass_area>0, NA,
                                                  garden_features_final$garden_long_grass)
garden_features_final$garden_long_grass_area <- ifelse(garden_features_final$garden_long_grass==2 & garden_features_final$garden_long_grass_area>0, NA,
                                                  garden_features_final$garden_long_grass_area)

# Put NAs where garden long grass == 1 and area is zero
garden_features_final$garden_long_grass <- ifelse(garden_features_final$garden_long_grass==1 & garden_features_final$garden_long_grass_area==0, NA,
                                                  garden_features_final$garden_long_grass)
garden_features_final$garden_long_grass_area <- ifelse(garden_features_final$garden_long_grass==1 & garden_features_final$garden_long_grass_area==0, NA,
                                                       garden_features_final$garden_long_grass_area)


## Take maximum value across duplicates of hard surface, area of long grass and take the yes of ivy if duplicated (which is the minimum)
garden_features_final$date <- NULL
garden_features_final2 <- garden_features_final %>% group_by(grid_reference, year) %>% summarise(garden_long_grass_area=max(garden_long_grass_area),
                                                                                                 garden_hard_surface=max(garden_hard_surface),
                                                                                                 garden_ivy=min(garden_ivy),
                                                                                                 garden_long_grass=min(garden_long_grass), across())
garden_features_final2 <- unique(garden_features_final2)
length(unique(garden_features_final2$grid_reference)) # 810
duplicated_sites <- garden_features_final2 %>%
  group_by(grid_reference, year) %>% 
  filter(n() > 1) %>%
  ungroup()  
# still some duplications in garden_long_grass - sites that had 0 long grass area and 2 for long grass (no) but this changed
# to a positive value for area and 1 (yes) within the same year
# solve this by taking the min of garden long grass to change them to 1 (yes) 
# We've already taken the maximum value of area so it needs to match up with a 1 
# Duplicated_sites is now zero after this change

# Still issues with gardens that have said yes to long grass (1) but zero in area (nothing to do with duplications)


summary(garden_features_final2)
library(rnrfa)
east_north <- as.data.frame(osg_parse(grid_refs=garden_features_final2$grid_reference))
garden_features_final2 <- cbind(garden_features_final2, east_north)
lon_lat <- as.data.frame(osg_parse(grid_refs=garden_features_final2$grid_reference, coord_system = "WGS84"))
garden_features_final2 <- cbind(garden_features_final2, lon_lat)

ggplot(garden_features_final2, aes(easting, northing))+
  geom_point()+coord_fixed()

## save garden features
write.csv(garden_features_final2, file="Data/Raw GBS data/GBS_garden_features_cleaned.csv", row.names=FALSE)
