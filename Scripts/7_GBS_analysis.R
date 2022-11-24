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
gbs_data <- read.csv("Data/Raw GBS data/GBS_raw_2016_2021.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/Land_cover_GBS_100m.csv", header=TRUE)

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

# 1. Which gardens have zeros in the columns? This is only possible for area of long grass
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
check <- garden_features3[,c("grid_reference", "garden_long_grass", "garden_long_grass_area")]
# NH562225
# SN656896 
# SO791457
# SP2574077407
# SU4037990084
# TQ121959 (negative value of area)
# SJ9213172288 (negative value of area)
# SU174653 (negative value of area)
# SS825128 (negative value of area)
# SX817557 (negative value of area)

# 3. Does area of long grass and garden size match up?

# Garden size 1 category is correct (all long grass <50)
# Garden size 2 category:
  # TQ121959 negative long grass area
  # NJ9209908413 long grass area is 300
# Garden size 3 category:
  # SS825128 negative long grass area
  # SJ9213172288 negative long grass area
  # TM253564 long grass area is 600
  # TG075067 long grass area is 3000
# Garden size 4 category: 
  # SX817557 negative long grass area
  # SU174653 negative long grass area


# Negative values in long grass area and zeros in all other columns get NAs
# These can be removed easily if needed
garden_features_final <- garden_features3 %>% 
  mutate(across(c(garden_direction, garden_hard_surface, garden_ivy, garden_long_grass, garden_size), ~ na_if(., 0)))
garden_features_final <- garden_features_final %>% mutate(garden_long_grass_area = replace(garden_long_grass_area, 
                         which(garden_long_grass_area<0), NA))




















##### Landscape and abundance
gbs_landscape[,3:24] <- sapply(gbs_landscape[,3:24],as.numeric)

gbs_abund_land <- merge(gbs_abund, gbs_landscape, by.x="grid_reference", by.y="site")
plot(gbs_abund_land$area.Urban, gbs_abund_land$SINDEX) # possible negative relationship
plot(gbs_abund_land$area.Suburban, gbs_abund_land$SINDEX) # possible negative relationship too

plot(gbs_abund_land$area.Broadleaf.woodland, gbs_abund_land$SINDEX) # positive
plot(gbs_abund_land$area.Arable.and.horticulture, gbs_abund_land$SINDEX) # positive
plot(gbs_abund_land$area.Improved.grassland, gbs_abund_land$SINDEX) # positive

library(lme4)
library(lmerTest)
mod <- lmer(log(SINDEX) ~ area.Suburban + (1|M_YEAR), data=gbs_abund_land)
summary(mod)

hist(resid(mod)) # log skewed is normally distributed
qqnorm(resid(mod))
qqline(resid(mod))

##### Landscape and richness
gbs_rich_land <- merge(gbs_richness, gbs_landscape, by.x="grid_reference", by.y="site")
plot(gbs_rich_land$area.Urban, gbs_rich_land$site_rel_SR) # possible negative relationship
plot(gbs_rich_land$area.Suburban, gbs_rich_land$site_rel_SR) # possible negative relationship too

plot(gbs_rich_land$area.Broadleaf.woodland, gbs_rich_land$site_rel_SR) # positive
plot(gbs_rich_land$area.Arable.and.horticulture, gbs_rich_land$site_rel_SR) # positive
plot(gbs_rich_land$area.Improved.grassland, gbs_rich_land$site_rel_SR) # positive

mod <- lmer(site_rel_SR ~ area.Improved.grassland + (1|year), data=gbs_rich_land)
summary(mod)

hist(resid(mod)) # looks good
qqnorm(resid(mod))
qqline(resid(mod))

