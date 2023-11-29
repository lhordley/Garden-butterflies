##########################
#### user: Lisbeth Hordley
#### date: November 2022
#### info: Put datasets together for analyses 
options(scipen = 100)

# libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

# 1. All sites and all species

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day.csv", header=TRUE)

# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))
# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 823

write.csv(gbs_analysis, file="Data/GBS_analysis_final.csv", row.names=FALSE)


# 2. Grass feeding species only

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS_grass.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
gbs_species <- read.csv("Data/Raw GBS data/GBS_species_list_grass.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day.csv", header=TRUE)

# Subset to grassland species only
gbs_species <- gbs_species[gbs_species$grass=="Y",]
gbs_abund <- gbs_abund[gbs_abund$SPECIES %in% gbs_species$common_name,]
length(unique(gbs_abund$SPECIES)) # 8 (no essex skipper or wall)
# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))
# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 801

write.csv(gbs_analysis, file="Data/GBS_analysis_final_grass.csv", row.names=FALSE)


# 2. Non-grass feeding species only

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS_not_grass.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
gbs_species <- read.csv("Data/Raw GBS data/GBS_species_list_grass.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day.csv", header=TRUE)

# Subset to grassland species only
gbs_species <- gbs_species[gbs_species$grass=="N",]
gbs_abund <- gbs_abund[gbs_abund$SPECIES %in% gbs_species$common_name,]
length(unique(gbs_abund$SPECIES)) # 15
# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))
# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 823

write.csv(gbs_analysis, file="Data/GBS_analysis_final_not_grass.csv", row.names=FALSE)


# 3. Climate zone only sites

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily_GEnS.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day.csv", header=TRUE)

gbs_abund <- unique(gbs_abund)
# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))
# subset richness to match the same sites as abundance (which are already in climate zone)
sites <- unique(gbs_abund[,1])
gbs_richness <- gbs_richness[gbs_richness$grid_reference %in% sites$grid_reference,]
length(unique(gbs_richness$grid_reference)) 
# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 785

write.csv(gbs_analysis, file="Data/GBS_analysis_final_GEnS.csv", row.names=FALSE)

# 4. Top 50% recorded sites

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day_top_rec.csv", header=TRUE)

gbs_abund <- unique(gbs_abund)
# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))
# subset richness to match the same sites as abundance (which are already in climate zone)
sites <- unique(gbs_abund[,1])
gbs_richness <- gbs_richness[gbs_richness$grid_reference %in% sites$grid_reference,]
length(unique(gbs_richness$grid_reference)) 
# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.y=TRUE)

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 397

write.csv(gbs_analysis, file="Data/GBS_analysis_final_top_rec.csv", row.names=FALSE)


# 5. Ivy autumn all species

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned_ivy.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily_ivy.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS_ivy.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day_ivy.csv", header=TRUE)

# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))

# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 746

write.csv(gbs_analysis, file="Data/GBS_analysis_final_ivy.csv", row.names=FALSE)

# 6. Ivy autumn GEnS sites

garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned_ivy.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily_ivy_GEnS.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS_ivy.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day_ivy.csv", header=TRUE)

# match sites from abundance to richness - i.e. take out sites not in those groups
sites <- unique(gbs_abund[,6])
gbs_richness <- gbs_richness[gbs_richness$grid_reference %in% sites,]
length(unique(gbs_richness$grid_reference)) # 589

# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))

# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 589

write.csv(gbs_analysis, file="Data/GBS_analysis_final_ivy_GEnS.csv", row.names=FALSE)

# 7. Ivy autumn top 50% sites

garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned_ivy.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily_ivy_top_rec.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS_ivy.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day_ivy_top_rec.csv", header=TRUE)

# match sites from abundance to richness - i.e. take out sites not in those groups
sites <- unique(gbs_abund[,6])
gbs_richness <- gbs_richness[gbs_richness$grid_reference %in% sites,]
length(unique(gbs_richness$grid_reference)) # 389

# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))

# put abundance and richness data together
gbs_analysis <- merge(gbs_abund, gbs_richness, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))
# remove unneccessary columns
gbs_analysis <- gbs_analysis[,-which(names(gbs_analysis) %in% c("site_SR", "reg_SR", "no_reg_sites"))]
# add in garden features
gbs_analysis <- merge(gbs_analysis, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.y=TRUE)

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 389

write.csv(gbs_analysis, file="Data/GBS_analysis_final_ivy_top_rec.csv", row.names=FALSE)

# 10. Ivy Red Admiral and Comma combined

garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned_ivy.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily_ivy.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day_ivy.csv", header=TRUE)

# select Red Admiral from abundance data
gbs_abund <- gbs_abund[gbs_abund$SPECIES=="Comma" | gbs_abund$SPECIES=="Red Admiral",]
# Take sum of species abundance indices per site across species 
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))

# add in garden features
gbs_analysis <- merge(gbs_abund, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 732

write.csv(gbs_analysis, file="Data/GBS_analysis_final_ivy_Comma_Red_Admiral.csv", row.names=FALSE)


# 11. Ivy summer Holly Blue

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily_ivy_HollyBlue.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/GBS_landscape_buffers_new.csv", header=TRUE)
gbs_landscape_dist <- read.csv("Data/Landscape data/Woodland_grassland_patch_distance_new.csv", header=TRUE)
rec_days <- read.csv("Data/Raw GBS data/Recording_effort_day.csv", header=TRUE)

# add in garden features
gbs_analysis <- merge(gbs_abund, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# add in landscape buffers
gbs_analysis <- merge(gbs_analysis, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
# add in landscape distances
gbs_analysis <- merge(gbs_analysis, gbs_landscape_dist, by="grid_reference", all.x=TRUE)
# add in n_days
gbs_analysis <- merge(gbs_analysis, rec_days, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"))

## put in lat/lon centre values to use to check spatial autocorrelation in models
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_lat_lon <- unique(gbs_data[,c("grid_reference", "lat_centre", "lon_centre")])
gbs_analysis <- merge(gbs_analysis, gbs_lat_lon, by="grid_reference", all.x=TRUE)
length(unique(gbs_analysis$grid_reference)) # 645

write.csv(gbs_analysis, file="Data/GBS_analysis_final_HollyBlue.csv", row.names=FALSE)

