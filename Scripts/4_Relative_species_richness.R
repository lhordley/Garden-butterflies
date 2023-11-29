##########################
#### user: Lisbeth Hordley
#### date: July 2022
#### info: Calculate relative species richness

rm(list = ls())
# load packages
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sp)
library(rgdal)
library(geosphere)
options(scipen=999)

####################### Richness calculation for all sites and species ###################################

gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_raw <-read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)

# save species list
sp_list <- unique(gbs_data[,c("common_name", "species")])
write.csv(sp_list, file="Data/Raw GBS data/GBS_species_list.csv", row.names=FALSE)

sp_list_grass <- read.csv("Data/Raw GBS data/GBS_species_list_grass.csv", header=TRUE)

# generate site vector
filtered_sites <- unique(gbs_data$grid_reference) # these are the cleaned final sites, n=823
year <- unique(gbs_data$year) # 6 years

# # Subset to grassland species only
# sp_list_grass <- sp_list_grass[sp_list_grass$grass=="N",]
# gbs_data <- gbs_data[gbs_data$common_name %in% sp_list_grass$common_name,]
# length(unique(gbs_data$common_name)) # 10
# 
# gbs_raw <- gbs_raw[gbs_raw$common_name %in% sp_list_grass$common_name,]
# length(unique(gbs_raw$common_name)) # 10

##############
## Method 1 ##
##############

# Relative species richness is the proportion of species within a focal site compared to the 
# regional species richness in the sites within 100km of the focal site

# Concern was that isolated sites wouldn't have enough regional sites to represent regional richness
# Tested different distances to ensure regional richness was represented by at least 10 sites
# This distance is 140km 

relative_SR <- NULL
for (x in filtered_sites){  print(x)
  for(i in year){print(i)
  # subset data by site and year of interest
  site <- gbs_data[which(gbs_data$grid_reference == x & gbs_data$year==i), ] 
  if(dim(site)[1]==0) { 
    next 
    }
  # only take spatial data of site now (don't need SR yet)
  site_info <- site %>% distinct(grid_reference, .keep_all = TRUE) # unique doesn't always work for some reason
  # pick out details of focal site
  candidates <- gbs_raw[which(gbs_raw$grid_reference != x), ] # pick out details of all others - from cleaned data, not filtered
  candidates_info <- unique(candidates[,c("grid_reference", "lat_centre","lon_centre")]) # only need spatial data again
  
  distances <- distm(x = site_info[, c('lon_centre', 'lat_centre')], 
                  y = candidates_info[, c('lon_centre', 'lat_centre')],
                  fun = distHaversine)
  distances = t(distances)
  candidates_info$distance <- distances
  
  # pick out those that are within 100km (100000m) of the focal site
  closest <- candidates_info[candidates_info$distance<=140000,]
  
  site_SR <- length(unique(site$species))       # calculate species richness for focal site and year
  
  reg_recs <- gbs_raw[which(gbs_raw$grid_reference %in% closest$grid_reference), ]  # pull out region records
  #site$year <- NULL
  reg_recs <- rbind(reg_recs,site)    # add in focal site records (they're part of the regional richness too)
  reg_SR <- length(unique(reg_recs$species))
  
  site_rel_SR <- site_SR/reg_SR
  no_reg_sites <- length(unique(closest$grid_reference))
  out <- cbind(x,i,site_SR,reg_SR,site_rel_SR,no_reg_sites)
  relative_SR <- data.frame(rbind(relative_SR,out))
  
  }

}
relative_SR$site_rel_SR <- as.numeric(relative_SR$site_rel_SR)
relative_SR$no_reg_sites <- as.numeric(relative_SR$no_reg_sites)
relative_SR$reg_SR <- as.numeric(relative_SR$reg_SR)

colnames(relative_SR) <- c("grid_reference", "year", "site_SR", "reg_SR", "site_rel_SR", "no_reg_sites")
length(unique(relative_SR$grid_reference)) # 823
# save file
write.csv(relative_SR, file="Data/Richness data/Relative_SR_GBS.csv", row.names=FALSE)
#


####################### Abundance calculation for autumn ivy sites ###################################


## Calculate richness for species in Sept-Nov
# Just do this for all sites and filter down later on for GeNS and top 50% sites for analysis

library(geosphere)

gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered_ivy.csv", header=TRUE)
gbs_raw <-read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned.csv", header=TRUE)

# filter GBS raw to be only months sept-nov so regional richness is based on these months
gbs_raw <- subset(gbs_raw, subset = month %in% c(9,10,11))

gbs_sites <- unique(gbs_data[,c("grid_reference", "year", "species", "lon_centre", "lat_centre")]) 
gbs_raw <- unique(gbs_raw[,c("grid_reference", "species", "lon_centre", "lat_centre")])
# removes day/month info - only need to know which species were recorded at each site in each year
# for gbs_raw we don't need to know any time info - regional richness is based on entire time period

# generate site vector
filtered_sites <- unique(gbs_data$grid_reference) # these are the cleaned final sites, n=823
year <- unique(gbs_sites$year) # 6 years

relative_SR <- NULL
for (x in filtered_sites){  print(x)
  for(i in year){print(i)
    # subset data by site and year of interest
    site <- gbs_sites[which(gbs_sites$grid_reference == x & gbs_sites$year==i), ] 
    if(dim(site)[1]==0) { 
      next 
    }
    # only take spatial data of site now (don't need SR yet)
    site_info <- site %>% distinct(grid_reference, .keep_all = TRUE) # unique doesn't always work for some reason
    # pick out details of focal site
    candidates <- gbs_raw[which(gbs_raw$grid_reference != x), ] # pick out details of all others - from cleaned data, not filtered
    candidates_info <- unique(candidates[,c("grid_reference", "lat_centre","lon_centre")]) # only need spatial data again
    
    distances <- distm(x = site_info[, c('lon_centre', 'lat_centre')], 
                       y = candidates_info[, c('lon_centre', 'lat_centre')],
                       fun = distHaversine)
    distances = t(distances)
    candidates_info$distance <- distances
    
    # pick out those that are within 170km (100000m) of the focal site
    closest <- candidates_info[candidates_info$distance<=170000,]
    
    site_SR <- length(unique(site$species))       # calculate species richness for focal site and year
    
    reg_recs <- gbs_raw[which(gbs_raw$grid_reference %in% closest$grid_reference), ]  # pull out region records
    site$year <- NULL
    reg_recs <- rbind(reg_recs,site)    # add in focal site records (they're part of the regional richness too)
    reg_SR <- length(unique(reg_recs$species))
    
    site_rel_SR <- site_SR/reg_SR
    no_reg_sites <- length(unique(closest$grid_reference))
    out <- cbind(x,i,site_SR,reg_SR,site_rel_SR,no_reg_sites)
    relative_SR <- data.frame(rbind(relative_SR,out))
    
  }
  
}
relative_SR$site_rel_SR <- as.numeric(relative_SR$site_rel_SR)
relative_SR$no_reg_sites <- as.numeric(relative_SR$no_reg_sites)
relative_SR$site_SR <- as.numeric(relative_SR$site_SR)
relative_SR$reg_SR <- as.numeric(relative_SR$reg_SR)

colnames(relative_SR) <- c("grid_reference", "year", "site_SR", "reg_SR", "site_rel_SR", "no_reg_sites")
length(unique(relative_SR$grid_reference)) # 753
# save file
write.csv(relative_SR, file="Data/Richness data/Relative_SR_GBS_ivy.csv", row.names=FALSE)
#
