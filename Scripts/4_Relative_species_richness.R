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

gbs_data <- read.csv("../Data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
gbs_raw <-read.csv("../Data/GBS_2016_2021_cleaned.csv", header=TRUE)

# number of sites each species has been recorded at across years
species_sites <- gbs_raw %>% group_by(common_name, species) %>% summarise(n_sites=n_distinct(grid_reference))
species_sites$total_sites <- length(unique(gbs_raw$grid_reference))
species_sites$percentage <- (species_sites$n_sites/species_sites$total_sites)*100
write.csv(species_sites, file="../Data/GBS_species_sites.csv", row.names=FALSE)

gbs_sites <- unique(gbs_data[,c("grid_reference", "year", "species", "lon", "lat")]) 
gbs_raw <- unique(gbs_raw[,c("grid_reference", "species", "lon", "lat")])
# removes day/month info - only need to know which species were recorded at each site in each year
# for gbs_raw we don't need to know any time info - regional richness is based on entire time period

# generate site vector
filtered_sites <- unique(gbs_data$grid_reference) # these are the cleaned final sites, n=740
year <- unique(gbs_sites$year)

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
  candidates_info <- unique(candidates[,c("grid_reference", "lat","lon")]) # only need spatial data again
  
  distances <- distm(x = site_info[, c('lon', 'lat')], 
                  y = candidates_info[, c('lon', 'lat')],
                  fun = distHaversine)
  distances = t(distances)
  candidates_info$distance <- distances
  
  # pick out those that are within 100km (100000m) of the focal site
  closest <- candidates_info[candidates_info$distance<=140000,]
  
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

#### Try again but take the closest 100 sites and also export the maximum latitude difference

# generate site vector
filtered_sites <- unique(gbs_data$grid_reference) # these are the cleaned final sites, n=740
year <- unique(gbs_data$year)
relative_SR2 <- NULL
for (x in filtered_sites){ # take each good hectad
  print(x)
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
    candidates_info <- unique(candidates[,c("grid_reference", "lat","lon")]) # only need spatial data again
    
  
  distances <- distm(x = site_info[, c('lon', 'lat')], 
                     y = candidates_info[, c('lon', 'lat')],
                     fun = distHaversine)
  distances = t(distances)
  candidates_info$distance <- distances
  
  candidates_info <- candidates_info[order(candidates_info$distance),] # sort by distance ascending
  closest <- candidates_info[1:100,] # select out closest 100
  
  site_SR <- length(unique(site$species))       # calculate species richness for focal site and year
  
  reg_recs <- gbs_raw[which(gbs_raw$grid_reference %in% closest$grid_reference), ]  # pull out region records
  site$year <- NULL
  reg_recs <- rbind(reg_recs,site)    # add in focal site records (they're part of the regional richness too)
  reg_SR <- length(unique(reg_recs$species))
  
  site_rel_SR <- site_SR/reg_SR
  lat_diff <- site_info$lat-closest$lat
  max_lat_diff <- lat_diff[which.max(abs(lat_diff))]
  out <- cbind(x,i,site_SR,reg_SR,site_rel_SR,max_lat_diff)
  relative_SR2 <- data.frame(rbind(relative_SR2,out))
  
  }
  
}


relative_SR2$site_rel_SR <- as.numeric(relative_SR2$site_rel_SR)
relative_SR2$max_lat_diff <- as.numeric(relative_SR2$max_lat_diff)
cor.test(relative_SR$site_rel_SR, relative_SR2$site_rel_SR) # 0.86

plot(relative_SR$site_rel_SR, relative_SR2$site_rel_SR)
# high correlation
# but as the relative SR gets higher, the difference between the two methods is also higher


hist(relative_SR2$max_lat_diff) # most sites have small differences in latitude
hist(relative_SR$no_reg_sites) # this is very varied up to ~1500 sites, quite a few sites use <500 candidates 
# main issue with above method is that some sites near the coast won't have as many sites to choose from - probably not an issue tho
# as it's still able to capture regional richness


# calculate species richness in each 100km square and plot this on a map 
# just out of interest to see how much it does vary over space

gbs_data$myriad <- substr(gbs_data$grid_reference, start = 1, stop = 2)
gbs_raw$myriad <- substr(gbs_raw$grid_reference, start = 1, stop = 2)

myriad_richness <- gbs_data %>% group_by(myriad) %>% mutate(nspp=n_distinct(species))
myriad_richness <- unique(myriad_richness[,c("myriad","nspp")])
myriad_richness2 <- gbs_raw %>% group_by(myriad) %>% mutate(nspp=n_distinct(species))
myriad_richness2 <- unique(myriad_richness2[,c("myriad","nspp")])

grid_refs <- myriad_richness2$myriad
# try lat/lon again
lon_lat <- as.data.frame(osg_parse(grid_refs=grid_refs, coord_system = "WGS84"))
myriad_richness2 <- cbind(myriad_richness2, lon_lat)

worldmap = map_data('world')
myriad_richness <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = myriad_richness2, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour=nspp), size=12, shape=15) + 
  scale_color_viridis_c(name="Species richness") + 
  theme_void() +
  theme(title = element_text(size = 12))
myriad_richness
