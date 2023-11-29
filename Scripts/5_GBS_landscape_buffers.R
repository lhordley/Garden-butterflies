##########################
#### user: Lisbeth Hordley
#### date: November 2022
#### info: Calculate landscape buffers around GBS sites

# libraries 
library(mapview)
library(sf)
library(sp)
library(rnrfa)
library(dplyr)
library(stringr)
options(scipen = 100)

# LCM 2015 vector data
# This has been clipped to 5km around each GBS site in QGIS so it can be read into R
lcm <- st_read("Data/Landscape data/LCM_GBS_1km_final_transformed.gpkg")

# get GBS site info
gbs <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
sites <- unique(gbs[,c("grid_reference", "lon_centre", "lat_centre")])
gbs2 <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered_ivy.csv", header=TRUE)
sites_2 <- unique(gbs2[,c("grid_reference", "lon_centre", "lat_centre")])
sites_final <- rbind(sites,sites_2)
sites_final <- unique(sites_final)
write.csv(sites_final, file="Data/GBS_sites.csv", row.names=FALSE)

## Step 1: Calculate proportion of land cover in buffers around GBS sites
# buffers are 100m, 250m and 500m

sites2 = st_as_sf(sites_final,coords=c("lon_centre","lat_centre"), dim = "XY", crs=4326)
st_crs(lcm)

gridref <- sites_final$grid_reference
x_final <- NULL
scales <- c(100,250,500)
# change projection to a UK based projection
lcm <- sf::st_transform(lcm, "epsg:32630") 
st_write(lcm, "Data/Landscape data/LCM_GBS_1km_final_transformed.gpkg")
sites2 <- sf::st_transform(sites2, "epsg:32630")
st_crs(lcm)
st_crs(sites2)

for(i in gridref[1:10]){print(i)
  sites_temp <- sites2[sites2$grid_reference==i,]
  
  for(f in scales){
    
    x <- st_intersection(lcm, sites_temp |> st_geometry() |> st_buffer(dist = f))
    x["area"] <- st_area(x)
    x_area <- aggregate(area ~ bhab, x, sum)
    x_area$site <- i
    x_area$scale <- f
    x_final <- rbind(x_area, x_final)
  }
}
write.csv(x_final, file="Data/Landscape data/gbs_land_100_250_500_new.csv", row.names=FALSE)
Sys.time() 

x_final2 <- read.csv("Data/Landscape data/gbs_land_100_250_500_new.csv", header=TRUE)
x_final$area <- as.numeric(x_final$area)

land_100 <- x_final[x_final$scale==100,]
land_250 <- x_final[x_final$scale==250,]
land_500 <- x_final[x_final$scale==500,]

land_100 <- reshape(land_100, idvar = c("site", "scale"), timevar = "bhab", direction = "wide")
land_250 <- reshape(land_250, idvar = c("site", "scale"), timevar = "bhab", direction = "wide")
land_500 <- reshape(land_500, idvar = c("site", "scale"), timevar = "bhab", direction = "wide")

land_100[is.na(land_100)] <- 0
land_250[is.na(land_250)] <- 0
land_500[is.na(land_500)] <- 0

land_100$total_count <- rowSums(land_100[3:19]) # create total_count column
land_250$total_count <- rowSums(land_250[3:22]) # create total_count column
land_500$total_count <- rowSums(land_500[3:23]) # create total_count column

# only interested in 4 land cover types:
# 1. Arable
# 2. Urban + suburban
# 3. Woodland (broadleaf + conifer)
# 4. Grassland (acid, improved, neutral, calcareous and heather)

land_100$Urban_100m <- land_100$area.Urban + land_100$area.Suburban
land_100$Arable_100m <- land_100$`area.Arable and horticulture`
land_100$Woodland_100m <- land_100$`area.Broadleaf woodland` + land_100$`area.Coniferous woodland`
land_100$Grassland_100m <- land_100$`area.Acid grassland` + land_100$`area.Heather grassland` +
  land_100$`area.Neutral grassland` + land_100$`area.Improved grassland` + land_100$`area.Calcareous grassland`
land_100 <- land_100[,c("site", "Urban_100m", "Arable_100m", "Woodland_100m", "Grassland_100m", "total_count")]

land_250$Urban_250m <- land_250$area.Urban + land_250$area.Suburban
land_250$Arable_250m <- land_250$`area.Arable and horticulture`
land_250$Woodland_250m <- land_250$`area.Broadleaf woodland` + land_250$`area.Coniferous woodland`
land_250$Grassland_250m <- land_250$`area.Acid grassland` + land_250$`area.Heather grassland` +
  land_250$`area.Neutral grassland` + land_250$`area.Improved grassland` + land_250$`area.Calcareous grassland`
land_250 <- land_250[,c("site", "Urban_250m", "Arable_250m", "Woodland_250m", "Grassland_250m", "total_count")]

land_500$Urban_500m <- land_500$area.Urban + land_500$area.Suburban
land_500$Arable_500m <- land_500$`area.Arable and horticulture`
land_500$Woodland_500m <- land_500$`area.Broadleaf woodland` + land_500$`area.Coniferous woodland`
land_500$Grassland_500m <- land_500$`area.Acid grassland` + land_500$`area.Heather grassland` +
  land_500$`area.Neutral grassland` + land_500$`area.Improved grassland` + land_500$`area.Calcareous grassland`
land_500 <- land_500[,c("site", "Urban_500m", "Arable_500m", "Woodland_500m", "Grassland_500m", "total_count")]

# divide all columns by total land area
habs <- colnames(land_100)[c(2:5)]

for (i in habs){
  land_100[,i] <- land_100[,i]/land_100[,"total_count"]
}

habs2 <- colnames(land_250)[c(2:5)]

for (i in habs2){
  land_250[,i] <- land_250[,i]/land_250[,"total_count"]
}

habs3 <- colnames(land_500)[c(2:5)]

for (i in habs3){
  land_500[,i] <- land_500[,i]/land_500[,"total_count"]
}

land_100$total_count <- NULL
land_250$total_count <- NULL
land_500$total_count <- NULL

# Check for high (>0.7) correlations between urban and grassland metrics
round(cor(land_100[,c("area.Improved grassland", "area.Suburban", "area.Urban", "area.Calcareous grassland", 
                      "area.Neutral grassland", "area.Acid grassland", "area.Heather grassland",
                      "Urban_100m", "Grassland_100m")]),3)
# Urban/suburban + improved grassland -0.76
# Suburban + grassland -0.7
# Urban/suburban + grassland -0.78

# Issue is with improved grassland, other grassland types are not correlated with urban/suburban
# Urban on it's own is not an issue either, so we could have urban only + grassland in the models
# but we couldn't have suburban + grassland in 

round(cor(land_250[,c("Urban_500m", "Arable_500m", "Grassland_500m", "Woodland_500m", 
                                   "woodland_dist_m", "grassland_dist_m", "Urban_vs_Grassland_500m")]),3)

round(cor(land_500[,c("Urban_500m", "Arable_500m", "Grassland_500m", "Woodland_500m", 
                                   "woodland_dist_m", "grassland_dist_m", "Urban_vs_Grassland_500m")]),3)


gbs_landscape_buffers <- merge(land_100, land_250, by="site")
gbs_landscape_buffers <- merge(gbs_landscape_buffers, land_500, by="site")

write.csv(gbs_landscape_buffers, file="Data/Landscape data/GBS_landscape_buffers_new.csv", row.names=FALSE)

## Step 2: Calculate Euclidean distance from each site to the nearest patch of woodland and grassland

lcm <- st_read("Data/Landscape data/LCM_GBS_5km_final.gpkg") # 5km buffer

st_crs(sites2)

gridref <- sites_final$grid_reference
x_final <- NULL
# change projection to a UK based projection
lcm <- sf::st_transform(lcm, "epsg:32630") # max 2 hours
sites2 <- sf::st_transform(sites2, "epsg:32630") 
st_write(lcm, "Data/Landscape data/LCM_GBS_5km_final_transformed.gpkg")

## Woodland (broadleaf and conifer)
woodland <- lcm[lcm$bhab=="Broadleaf woodland" | lcm$bhab=="Coniferous woodland",]

nearest <- st_nearest_feature(sites2,woodland)
dist = st_distance(sites2, woodland[nearest,], by_element=TRUE)

woodland_distance = cbind(sites2, st_drop_geometry(woodland)[nearest,])
woodland_distance$woodland_dist_m = dist
woodland_distance <- data.frame(woodland_distance)
woodland_distance <- woodland_distance[,c("grid_reference", "woodland_dist_m")]
woodland_distance$woodland_dist_m <- as.numeric(woodland_distance$woodland_dist_m)

## Grassland (improved, calcareous, neutral and heather)
grassland <- lcm[lcm$bhab=="Improved grassland" | lcm$bhab=="Calcareous grassland" | lcm$bhab=="Neutral grassland" | 
                   lcm$bhab=="Heather grassland" | lcm$bhab=="Acid grassland",]

nearest <- st_nearest_feature(sites2,grassland)
dist = st_distance(sites2, grassland[nearest,], by_element=TRUE)

grassland_distance = cbind(sites2, st_drop_geometry(grassland)[nearest,])
grassland_distance$grassland_dist_m = dist
grassland_distance <- data.frame(grassland_distance)
grassland_distance <- grassland_distance[,c("grid_reference", "grassland_dist_m")]
grassland_distance$grassland_dist_m <- as.numeric(grassland_distance$grassland_dist_m)

wood_grass_distance <- merge(woodland_distance, grassland_distance, by="grid_reference")
write.csv(wood_grass_distance, file="Data/Landscape data/Woodland_grassland_patch_distance_new.csv", row.names=FALSE)
