##########################
#### user: Lisbeth Hordley
#### date: November 2022
#### info: Calculate landscape buffers around GBS sites

# libraries 
require(rgdal)
library(mapview)
library(sf)
library(sp)
library(rnrfa)
library(dplyr)
options(scipen = 100)

lcm <- st_read("Data/LCM_GBS_1000m_clip.gpkg")


# get GBS sites and change to spatial something?
sites <- read.csv("Data/GBS_site_info.csv", header=TRUE)

#coordinates(sites) <- c("lon", "lat")
sites2 = st_as_sf(sites,coords=c("lon","lat"), dim = "XY", crs=4326)
coords <- st_geometry(sites2)
st_crs(coords)

buff_coords <- st_buffer(coords, dist = 1000)

mapview(coords) + mapview(buff_coords)
mapview(lcm) + mapview(buff_coords)



gridref <- sites$grid_reference
x_final <- NULL
# change projection to a UK based projection
lcm <- sf::st_transform(lcm, "epsg:32630")
sites2 <- sf::st_transform(sites2, "epsg:32630")

for(i in gridref[1:45]){print(i)
  sites_temp <- sites2[sites2$grid_reference==i,]
    x <- st_intersection(lcm, sites_temp |> st_geometry() |> st_buffer(dist = 100))
    x["area"] <- st_area(x)
    x_area <- aggregate(area ~ bhab, x, sum)
    x_area$site <- i
    x_final <- rbind(x_area, x_final)
  }

x_final$area <- as.numeric(x_final$area)
x_final2 <- reshape(x_final, idvar = c("site"), timevar = "bhab", direction = "wide")
x_final2[is.na(x_final2)] <- 0
x_final2$total_count <- rowSums(x_final2[2:13]) # create total_count column
x_final3 <- x_final2[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)]
x_final3$total_count <- rowSums(x_final3[2:16]) # create total_count column

mapview(lcm) + mapview(buff_coords)
