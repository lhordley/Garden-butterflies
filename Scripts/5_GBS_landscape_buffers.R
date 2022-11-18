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
library(stringr)
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
scales <- c(100,250,500)
# change projection to a UK based projection
lcm <- sf::st_transform(lcm, "epsg:32630")
sites2 <- sf::st_transform(sites2, "epsg:32630")

for(i in gridref){print(i)
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
write.csv(x_final, file="Data/gbs_land_100_250_500.csv", row.names=FALSE)

x_final <- read.csv("Data/gbs_land_100_250_500.csv", header=TRUE)
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

# divide all columns by total land area
habs <- colnames(land_100)
habs <- habs[-c(1:2)]

for (i in habs){
  land_100[,i] <- land_100[,i]/land_100[,"total_count"]
}

habs2 <- colnames(land_250)
habs2 <- habs2[-c(1:2)]

for (i in habs2){
  land_250[,i] <- land_250[,i]/land_250[,"total_count"]
}

habs3 <- colnames(land_500)
habs3 <- habs3[-c(1:2)]

for (i in habs3){
  land_500[,i] <- land_500[,i]/land_500[,"total_count"]
}

write.csv(land_100, file="Data/Land_cover_GBS_100m.csv", row.names=FALSE)
write.csv(land_250, file="Data/Land_cover_GBS_250m.csv", row.names=FALSE)
write.csv(land_500, file="Data/Land_cover_GBS_500m.csv", row.names=FALSE)




