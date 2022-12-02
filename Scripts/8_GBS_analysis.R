##########################
#### user: Lisbeth Hordley
#### date: November 2022
#### info: Initial GBS abundance/richness and landscape analysis
options(scipen = 100)

# libraries 
library(ggplot2)
library(dplyr)

# load data
garden_features <- read.csv("Data/Raw GBS data/GBS_garden_features_cleaned.csv", header=TRUE)
gbs_abund <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)
gbs_richness <- read.csv("Data/Richness data/Relative_SR_GBS.csv", header=TRUE)
gbs_landscape <- read.csv("Data/Landscape data/Land_cover_GBS_100m.csv", header=TRUE)

# Need to take sum of species abundance indices per site
gbs_abund <- gbs_abund %>% group_by(grid_reference, M_YEAR) %>% summarise(total_abund=sum(SINDEX))
# remove outlier 
gbs_abund <- gbs_abund[!gbs_abund$total_abund>2000,]

# merge garden and landscape data with abundance and richness
gbs_abund <- merge(gbs_abund, garden_features, by.x=c("grid_reference", "M_YEAR"), by.y=c("grid_reference", "year"), all.x=TRUE)
# gbs_abund <- merge(gbs_abund, gbs_landscape, by.x="grid_reference", by.y="site", all.x=TRUE)
gbs_abund$M_YEAR <- as.factor(gbs_abund$M_YEAR)
gbs_abund$garden_long_grass_area <- as.numeric(gbs_abund$garden_long_grass_area)

# some plots with garden features
ggplot(gbs_abund[!is.na(gbs_abund$garden_size),], aes(x=garden_size, y=total_abund))+
  geom_point()+
  theme_classic()
# larger gardens = higher site abundance?

ggplot(gbs_abund[!is.na(gbs_abund$garden_direction),], aes(x=garden_direction, y=total_abund))+
  geom_point()+
  theme_classic()
# gardens facing south have highest (or at least greatest variation) in abundance 

ggplot(gbs_abund[!is.na(gbs_abund$garden_long_grass),], aes(x=garden_long_grass, y=total_abund))+
  geom_point()+
  theme_classic()
# gardens with long grass have higher abundance

x <- gbs_abund[gbs_abund$garden_long_grass_area>0,]
x <- x[!is.na(x$garden_long_grass_area),]
ggplot(x, aes(x=garden_long_grass_area, y=total_abund))+
  geom_point()+
  theme_classic()
# possible positive relationship but with a few outliers

ggplot(gbs_abund[!is.na(gbs_abund$garden_hard_surface),], aes(x=garden_hard_surface, y=total_abund))+
  geom_point()+
  theme_classic()
# all gardens put down 1 - less than a quarter. Odd or just a coincidence? 

ggplot(gbs_abund[!is.na(gbs_abund$garden_ivy),], aes(x=garden_ivy, y=total_abund))+
  geom_point()+
  theme_classic()
# higher abundance with flowering ivy (but this will be more useful as a species-specific analysis)

### Richness and garden features
gbs_richness <- merge(gbs_richness, garden_features, by=c("grid_reference", "year"), all.x=TRUE)

# some plots with garden features
ggplot(gbs_richness[!is.na(gbs_richness$garden_size),], aes(x=garden_size, y=site_rel_SR))+
  geom_point()+
  theme_classic()
# possible pattern: larger gardens = higher richness

ggplot(gbs_richness[!is.na(gbs_richness$garden_direction),], aes(x=garden_direction, y=site_rel_SR))+
  geom_point()+
  theme_classic()
# no real pattern here 

ggplot(gbs_richness[!is.na(gbs_richness$garden_long_grass),], aes(x=garden_long_grass, y=site_rel_SR))+
  geom_point()+
  theme_classic()
# gardens with long grass possibly higher SR

ggplot(gbs_richness[!is.na(gbs_richness$garden_long_grass_area),], aes(x=garden_long_grass_area, y=site_rel_SR))+
  geom_point()+
  theme_classic()
# positive relationship

ggplot(gbs_richness[!is.na(gbs_richness$garden_hard_surface),], aes(x=garden_hard_surface, y=site_rel_SR))+
  geom_point()+
  theme_classic()
# all gardens put down 1 - less than a quarter. Odd or just a coincidence? 

ggplot(gbs_richness[!is.na(gbs_richness$garden_ivy),], aes(x=garden_ivy, y=site_rel_SR))+
  geom_point()+
  theme_classic()
# higher abundance with flowering ivy (but this will be more useful as a species-specific analysis)






##### Landscape and abundance
gbs_landscape[,3:20] <- sapply(gbs_landscape[,3:20],as.numeric)

gbs_abund_land <- merge(gbs_abund, gbs_landscape, by.x="grid_reference", by.y="site")
plot(gbs_abund_land$area.Urban, gbs_abund_land$total_abund) 
plot(gbs_abund_land$area.Suburban, gbs_abund_land$total_abund) 

plot(gbs_abund_land$area.Broadleaf.woodland, gbs_abund_land$total_abund) 
plot(gbs_abund_land$area.Arable.and.horticulture, gbs_abund_land$total_abund) 
plot(gbs_abund_land$area.Improved.grassland, gbs_abund_land$total_abund) 

library(lme4)
library(lmerTest)
mod <- lmer(log(total_abund) ~ area.Urban + (1|M_YEAR), data=gbs_abund_land)
summary(mod)

hist(resid(mod)) # log skewed is normally distributed
qqnorm(resid(mod))
qqline(resid(mod))

library(ggeffects)
ggpredict(mod, terms = "area.Urban", back.transform = TRUE) %>% plot(add.data = TRUE)

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






















london_gardens <- garden_features_final2[garden_features_final2$lat >= 51 & garden_features_final2$lat <= 52 &
                                           garden_features_final2$lon >= -0.9 & garden_features_final2$lon <= 0.5,]

library(ggmap)
lats<-c(51,52)
lons<-c(-0.9,0.5)
bb<-make_bbox(lon=lons,lat=lats,f=0.05)
cda<-get_map(bb,source="osm")
london_map <- ggmap(cda)+geom_point(data=london_gardens, aes(x=lon, y=lat), color="red", size=2, alpha=0.5)
london_map
