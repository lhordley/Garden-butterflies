##########################
#### user: Lisbeth Hordley
#### date: April 2023
#### info: Surrounding landscape analysis 

options(scipen = 100)

# libraries 
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(DHARMa)
library(glmmTMB)
library(data.table)
library(tidyr)
library(ggeffects)
library(MuMIn)

gbs_analysis <- read.csv("Data/GBS_analysis_final.csv", header=TRUE)
length(unique(gbs_analysis$grid_reference)) # 823 // 785 // 397

# map of the 823 sites for supplementary info
lat_lon <- unique(gbs_analysis[,c("grid_reference","lat_centre","lon_centre"),]) # lat/lon from all GBS sites

worldmap = map_data('world')
all_sites <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = lat_lon, 
             aes(x = as.numeric(lon_centre), 
                 y = as.numeric(lat_centre)), shape=20, size=2) + 
  theme_void() +
  theme(title = element_text(size = 12))
all_sites
ggsave(all_sites, file="Graphs/FigureS1.png")

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden size
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_size)
length(unique(gbs_analysis$grid_reference)) # 805 // 768 // 387

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,] # 142 gardens in size category 4
length(unique(gbs_analysis$grid_reference)) # 665 // 638 // 310

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- subset(gbs_analysis, garden_long_grass_area < 450 | is.na(garden_long_grass_area))
length(unique(gbs_analysis$grid_reference)) # 663 // 636 // 310


##############
### Models ###
##############

###############
## Abundance ##
###############

#### 100m landscape ####

# Check correlation between all variables
round(cor(gbs_analysis[,c("Grassland_100m", "Urban_100m", "Arable_100m", "Woodland_100m", "grassland_dist_m",
                          "woodland_dist_m")]),3)
# urban and grassland -0.8
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_100m) # 0.23
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_100m) # -0.28
# remove grassland from analysis 
abund_pres_100_land <- lmer(log(total_abund) ~ garden_size + n_days + scale(Urban_100m) + scale(Arable_100m) + scale(Woodland_100m) + 
                               scale(grassland_dist_m) + scale(woodland_dist_m) + 
                               scale(I(Urban_100m^2)) + scale(I(Arable_100m^2)) + scale(I(Woodland_100m^2)) + 
                               scale(I(grassland_dist_m^2)) + scale(I(woodland_dist_m^2)) + (1|M_YEAR),
                               na.action = "na.fail", data=gbs_analysis) 
summary(abund_pres_100_land)
# dredge to determine whether each land cover variable is linear or quadratic
abund_pres_100_land_d <- dredge(abund_pres_100_land, subset="garden_size" & "n_days" & "scale(Arable_100m)" & "scale(Urban_100m)" & "scale(Woodland_100m)" & 
                                  "scale(grassland_dist_m)" & "scale(woodland_dist_m)")
abund_pres_100_land_topmods <- subset(abund_pres_100_land_d, delta <2) # take top model set
# model #1 has the lowest number of predictors (df=11) (second model in top rec)
abund_pres_100_land_final <- get.models(abund_pres_100_land_topmods, subset=2)[[1]]
summary(abund_pres_100_land_final) # Urban, grassland distance and woodland quadratic significant (not the linear term though)
# woodland quadratic
# same for climate zone 
# urban linear only for top rec

# save model results
abund_pres_100_land_d <- as.data.frame(coef(summary(abund_pres_100_land_final))) #selecting full model coefficient averages
abund_pres_100_land_d$parameters <- row.names(abund_pres_100_land_d)
row.names(abund_pres_100_land_d) <- 1:nrow(abund_pres_100_land_d)
write.csv(abund_pres_100_land_d, file="Results/Abundance/Abund_land_cover_100m_model.csv", row.names=FALSE)

car::vif(abund_pres_100_land_final) # all <3
r.squaredGLMM(abund_pres_100_land_final) # 14.8%
AIC(abund_pres_100_land_final) # 2799.27

testDispersion(abund_pres_100_land_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_100_land_final, plot = F)
plot(simulationOutput) # no assumptions violated

### Plot graphs ###
abund_pres_urban <- ggpredict(abund_pres_100_land_final, terms="Urban_100m [all]")
abund_pres_urban$landscape <- "Proportion \nof urban \nwithin buffer"
abund_pres_urban$scale <- "100m"
abund_pres_urban$significant <- "yes"

abund_pres_urban_100 <- ggplot() +
  geom_line(data=abund_pres_urban, aes(x = x, y = predicted)) +
  geom_ribbon(data=abund_pres_urban, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  geom_point(data=gbs_analysis, aes(x=Urban_100m, y=total_abund), colour="lightgrey", alpha=0.4)+
  labs(x="Proportion of urban in 100m", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()
abund_pres_urban_100

abund_pres_wood <- ggpredict(abund_pres_100_land_final, terms="Woodland_100m")
abund_pres_wood$landscape <- "Proportion \nof woodland \nwithin buffer"
abund_pres_wood$scale <- "100m"
abund_pres_wood$significant <- "yes"

abund_pres_wood_100 <- ggplot() +
  geom_line(data=abund_pres_wood, aes(x = x, y = predicted)) +
  geom_ribbon(data=abund_pres_wood, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  geom_point(data=gbs_analysis, aes(x=Woodland_100m, y=total_abund), colour="lightgrey", alpha=0.4)+
  labs(x="Proportion of woodland in 100m", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()
abund_pres_wood_100

abund_pres_grass_dist <- ggpredict(abund_pres_100_land_final, terms="grassland_dist_m [all]")
abund_pres_grass_dist$landscape <- "Distance \nto nearest \ngrassland (m)"
abund_pres_grass_dist$scale <- "100m"
abund_pres_grass_dist$significant <- "yes"

abund_pres_grass_dist_100 <- ggplot() +
  geom_line(data=abund_pres_grass_dist, aes(x = x, y = predicted)) +
  geom_ribbon(data=abund_pres_grass_dist, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  geom_point(data=gbs_analysis, aes(x=grassland_dist_m, y=total_abund), colour="lightgrey", alpha=0.4)+
  labs(x="Distance to nearest grassland (m)", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()
abund_pres_grass_dist_100

abund_pres_arable <- ggpredict(abund_pres_100_land_final, terms="Arable_100m [all]")
abund_pres_arable$landscape <- "Proportion \nof arable \nwithin buffer"
abund_pres_arable$scale <- "100m"
abund_pres_arable$significant <- "no"

abund_pres_wood_dist <- ggpredict(abund_pres_100_land_final, terms="woodland_dist_m [all]")
abund_pres_wood_dist$landscape <- "Distance \nto nearest \nwoodland (m)"
abund_pres_wood_dist$scale <- "100m"
abund_pres_wood_dist$significant <- "no"

abund_pres_land_100 <- rbind(abund_pres_urban, abund_pres_arable, abund_pres_wood, abund_pres_grass_dist, abund_pres_wood_dist)



#### 250m landscape ####

# Check correlation between all variables
round(cor(gbs_analysis[,c("Grassland_250m", "Urban_250m", "Arable_250m", "Woodland_250m", "grassland_dist_m",
                          "woodland_dist_m")]),3)
# urban and grassland -0.79
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_250m) # 0.23
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_250m) # -0.31
# remove grassland from analysis 

abund_pres_250_land <- lmer(log(total_abund) ~ garden_size + n_days + scale(Urban_250m) + scale(Arable_250m) + scale(Woodland_250m) + 
                               scale(grassland_dist_m) + scale(woodland_dist_m) + 
                               scale(I(Urban_250m^2)) + scale(I(Arable_250m^2)) + scale(I(Woodland_250m^2)) + 
                               scale(I(grassland_dist_m^2)) + scale(I(woodland_dist_m^2)) + (1|M_YEAR), 
                             na.action = "na.fail", data=gbs_analysis) 

# dredge to determine whether each land cover variable is linear or quadratic
abund_pres_250_land2_d <- dredge(abund_pres_250_land, subset="garden_size" & "n_days" & "scale(Arable_250m)" & "scale(Urban_250m)" & 
                                   "scale(Woodland_250m)" & "scale(grassland_dist_m)" & "scale(woodland_dist_m)")
abund_pres_250_land2_topmods <- subset(abund_pres_250_land2_d, delta <2) # only one top model (AIC<2)
abund_pres_250_land_final <- get.models(abund_pres_250_land2_d, subset=1)[[1]]
summary(abund_pres_250_land_final) # Urban linear significant 
# linear terms only
# arable and urban significant in climate zone subset
# urban only in top rec dataset

# save model results
abund_pres_250_land_d <- as.data.frame(coef(summary(abund_pres_250_land_final))) #selecting full model coefficient averages
abund_pres_250_land_d$parameters <- row.names(abund_pres_250_land_d)
row.names(abund_pres_250_land_d) <- 1:nrow(abund_pres_250_land_d)
write.csv(abund_pres_250_land_d, file="Results/Abundance/Abund_land_cover_250m_model.csv", row.names=FALSE)

car::vif(abund_pres_250_land_final) # all <3
r.squaredGLMM(abund_pres_250_land_final) # 15.9%
AIC(abund_pres_250_land_final) # 2767.44

testDispersion(abund_pres_250_land_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_250_land_final, plot = F)
plot(simulationOutput) # quantile deviations

### Plot graphs ###
abund_pres_urban <- ggpredict(abund_pres_250_land_final, terms="Urban_250m")
abund_pres_urban$landscape <- "Proportion \nof urban \nwithin buffer"
abund_pres_urban$scale <- "250m"
abund_pres_urban$significant <- "yes"

abund_urban_250_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Urban_250m, y=total_abund), colour="lightgrey", alpha=0.4)+
  geom_line(data=abund_pres_urban, aes(x = x, y = predicted)) +
  geom_ribbon(data=abund_pres_urban, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \nurban in 250m", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_urban_250_p

abund_pres_arable <- ggpredict(abund_pres_250_land_final, terms="Arable_250m")
abund_pres_arable$landscape <- "Proportion \nof arable \nwithin buffer"
abund_pres_arable$scale <- "250m"
abund_pres_arable$significant <- "no"

abund_pres_wood <- ggpredict(abund_pres_250_land_final, terms="Woodland_250m [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]")
abund_pres_wood$landscape <- "Proportion \nof woodland \nwithin buffer"
abund_pres_wood$scale <- "250m"
abund_pres_wood$significant <- "no"

abund_pres_grass_dist <- ggpredict(abund_pres_250_land_final, terms="grassland_dist_m")
abund_pres_grass_dist$landscape <- "Distance \nto nearest \ngrassland (m)"
abund_pres_grass_dist$scale <- "250m"
abund_pres_grass_dist$significant <- "no"

abund_pres_wood_dist <- ggpredict(abund_pres_250_land_final, terms="woodland_dist_m")
abund_pres_wood_dist$landscape <- "Distance \nto nearest \nwoodland (m)"
abund_pres_wood_dist$scale <- "250m"
abund_pres_wood_dist$significant <- "no"

abund_pres_land_250 <- rbind(abund_pres_urban, abund_pres_arable, abund_pres_wood, abund_pres_grass_dist, abund_pres_wood_dist)


#### 500m landscape ####

# Check correlation between all variables
round(cor(gbs_analysis[,c("Grassland_500m", "Urban_500m", "Arable_500m", "Woodland_500m", "grassland_dist_m",
                          "woodland_dist_m")]),3)
# urban and grassland -0.76
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_500m) # 0.19
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_500m) # -0.30
# remove grassland from analysis 
abund_pres_500_land <- lmer(log(total_abund) ~ garden_size + n_days + scale(Urban_500m) + scale(Arable_500m) + scale(Woodland_500m) + 
                               scale(grassland_dist_m) + scale(woodland_dist_m) +
                               scale(I(Urban_500m^2)) + scale(I(Arable_500m^2)) + scale(I(Woodland_500m^2)) + 
                               scale(I(grassland_dist_m^2)) + scale(I(woodland_dist_m^2)) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis) 

# dredge to determine whether each land cover variable is linear or quadratic
abund_pres_500_land2_d <- dredge(abund_pres_500_land, subset="garden_size" & "n_days" & "scale(Arable_500m)" & "scale(Urban_500m)" & 
                                   "scale(Woodland_500m)" & "scale(grassland_dist_m)" & "scale(woodland_dist_m)")
abund_pres_500_land2_topmods <- subset(abund_pres_500_land2_d, delta <2)
# only 2 models, second model has lowest number of parameters (df=10) (top model for climate zones)
abund_pres_500_land_final <- get.models(abund_pres_500_land2_topmods, subset=1)[[1]]
summary(abund_pres_500_land_final) # arable and urban significant
# no linear terms
# arable, urban (and quadratic term) and woodland distance significant in climate zone 
# urban linear term only in top rec dataset

# save model results
abund_pres_500_land_d <- as.data.frame(coef(summary(abund_pres_500_land_final))) 
abund_pres_500_land_d$parameters <- row.names(abund_pres_500_land_d)
row.names(abund_pres_500_land_d) <- 1:nrow(abund_pres_500_land_d)
write.csv(abund_pres_500_land_d, file="Results/Abundance/Abund_land_cover_500m_model.csv", row.names=FALSE)

car::vif(abund_pres_500_land_final) 
r.squaredGLMM(abund_pres_500_land_final) # 16.2%
AIC(abund_pres_500_land_final) # 2761.46 (LOWEST AIC)

testDispersion(abund_pres_500_land_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_500_land_final, plot = F)
plot(simulationOutput) # quantile deviations

### Plot graphs ###
abund_urban_500 <- ggpredict(abund_pres_500_land_final, terms="Urban_500m")
abund_urban_500$landscape <- "Proportion \nof urban \nwithin buffer"
abund_urban_500$scale <- "500m"
abund_urban_500$significant <- "yes"

abund_urban_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Urban_500m, y=total_abund), colour="lightgrey", alpha=0.4)+
  geom_line(data=abund_urban_500, aes(x = x, y = predicted)) +
  geom_ribbon(data=abund_urban_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \nurban within 500m", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_urban_500_p

abund_arable_500 <- ggpredict(abund_pres_500_land_final, terms="Arable_500m")
abund_arable_500$landscape <- "Proportion \nof arable \nwithin buffer"
abund_arable_500$scale <- "500m"
abund_arable_500$significant <- "yes"

abund_arable_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Arable_500m, y=total_abund), colour="lightgrey", alpha=0.4)+
  geom_line(data=abund_arable_500, aes(x = x, y = predicted)) +
  geom_ribbon(data=abund_arable_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \narable within 500m", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_arable_500_p

abund_wood_dist_500 <- ggpredict(abund_pres_500_land_final, terms="woodland_dist_m")
abund_wood_dist_500$landscape <- "Distance \nto nearest \nwoodland (m)"
abund_wood_dist_500$scale <- "500m"
abund_wood_dist_500$significant <- "no"

abund_wood_dist_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=woodland_dist_m, y=total_abund), colour="lightgrey", alpha=0.4)+
  geom_line(data=abund_wood_dist_500, aes(x = x, y = predicted), linetype=2) +
  geom_ribbon(data=abund_wood_dist_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Distance to \nnearest woodland (m)", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_wood_dist_500_p

abund_grass_dist_500 <- ggpredict(abund_pres_500_land_final, terms="grassland_dist_m [all]")
abund_grass_dist_500$landscape <- "Distance \nto nearest \ngrassland (m)"
abund_grass_dist_500$scale <- "500m"
abund_grass_dist_500$significant <- "no"

abund_grass_dist_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=grassland_dist_m, y=total_abund), colour="lightgrey", alpha=0.4)+
  geom_line(data=abund_grass_dist_500, aes(x = x, y = predicted), linetype=2) +
  geom_ribbon(data=abund_grass_dist_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Distance to \nnearest grassland (m)", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_grass_dist_500_p

abund_wood_500 <- ggpredict(abund_pres_500_land_final, terms="Woodland_500m [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7]")
abund_wood_500$landscape <- "Proportion \nof woodland \nwithin buffer"
abund_wood_500$scale <- "500m"
abund_wood_500$significant <- "no"

abund_wood_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Woodland_500m, y=total_abund), colour="lightgrey", alpha=0.4)+
  geom_line(data=abund_wood_500, aes(x = x, y = predicted), linetype=2) +
  geom_ribbon(data=abund_wood_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \nwoodland within 500m", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_wood_500_p

abund_pres_land_500 <- rbind(abund_urban_500, abund_arable_500, abund_wood_500, abund_grass_dist_500, abund_wood_dist_500)



##############
## Richness ##
##############

#### 100m landscape ####

# Check correlation between all variables
round(cor(gbs_analysis[,c("Grassland_100m", "Urban_100m", "Arable_100m", "Woodland_100m", "grassland_dist_m",
                          "woodland_dist_m")]),3)
# urban and grassland -0.8
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_100m) # 0.23
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_100m) # -0.29

# remove grassland from analysis 
rich_pres_100_land <- lmer(site_rel_SR ~ garden_size + n_days + scale(Urban_100m) + scale(Arable_100m) + scale(Woodland_100m) + 
                              scale(grassland_dist_m) + scale(woodland_dist_m) + 
                              scale(I(Urban_100m^2)) + scale(I(Arable_100m^2)) + scale(I(Woodland_100m^2)) + 
                              scale(I(grassland_dist_m^2)) + scale(I(woodland_dist_m^2)) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis) 

# dredge to determine whether each land cover variable is linear or quadratic
rich_pres_100_land2_d <- dredge(rich_pres_100_land, subset="garden_size" & "n_days" & "scale(Arable_100m)" & "scale(Urban_100m)" & 
                                  "scale(Woodland_100m)" & "scale(grassland_dist_m)" & "scale(woodland_dist_m)")
rich_pres_100_land2_topmods <- subset(rich_pres_100_land2_d, delta <2)
# top set only has one model
rich_pres_100_land_final <- get.models(rich_pres_100_land2_d, subset=1)[[1]]
summary(rich_pres_100_land_final) # urban, woodland and woodland distance
# all linear terms
# urban, woodland, woodland distance and grassland distance significant in climate zone data
# urban, woodland and woodland distance in top rec dataset

# save model results
rich_pres_100_land_d <- as.data.frame(coef(summary(rich_pres_100_land_final))) #selecting full model coefficient averages
rich_pres_100_land_d$parameters <- row.names(rich_pres_100_land_d)
row.names(rich_pres_100_land_d) <- 1:nrow(rich_pres_100_land_d)
write.csv(rich_pres_100_land_d, file="Results/Richness/Richness_land_cover_100m_model.csv", row.names=FALSE)

car::vif(rich_pres_100_land_final) # all <3 
r.squaredGLMM(rich_pres_100_land_final) # 27.6%
AIC(rich_pres_100_land_final) # -4764.91

testDispersion(rich_pres_100_land_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_100_land_final, plot = F)
plot(simulationOutput) # no assumptions violated

### Plot graphs ###

rich_urban_100 <- ggpredict(rich_pres_100_land_final, terms="Urban_100m")
rich_urban_100$landscape <- "Proportion \nof urban \nwithin buffer"
rich_urban_100$scale <- "100m"
rich_urban_100$significant <- "yes"

rich_urban_100_p <- ggplot() +
  geom_line(data=rich_urban_100, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_urban_100, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  geom_point(data=gbs_analysis, aes(x=Urban_100m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  labs(x="Proportion of urban within 100m", y="Relative richness")+
  theme_classic()
rich_urban_100_p

rich_wood_100 <- ggpredict(rich_pres_100_land_final, terms="Woodland_100m [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]")
rich_wood_100$landscape <- "Proportion \nof woodland \nwithin buffer"
rich_wood_100$scale <- "100m"
rich_wood_100$significant <- "yes"

rich_arable_100 <- ggpredict(rich_pres_100_land_final, terms="Arable_100m")
rich_arable_100$landscape <- "Proportion \nof arable \nwithin buffer"
rich_arable_100$scale <- "100m"
rich_arable_100$significant <- "no"

rich_grass_dist_100 <- ggpredict(rich_pres_100_land_final, terms="grassland_dist_m [all]")
rich_grass_dist_100$landscape <- "Distance \nto nearest \ngrassland (m)"
rich_grass_dist_100$scale <- "100m"
rich_grass_dist_100$significant <- "no"

rich_wood_dist_100 <- ggpredict(rich_pres_100_land_final, terms="woodland_dist_m [all]")
rich_wood_dist_100$landscape <- "Distance \nto nearest \nwoodland (m)"
rich_wood_dist_100$scale <- "100m"
rich_wood_dist_100$significant <- "yes"

rich_pres_land_100 <- rbind(rich_urban_100, rich_arable_100, rich_wood_100, rich_grass_dist_100, rich_wood_dist_100)

#### 250m landscape ####

# Check correlation between all variables
round(cor(gbs_analysis[,c("Grassland_250m", "Urban_250m", "Arable_250m", "Woodland_250m", "grassland_dist_m",
                          "woodland_dist_m")]),3)
# urban and grassland -0.79
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_250m) # 0.21
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_250m) # -0.32
# remove grassland from analysis 
rich_pres_250_land <- lmer(site_rel_SR ~ garden_size + n_days + scale(Urban_250m) + scale(Arable_250m) + scale(Woodland_250m) + 
                              scale(grassland_dist_m) + scale(woodland_dist_m) + 
                              scale(I(Urban_250m^2)) + scale(I(Arable_250m^2)) + scale(I(Woodland_250m^2)) + 
                              scale(I(grassland_dist_m^2)) + scale(I(woodland_dist_m^2)) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis)

# dredge to determine whether each land cover variable is linear or quadratic
rich_pres_250_land2_d <- dredge(rich_pres_250_land, subset="garden_size" & "n_days" & "scale(Arable_250m)" & "scale(Urban_250m)" & "scale(Woodland_250m)" & 
                                  "scale(grassland_dist_m)" & "scale(woodland_dist_m)")
rich_pres_250_land2_topmods <- subset(rich_pres_250_land2_d, delta <2)
# top model has lowest number of predictors
rich_pres_250_land_final <- get.models(rich_pres_250_land2_d, subset=1)[[1]]
summary(rich_pres_250_land_final) # arable, urban, woodland and woodland distance significant
# all linear terms
# same results for climate zone data
# urban, woodland and woodland distance for top rec data

# save model results
rich_pres_250_land_d <- as.data.frame(coef(summary(rich_pres_250_land_final))) 
rich_pres_250_land_d$parameters <- row.names(rich_pres_250_land_d)
row.names(rich_pres_250_land_d) <- 1:nrow(rich_pres_250_land_d)
write.csv(rich_pres_250_land_d, file="Results/Richness/Richness_land_cover_250m_model.csv", row.names=FALSE)

car::vif(rich_pres_250_land_final) # all <3
r.squaredGLMM(rich_pres_250_land_final) # 29%
AIC(rich_pres_250_land_final) # -4803.02 (LOWEST AIC)

testDispersion(rich_pres_250_land_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_250_land_final, plot = F)
plot(simulationOutput) # quantile deviations

### Plot graphs ###
rich_arable_250 <- ggpredict(rich_pres_250_land_final, terms="Arable_250m [all]")
rich_arable_250$landscape <- "Proportion \nof arable \nwithin buffer"
rich_arable_250$scale <- "250m"
rich_arable_250$significant <- "yes"

rich_arable_250_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Arable_250m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_arable_250, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_arable_250, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \narable within 250m", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_arable_250_p

rich_urban_250 <- ggpredict(rich_pres_250_land_final, terms="Urban_250m")
rich_urban_250$landscape <- "Proportion \nof urban \nwithin buffer"
rich_urban_250$scale <- "250m"
rich_urban_250$significant <- "yes"

rich_urban_250_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Urban_250m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_urban_250, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_urban_250, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \nurban within 250m", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_urban_250_p

rich_wood_250 <- ggpredict(rich_pres_250_land_final, terms="Woodland_250m [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]")
rich_wood_250$landscape <- "Proportion \nof woodland \nwithin buffer"
rich_wood_250$scale <- "250m"
rich_wood_250$significant <- "yes"

rich_wood_250_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Woodland_250m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_wood_250, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_wood_250, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \nwoodland within 250m", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_wood_250_p

rich_wood_dist_250 <- ggpredict(rich_pres_250_land_final, terms="woodland_dist_m")
rich_wood_dist_250$landscape <- "Distance \nto nearest \nwoodland (m)"
rich_wood_dist_250$scale <- "250m"
rich_wood_dist_250$significant <- "yes"

rich_wood_dist_250_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=woodland_dist_m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_wood_dist_250, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_wood_dist_250, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Distance to \nnearest woodland (m)", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_wood_dist_250_p

rich_grass_dist_250 <- ggpredict(rich_pres_250_land_final, terms="grassland_dist_m")
rich_grass_dist_250$landscape <- "Distance \nto nearest \ngrassland (m)"
rich_grass_dist_250$scale <- "250m"
rich_grass_dist_250$significant <- "no"

rich_grass_dist_250_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=grassland_dist_m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_grass_dist_250, aes(x = x, y = predicted), linetype=2) +
  geom_ribbon(data=rich_grass_dist_250, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Distance to \nnearest grassland (m)", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_grass_dist_250_p

rich_pres_land_250 <- rbind(rich_urban_250, rich_arable_250, rich_wood_250, rich_grass_dist_250, rich_wood_dist_250)

#### 500m landscape ####

# Check correlation between all variables
round(cor(gbs_analysis[,c("Grassland_500m", "Urban_500m", "Arable_500m", "Woodland_500m", "grassland_dist_m",
                          "woodland_dist_m")]),3)
# urban and grassland -0.76
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_500m) # 0.20
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_500m) # -0.33
# remove grassland from analysis 
rich_pres_500_land <- lmer(site_rel_SR ~ garden_size + n_days + scale(Urban_500m) + scale(Arable_500m) + scale(Woodland_500m) + scale(grassland_dist_m) + 
                              scale(woodland_dist_m)+ scale(I(Urban_500m^2)) + scale(I(Arable_500m^2)) + 
                              scale(I(Woodland_500m^2)) + scale(I(grassland_dist_m^2)) + scale(I(woodland_dist_m^2)) + 
                              (1|M_YEAR), na.action = "na.fail", data=gbs_analysis) 

# dredge to determine whether each land cover variable is linear or quadratic
rich_pres_500_land2_d <- dredge(rich_pres_500_land, subset="garden_size" & "n_days" & "scale(Arable_500m)" & "scale(Urban_500m)" & "scale(Woodland_500m)" & 
                                  "scale(grassland_dist_m)" & "scale(woodland_dist_m)")
rich_pres_500_land2_topmods <- subset(rich_pres_500_land2_d, delta <2)
# model #1 in top set has lowest number of parameters
rich_pres_500_land_final <- get.models(rich_pres_500_land2_d, subset=1)[[1]]
summary(rich_pres_500_land_final) # arable, urban, woodland, woodland distance and grassland distance significant
# all linear terms
# arable, urban, woodland and woodland distance significant in climate zone data
# urban, woodland and woodland distance significant in top rec data

# save model results
rich_pres_500_land_d <- as.data.frame(coef(summary(rich_pres_500_land_final))) 
rich_pres_500_land_d$parameters <- row.names(rich_pres_500_land_d)
row.names(rich_pres_500_land_d) <- 1:nrow(rich_pres_500_land_d)
write.csv(rich_pres_500_land_d, file="Results/Richness/Richness_land_cover_500m_model.csv", row.names=FALSE)

car::vif(rich_pres_500_land_final) 
r.squaredGLMM(rich_pres_500_land_final) # 28.8%
AIC(rich_pres_500_land_final) # -4798.47

testDispersion(rich_pres_500_land_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_500_land_final, plot = F)
plot(simulationOutput) # quantile deviations

### Plot graphs ###

rich_arable_500 <- ggpredict(rich_pres_500_land_final, terms="Arable_500m")
rich_arable_500$landscape <- "Proportion \nof arable \nwithin buffer"
rich_arable_500$scale <- "500m"
rich_arable_500$significant <- "yes"

rich_arable_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Arable_500m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_arable_500, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_arable_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \narable within 500m", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_arable_500_p

rich_wood_500 <- ggpredict(rich_pres_500_land_final, terms="Woodland_500m [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7]")
rich_wood_500$landscape <- "Proportion \nof woodland \nwithin buffer"
rich_wood_500$scale <- "500m"
rich_wood_500$significant <- "yes"

rich_wood_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Woodland_500m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_wood_500, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_wood_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \nwoodland within 500m", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_wood_500_p

rich_urban_500 <- ggpredict(rich_pres_500_land_final, terms="Urban_500m")
rich_urban_500$landscape <- "Proportion \nof urban \nwithin buffer"
rich_urban_500$scale <- "500m"
rich_urban_500$significant <- "yes"

rich_urban_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=Urban_500m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_urban_500, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_urban_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Proportion of \nurban within 500m", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_urban_500_p

rich_grass_dist_500 <- ggpredict(rich_pres_500_land_final, terms="grassland_dist_m")
rich_grass_dist_500$landscape <- "Distance \nto nearest \ngrassland (m)"
rich_grass_dist_500$scale <- "500m"
rich_grass_dist_500$significant <- "yes"

rich_wood_dist_500 <- ggpredict(rich_pres_500_land_final, terms="woodland_dist_m")
rich_wood_dist_500$landscape <- "Distance \nto nearest \nwoodland (m)"
rich_wood_dist_500$scale <- "500m"
rich_wood_dist_500$significant <- "yes"

rich_wood_dist_500_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x=woodland_dist_m, y=site_rel_SR), colour="lightgrey", alpha=0.4)+
  geom_line(data=rich_wood_dist_500, aes(x = x, y = predicted)) +
  geom_ribbon(data=rich_wood_dist_500, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  labs(x="Distance to \nnearest woodland (m)", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_wood_dist_500_p

rich_pres_land_500 <- rbind(rich_urban_500, rich_arable_500, rich_wood_500, rich_grass_dist_500, rich_wood_dist_500)

# Lowest AIC models: 500m for abundance, 250m for richness
# put these plots together for main text
library(ggpubr)
landscape_all <- ggarrange(abund_arable_500_p, abund_urban_500_p, abund_wood_500_p, abund_wood_dist_500_p, abund_grass_dist_500_p, 
                           rich_arable_250_p, rich_urban_250_p, rich_wood_250_p, rich_wood_dist_250_p,
                           rich_grass_dist_250_p, 
                           labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)"),
                           ncol = 5, nrow = 2)
landscape_all
ggsave(landscape_all, file="Graphs/Landscape/Landsacpe_plots_Figure2.png", height=7, width=15)

## Plot all landscape results for supplementary material 

# Abundance

abund_pres_land <- rbind(abund_pres_land_100,abund_pres_land_250,abund_pres_land_500)
abund_pres_land$significant  <- factor(abund_pres_land$significant,levels = c("yes", "no"))
abund_pres_land$landscape  <- factor(abund_pres_land$landscape,levels = c("Proportion \nof arable \nwithin buffer", 
                                                                          "Proportion \nof urban \nwithin buffer", "Proportion \nof woodland \nwithin buffer", 
                                                                          "Distance \nto nearest \nwoodland (m)", "Distance \nto nearest \ngrassland (m)"))


equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    if(seq[2]-seq[1] < 10^(-r)) seq else round(seq, r)
  }
}

plot <- ggplot()+
  geom_line(data=abund_pres_land, aes(x=x, y=predicted, linetype=significant), lwd=1)+
  geom_ribbon(data=abund_pres_land, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  #geom_point(data=gbs_land_data2, aes(x=value, y=total_abund))+
  facet_grid(scale ~ landscape, scales = "free", switch = 'x')+
  #scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  labs(y="Total abundance")+
  scale_x_continuous(breaks=equal_breaks(n=3, s=0.05, r=0))+
  theme_bw()+
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=24),
        panel.spacing = unit(1.5, "lines"),
        legend.position = "none")
plot
ggsave(plot, file="Graphs/Landscape/Abundance_landscape.png", height=12, width=14)

# Richness
rich_pres_land <- rbind(rich_pres_land_100,rich_pres_land_250,rich_pres_land_500)
rich_pres_land$significant  <- factor(rich_pres_land$significant,levels = c("yes", "no"))
rich_pres_land$landscape  <- factor(rich_pres_land$landscape,levels = c("Proportion \nof arable \nwithin buffer", 
                                                                        "Proportion \nof urban \nwithin buffer", "Proportion \nof woodland \nwithin buffer", 
                                                                        "Distance \nto nearest \nwoodland (m)", "Distance \nto nearest \ngrassland (m)"))


equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    if(seq[2]-seq[1] < 10^(-r)) seq else round(seq, r)
  }
}

plot <- ggplot()+
  geom_line(data=rich_pres_land, aes(x=x, y=predicted, linetype=significant), lwd=1)+
  geom_ribbon(data=rich_pres_land, aes(x = x, y = predicted, ymin=conf.low, ymax=conf.high), alpha=.3, linetype = 0)+
  #geom_point(data=gbs_land_data2, aes(x=value, y=total_abund))+
  facet_grid(scale ~ landscape, scales = "free", switch = 'x')+
  #scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  labs(y="Relative richness")+
  scale_x_continuous(n.breaks = 3)+
  theme_bw()+
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        text=element_text(size=24),
        panel.spacing = unit(1.5, "lines"),
        legend.position = "none")
plot
ggsave(plot, file="Graphs/Landscape/Richness_landscape.png", height=12, width=14)
