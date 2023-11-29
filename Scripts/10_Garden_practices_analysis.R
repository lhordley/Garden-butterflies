##########################
#### user: Lisbeth Hordley
#### date: November 2022
#### info: Wildlife garden practice analysis 
options(scipen = 100)

# libraries 
library(ggplot2)
library(broom)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(DHARMa)
library(ggeffects)
library(data.table)


###########################################################
######## Garden practice 1: Presence of long grass ########
###########################################################

####### 1a: all species, all sites ####### 

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden presence of long grass and garden size
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass, garden_size)
length(unique(gbs_analysis$grid_reference)) # 782 sites 

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 649 gardens 

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 647 gardens 

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed


abund_pres_garden <- lmer(log(total_abund) ~ garden_size + garden_long_grass + n_days + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis) 
summary(abund_pres_garden) # all significant positive

car::vif(abund_pres_garden) # all <3
r.squaredGLMM(abund_pres_garden) # 5.8%
AIC(abund_pres_garden) # 2844.591

testDispersion(abund_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_pres_garden_d <- as.data.frame(coef(summary(abund_pres_garden))) #selecting full model coefficient averages
abund_pres_garden_d$parameters <- row.names(abund_pres_garden_d)
row.names(abund_pres_garden_d) <- 1:nrow(abund_pres_garden_d)
abund_pres_garden_d$AIC <- AIC(abund_pres_garden)
write.csv(abund_pres_garden_d, file="Results/Abundance/Abund_presence_long_grass_model.csv")

# Plot graph
abund_long_grass <- ggpredict(abund_pres_garden, terms="garden_long_grass")
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==0, "Absent", "Present") # 0 = no long grass, 1 = long grass
abund_long_grass$x <- ifelse(abund_long_grass$x==0, "Absent", "Present") # 0 = no long grass, 1 = long grass

abund_long_grass_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x = garden_long_grass, y=total_abund), colour="lightgrey", alpha=0.4,
             position=position_jitter(w = 0.2, h = 0)) +
  geom_point(data = abund_long_grass, aes(x = x, y = predicted), 
             size = 1,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(data = abund_long_grass, 
                aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low), 
                position = position_dodge(width = 0.4),
                width = 0.1) +
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  labs(x="Long grass", y="Total abundance")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_long_grass_p
ggsave(abund_long_grass_p, file="Graphs/Abundance/Abund_long_grass.png", height=4, width=4)

## Richness model 

rich_pres_garden <- lmer(site_rel_SR ~ garden_size + garden_long_grass + n_days + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(rich_pres_garden) 

car::vif(rich_pres_garden) # all <3
r.squaredGLMM(rich_pres_garden) # 20%
AIC(rich_pres_garden) # -4539.91

testDispersion(rich_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_pres_garden_d <- as.data.frame(coef(summary(rich_pres_garden))) #selecting full model coefficient averages
rich_pres_garden_d$parameters <- row.names(rich_pres_garden_d)
row.names(rich_pres_garden_d) <- 1:nrow(rich_pres_garden_d)
rich_pres_garden_d$AIC <- AIC(rich_pres_garden)
write.csv(rich_pres_garden_d, file="Results/Richness/Richness_presence_long_grass_model.csv")

# Plot graph
# gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==0, "Absent", "Present") # 0 = no long grass, 1 = long grass
rich_long_grass <- ggpredict(rich_pres_garden, terms="garden_long_grass")

rich_long_grass_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x = garden_long_grass, y=site_rel_SR), colour="lightgrey", alpha=0.4,
             position=position_jitter(w = 0.2, h = 0)) +
  geom_point(data = rich_long_grass, aes(x = x, y = predicted), 
             size = 1,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(data = rich_long_grass, 
                aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low), 
                position = position_dodge(width = 0.4),
                width = 0.1) +
  labs(x="Long grass", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_long_grass_p
ggsave(rich_long_grass_p, file="Graphs/Richness/Rich_long_grass.png", height=4, width=4)


####### 1b: grassland species, all sites #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_grass.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden size and presence of long grass
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass, garden_size)
length(unique(gbs_analysis$grid_reference)) # 761

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$n_days <- as.numeric(gbs_analysis$n_days)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 630 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 628 gardens 

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_pres_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis) 
summary(abund_pres_garden) # both significant positive

car::vif(abund_pres_garden) # all <3
r.squaredGLMM(abund_pres_garden) # 5.8%
AIC(abund_pres_garden) # 2844.591

testDispersion(abund_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_pres_garden_d <- as.data.frame(coef(summary(abund_pres_garden))) #selecting full model coefficient averages
abund_pres_garden_d$parameters <- row.names(abund_pres_garden_d)
row.names(abund_pres_garden_d) <- 1:nrow(abund_pres_garden_d)
abund_pres_garden_d$AIC <- AIC(abund_pres_garden)
write.csv(abund_pres_garden_d, file="Results/Abundance/Abund_presence_long_grass_model_grass.csv")

## Richness model 

rich_pres_garden <- lmer(site_rel_SR ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(rich_pres_garden) 

car::vif(rich_pres_garden) # all <3
r.squaredGLMM(rich_pres_garden) # 9.8%
AIC(rich_pres_garden2) # -4245.485

testDispersion(rich_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_pres_garden_d <- as.data.frame(coef(summary(rich_pres_garden))) #selecting full model coefficient averages
rich_pres_garden_d$parameters <- row.names(rich_pres_garden_d)
row.names(rich_pres_garden_d) <- 1:nrow(rich_pres_garden_d)
rich_pres_garden_d$AIC <- AIC(rich_pres_garden)
write.csv(rich_pres_garden_d, file="Results/Richness/Richness_presence_long_grass_model_grass.csv")

####### 1c: non-grassland species, all sites ####### 

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_not_grass.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden size and presence of long grass
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass, garden_size)
length(unique(gbs_analysis$grid_reference)) # 783

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 649 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 647 gardens 

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_pres_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis) 
summary(abund_pres_garden) # both significant positive

car::vif(abund_pres_garden) # all <3
r.squaredGLMM(abund_pres_garden) # 5.8%
AIC(abund_pres_garden) # 2844.591

testDispersion(abund_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_pres_garden_d <- as.data.frame(coef(summary(abund_pres_garden))) #selecting full model coefficient averages
abund_pres_garden_d$parameters <- row.names(abund_pres_garden_d)
row.names(abund_pres_garden_d) <- 1:nrow(abund_pres_garden_d)
abund_pres_garden_d$AIC <- AIC(abund_pres_garden)
write.csv(abund_pres_garden_d, file="Results/Abundance/Abund_presence_long_grass_model_not_grass.csv")

## Richness model 

rich_pres_garden <- lmer(site_rel_SR ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(rich_pres_garden) 

car::vif(rich_pres_garden) # all <3
r.squaredGLMM(rich_pres_garden) # 9.8%
AIC(rich_pres_garden2) # -4245.485

testDispersion(rich_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_pres_garden_d <- as.data.frame(coef(summary(rich_pres_garden))) #selecting full model coefficient averages
rich_pres_garden_d$parameters <- row.names(rich_pres_garden_d)
row.names(rich_pres_garden_d) <- 1:nrow(rich_pres_garden_d)
rich_pres_garden_d$AIC <- AIC(rich_pres_garden)
write.csv(rich_pres_garden_d, file="Results/Richness/Richness_presence_long_grass_model_not_grass.csv")


####### 1d: all species, climate zone sites ####### 

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_GEnS.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden features (apart from ivy - only need that data later on)
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass, garden_size)
length(unique(gbs_analysis$grid_reference)) # 748

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 623 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 621 gardens 

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_pres_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis) 
summary(abund_pres_garden) # both significant positive

car::vif(abund_pres_garden) # all <3
r.squaredGLMM(abund_pres_garden) # 5.8%
AIC(abund_pres_garden) # 2844.591

testDispersion(abund_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_pres_garden_d <- as.data.frame(coef(summary(abund_pres_garden))) #selecting full model coefficient averages
abund_pres_garden_d$parameters <- row.names(abund_pres_garden_d)
row.names(abund_pres_garden_d) <- 1:nrow(abund_pres_garden_d)
abund_pres_garden_d$AIC <- AIC(abund_pres_garden)
write.csv(abund_pres_garden_d, file="Results/Abundance/Abund_presence_long_grass_model_GEnS.csv")

## Richness model 

rich_pres_garden <- lmer(site_rel_SR ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(rich_pres_garden) 

car::vif(rich_pres_garden) # all <3
r.squaredGLMM(rich_pres_garden) # 9.8%
AIC(rich_pres_garden2) # -4245.485

testDispersion(rich_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_pres_garden_d <- as.data.frame(coef(summary(rich_pres_garden))) #selecting full model coefficient averages
rich_pres_garden_d$parameters <- row.names(rich_pres_garden_d)
row.names(rich_pres_garden_d) <- 1:nrow(rich_pres_garden_d)
rich_pres_garden_d$AIC <- AIC(rich_pres_garden)
write.csv(rich_pres_garden_d, file="Results/Richness/Richness_presence_long_grass_model_GEnS.csv")


####### 1e: all species, top 50% sites ####### 

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_top_rec.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden size and presence of long grass
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass, garden_size)
length(unique(gbs_analysis$grid_reference)) # 382

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 305 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 305 gardens 

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_pres_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis) 
summary(abund_pres_garden) # both significant positive

car::vif(abund_pres_garden) # all <3
r.squaredGLMM(abund_pres_garden) # 5.8%
AIC(abund_pres_garden) # 2844.591

testDispersion(abund_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_pres_garden_d <- as.data.frame(coef(summary(abund_pres_garden))) #selecting full model coefficient averages
abund_pres_garden_d$parameters <- row.names(abund_pres_garden_d)
row.names(abund_pres_garden_d) <- 1:nrow(abund_pres_garden_d)
abund_pres_garden_d$AIC <- AIC(abund_pres_garden)
write.csv(abund_pres_garden_d, file="Results/Abundance/Abund_presence_long_grass_model_top_rec.csv")

## Richness model 

rich_pres_garden <- lmer(site_rel_SR ~ garden_size + n_days + garden_long_grass + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(rich_pres_garden) 

car::vif(rich_pres_garden) # all <3
r.squaredGLMM(rich_pres_garden) # 9.8%
AIC(rich_pres_garden2) # -4245.485

testDispersion(rich_pres_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_pres_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_pres_garden_d <- as.data.frame(coef(summary(rich_pres_garden))) #selecting full model coefficient averages
rich_pres_garden_d$parameters <- row.names(rich_pres_garden_d)
row.names(rich_pres_garden_d) <- 1:nrow(rich_pres_garden_d)
rich_pres_garden_d$AIC <- AIC(rich_pres_garden)
write.csv(rich_pres_garden_d, file="Results/Richness/Richness_presence_long_grass_model_top_rec.csv")



###########################################################
########## Garden practice 2: Area of long grass ########## 
###########################################################

####### 2a: all species, all sites ####### 

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from area of long grass and garden size
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass_area, garden_size)
length(unique(gbs_analysis$grid_reference)) # 782 sites 

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 649 gardens 

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 647 gardens 

# subset analysis to only sites with long grass present 
gbs_analysis2 <- gbs_analysis[gbs_analysis$garden_long_grass==1,]
length(unique(gbs_analysis2$grid_reference)) # 284

## Abundance model

abund_area_garden <- lmer(log(total_abund) ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis2) 
summary(abund_area_garden) 

car::vif(abund_area_garden) # all <3
r.squaredGLMM(abund_area_garden) # 6.8%
AIC(abund_area_garden2) # 1269.403

testDispersion(abund_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_area_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# Plot graph
abund_long_grass <- ggpredict(abund_area_garden, terms="garden_long_grass_area")

abund_long_grass_p2 <- ggplot() +
  geom_point(data=gbs_analysis2, aes(x = garden_long_grass_area, y=total_abund), colour="lightgrey", alpha=0.4) +
  geom_line(data = abund_long_grass, aes(x = x, y = predicted), lwd=1) +
  geom_ribbon(data = abund_long_grass, 
              aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low), alpha=0.1) +
  labs(x="Area of long grass", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_long_grass_p2
ggsave(abund_long_grass_p2, file="Graphs/Abundance/Abund_area_long_grass.png", height=4, width=4)


## Richness model

rich_area_garden <- lmer(site_rel_SR ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis2) 
summary(rich_area_garden) 

car::vif(rich_area_garden) # all <3
r.squaredGLMM(rich_area_garden) # 8.4%
AIC(rich_area_garden) # -1818.518

testDispersion(rich_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_area_garden, plot = F)
plot(simulationOutput) # quantile deviations

rich_long_grass <- ggpredict(rich_area_garden, terms="garden_long_grass_area")

rich_long_grass_p2 <- ggplot() +
  geom_point(data=gbs_analysis2, aes(x = garden_long_grass_area, y=site_rel_SR), colour="lightgrey", alpha=0.4) +
  geom_line(data = rich_long_grass, aes(x = x, y = predicted), lwd=1) +
  geom_ribbon(data = rich_long_grass, 
              aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low), alpha=0.1) +
  labs(x="Area of long grass", y="Relative richness")+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
rich_long_grass_p2
ggsave(rich_long_grass_p2, file="Graphs/Richness/Rich_area_long_grass.png", height=4, width=4)

## Put garden plots together
library(ggpubr)
garden_features <- ggarrange(abund_long_grass_p, rich_long_grass_p, abund_long_grass_p2, rich_long_grass_p2, 
                             labels = c("(a)", "(b)", "(c)", "(d)"),
                             ncol = 2, nrow = 2)
garden_features
ggsave(garden_features, file="Graphs/Garden_features_plot_Fig3.png", height=7, width=8)

####### 2b: grassland species, all sites ####### 

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_grass.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from area of long grass and garden size
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass_area, garden_size)
length(unique(gbs_analysis$grid_reference)) # 761

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 630 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 628 gardens 

# subset analysis to only sites with long grass present 
gbs_analysis2 <- gbs_analysis[gbs_analysis$garden_long_grass==1,]
length(unique(gbs_analysis2$grid_reference)) # 278

## Abundance model

abund_area_garden <- lmer(log(total_abund) ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis2) 
summary(abund_area_garden) 

car::vif(abund_area_garden) # all <3
r.squaredGLMM(abund_area_garden) # 6.8%
AIC(abund_area_garden2) # 1269.403

testDispersion(abund_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_area_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

## Richness model

rich_area_garden <- lmer(site_rel_SR ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis2) 
summary(rich_area_garden) 

car::vif(rich_area_garden) # all <3
r.squaredGLMM(rich_area_garden) # 8.4%
AIC(rich_area_garden) # -1818.518

testDispersion(rich_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_area_garden, plot = F)
plot(simulationOutput) # quantile deviations

####### 2c: non-grassland species, all sites #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_not_grass.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from area of long grass and garden size
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass_area, garden_size)
length(unique(gbs_analysis$grid_reference)) # 783

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 649 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 647 gardens 

# subset analysis to only sites with long grass present 
gbs_analysis2 <- gbs_analysis[gbs_analysis$garden_long_grass==1,]
length(unique(gbs_analysis2$grid_reference)) # 284

## Abundance model

abund_area_garden <- lmer(log(total_abund) ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis2) 
summary(abund_area_garden) 

car::vif(abund_area_garden) # all <3
r.squaredGLMM(abund_area_garden) # 6.8%
AIC(abund_area_garden2) # 1269.403

testDispersion(abund_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_area_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

## Richness model

rich_area_garden <- lmer(site_rel_SR ~ garden_size+ n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis2) 
summary(rich_area_garden) 

car::vif(rich_area_garden) # all <3
r.squaredGLMM(rich_area_garden) # 8.4%
AIC(rich_area_garden) # -1818.518

testDispersion(rich_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_area_garden, plot = F)
plot(simulationOutput) # quantile deviations


####### 2d: all species, climate zone sites #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_GEnS.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from area of long grass and garden size
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass_area, garden_size)
length(unique(gbs_analysis$grid_reference)) # 748

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 623 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 621 gardens 

# subset analysis to only sites with long grass present 
gbs_analysis2 <- gbs_analysis[gbs_analysis$garden_long_grass==1,]
length(unique(gbs_analysis2$grid_reference)) # 273

## Abundance model

abund_area_garden <- lmer(log(total_abund) ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis2) 
summary(abund_area_garden) 

car::vif(abund_area_garden) # all <3
r.squaredGLMM(abund_area_garden) # 6.8%
AIC(abund_area_garden2) # 1269.403

testDispersion(abund_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_area_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

## Richness model

rich_area_garden <- lmer(site_rel_SR ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis2) 
summary(rich_area_garden) 

car::vif(rich_area_garden) # all <3
r.squaredGLMM(rich_area_garden) # 8.4%
AIC(rich_area_garden) # -1818.518

testDispersion(rich_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_area_garden, plot = F)
plot(simulationOutput) # quantile deviations

####### 2e: all species, top 50% sites #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_top_rec.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from area of long grass and garden size
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_long_grass_area, garden_size)
length(unique(gbs_analysis$grid_reference)) # 382

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_long_grass_area <- as.numeric(gbs_analysis$garden_long_grass_area)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==2, 0, 1) # 0 = no long grass, 1 = long grass
gbs_analysis$garden_long_grass <- as.factor(gbs_analysis$garden_long_grass)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 305 gardens

# remove gardens with area long grass >450 (this is the largest a garden can be in category 3)
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_long_grass_area>=450,]
length(unique(gbs_analysis$grid_reference)) # 305 gardens 

# subset analysis to only sites with long grass present 
gbs_analysis2 <- gbs_analysis[gbs_analysis$garden_long_grass==1,]
length(unique(gbs_analysis2$grid_reference)) # 131

## Abundance model

abund_area_garden <- lmer(log(total_abund) ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                          na.action = "na.fail", data=gbs_analysis2) 
summary(abund_area_garden) 

car::vif(abund_area_garden) # all <3
r.squaredGLMM(abund_area_garden) # 6.8%
AIC(abund_area_garden2) # 1269.403

testDispersion(abund_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_area_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

## Richness model

rich_area_garden <- lmer(site_rel_SR ~ garden_size + n_days + scale(garden_long_grass_area) + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis2) 
summary(rich_area_garden) 

car::vif(rich_area_garden) # all <3
r.squaredGLMM(rich_area_garden) # 8.4%
AIC(rich_area_garden) # -1818.518

testDispersion(rich_area_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_area_garden, plot = F)
plot(simulationOutput) # quantile deviations


########################################################
########## Garden practice 2: Presence of ivy ########## 
########################################################

####### 3a: all species, all sites in autumn ####### 

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_ivy.csv", header=TRUE)

## remove NAs from garden features 
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_ivy, garden_size)
length(unique(gbs_analysis$grid_reference)) # 724 sites 

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==2, 0, 1) # 0 = no ivy, 1 = ivy
gbs_analysis$garden_ivy <- as.factor(gbs_analysis$garden_ivy)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 581

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_ivy_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(abund_ivy_garden) # ivy non-significant 

car::vif(abund_ivy_garden) # all <3
r.squaredGLMM(abund_ivy_garden) 

testDispersion(abund_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

## Richness model

rich_ivy_garden <- lmer(site_rel_SR ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                        na.action = "na.fail", data=gbs_analysis) 
summary(rich_ivy_garden) # garden size significant, ivy is not

car::vif(rich_ivy_garden) # all <3
r.squaredGLMM(rich_ivy_garden) 

testDispersion(rich_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 


####### 3b: all species, climate zone sites in autumn #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_ivy_GEnS.csv", header=TRUE)

## remove NAs from garden features 
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_ivy, garden_size)
length(unique(gbs_analysis$grid_reference)) # 572 sites 

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==2, 0, 1) # 0 = no ivy, 1 = ivy
gbs_analysis$garden_ivy <- as.factor(gbs_analysis$garden_ivy)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 461

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_ivy_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(abund_ivy_garden) # ivy non-significant 

car::vif(abund_ivy_garden) # all <3
r.squaredGLMM(abund_ivy_garden) 

testDispersion(abund_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

## Richness model

rich_ivy_garden <- lmer(site_rel_SR ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                        na.action = "na.fail", data=gbs_analysis) 
summary(rich_ivy_garden) # garden size significant, ivy is not

car::vif(rich_ivy_garden) # all <3
r.squaredGLMM(rich_ivy_garden) 

testDispersion(rich_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

####### 3c: all species, top 50% sites in autumn #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_ivy_top_rec.csv", header=TRUE)

## remove NAs from garden features 
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_ivy, garden_size)
length(unique(gbs_analysis$grid_reference)) # 374 sites 

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==2, 0, 1) # 0 = no ivy, 1 = ivy
gbs_analysis$garden_ivy <- as.factor(gbs_analysis$garden_ivy)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 301

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_ivy_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(abund_ivy_garden) # ivy non-significant 

car::vif(abund_ivy_garden) # all <3
r.squaredGLMM(abund_ivy_garden) 

testDispersion(abund_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

## Richness model

rich_ivy_garden <- lmer(site_rel_SR ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                        na.action = "na.fail", data=gbs_analysis) 
summary(rich_ivy_garden) # garden size significant, ivy is not

car::vif(rich_ivy_garden) # all <3
r.squaredGLMM(rich_ivy_garden) 

testDispersion(rich_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 


####### 3d: Red Admiral and Comma combined, all sites #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_ivy_Comma_Red_Admiral.csv", header=TRUE)

## remove NAs from garden features 
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_ivy, garden_size)
length(unique(gbs_analysis$grid_reference)) # 710 sites 

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==2, 0, 1) # 0 = no ivy, 1 = ivy
gbs_analysis$garden_ivy <- as.factor(gbs_analysis$garden_ivy)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 570

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_ivy_garden <- lmer(log(total_abund) ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(abund_ivy_garden) # ivy non-significant 

car::vif(abund_ivy_garden) # all <3
r.squaredGLMM(abund_ivy_garden) 

testDispersion(abund_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# plot result
abund_ivy_ra_c <- ggpredict(abund_ivy_garden, terms="garden_ivy")
gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==0, "Absent", "Present") # 0 = no long grass, 1 = long grass
abund_ivy_ra_c$x <- ifelse(abund_ivy_ra_c$x==0, "Absent", "Present") # 0 = no long grass, 1 = long grass

abund_ivy_ra_c_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x = garden_ivy, y=total_abund), colour="lightgrey", alpha=0.4,
             position=position_jitter(w = 0.2, h = 0)) +
  geom_point(data = abund_ivy_ra_c, aes(x = x, y = predicted), 
             size = 1,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(data = abund_ivy_ra_c, 
                aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low), 
                position = position_dodge(width = 0.4),
                width = 0.1) +
  labs(x="Flowering ivy", y="Total autumn abundance of \nRed Admiral and Comma")+
  scale_y_continuous(trans = 'log', breaks=c(2,10,50,150,300)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_ivy_ra_c_p

####### 3e: Holly blue, all sites in summer (2nd gen) #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_HollyBlue.csv", header=TRUE)

## remove NAs from garden features 
gbs_analysis <- gbs_analysis %>% 
  drop_na(garden_ivy, garden_size)
length(unique(gbs_analysis$grid_reference)) # 628 sites 

## make sure variables are in correct format
str(gbs_analysis)
gbs_analysis$M_YEAR <- as.factor(gbs_analysis$M_YEAR)
gbs_analysis$garden_size <- as.integer(gbs_analysis$garden_size)
gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==2, 0, 1) # 0 = no ivy, 1 = ivy
gbs_analysis$garden_ivy <- as.factor(gbs_analysis$garden_ivy)

# Remove gardens in size 4 category - these are skewing the results as gardens are v large
gbs_analysis <- gbs_analysis[!gbs_analysis$garden_size==4,]
length(unique(gbs_analysis$grid_reference)) # 527 sites

## Abundance model

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

abund_ivy_garden <- lmer(log(SINDEX) ~ garden_size + n_days + garden_ivy + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis) 
summary(abund_ivy_garden) # ivy non-significant 

car::vif(abund_ivy_garden) # all <3
r.squaredGLMM(abund_ivy_garden) 

testDispersion(abund_ivy_garden) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_ivy_garden, plot = F)
plot(simulationOutput) # no assumptions violated 

# plot result
abund_ivy_hb <- ggpredict(abund_ivy_garden, terms="garden_ivy")
gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==0, "Absent", "Present") # 0 = no long grass, 1 = long grass
abund_ivy_hb$x <- ifelse(abund_ivy_hb$x==0, "Absent", "Present") # 0 = no long grass, 1 = long grass

abund_ivy_hb_p <- ggplot() +
  geom_point(data=gbs_analysis, aes(x = garden_ivy, y=SINDEX), colour="lightgrey", alpha=0.4,
             position=position_jitter(w = 0.2, h = 0)) +
  geom_point(data = abund_ivy_hb, aes(x = x, y = predicted), 
             size = 1,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(data = abund_ivy_hb, 
                aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low), 
                position = position_dodge(width = 0.4),
                width = 0.1) +
  labs(x="Flowering ivy", y="Total abundance of summer \ngeneration Holly Blue")+
  scale_y_continuous(trans = 'log', breaks=c(2,10,50)) +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_ivy_hb_p

# put the two ivy graphs together 
library(ggpubr)
garden_features <- ggarrange(abund_ivy_ra_c_p, abund_ivy_hb_p, 
                             labels = c("(a)", "(b)"),
                             ncol = 2, nrow = 1)
garden_features
ggsave(garden_features, file="Graphs/Garden_features_plot_Fig4.png", height=5, width=8)

