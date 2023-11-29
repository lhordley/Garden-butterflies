##########################
#### user: Lisbeth Hordley
#### date: November 2022
#### info: Interactions between garden management and landscape 
options(scipen = 100)

library(ggplot2)
library(broom)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)
library(lme4)
library(lmerTest)
library(DHARMa)
library(ggeffects)
library(MuMIn)

# First want to check whether there are enough sites across a gradient of landscape in presence and absence of long grass

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


# split area and landscape into 4 categories and count number of sites

gbs_analysis$Urban_100m2 <- as.factor(ifelse(gbs_analysis$Urban_100m<=0.25, "0-0.25",
                                              ifelse(gbs_analysis$Urban_100m>0.25 & gbs_analysis$Urban_100m<=0.5, "0.25-0.5",
                                                     ifelse(gbs_analysis$Urban_100m>0.5 & gbs_analysis$Urban_100m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Arable_100m2 <- as.factor(ifelse(gbs_analysis$Arable_100m<=0.25, "0-0.25",
                                               ifelse(gbs_analysis$Arable_100m>0.25 & gbs_analysis$Arable_100m<=0.5, "0.25-0.5",
                                                      ifelse(gbs_analysis$Arable_100m>0.5 & gbs_analysis$Arable_100m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Woodland_100m2 <- as.factor(ifelse(gbs_analysis$Woodland_100m<=0.25, "0-0.25",
                                                 ifelse(gbs_analysis$Woodland_100m>0.25 & gbs_analysis$Woodland_100m<=0.5, "0.25-0.5",
                                                        ifelse(gbs_analysis$Woodland_100m>0.5 & gbs_analysis$Woodland_100m<=0.75, "0.5-0.75", "0.75-1"))))

gbs_analysis$Urban_250m2 <- as.factor(ifelse(gbs_analysis$Urban_250m<=0.25, "0-0.25",
                                              ifelse(gbs_analysis$Urban_250m>0.25 & gbs_analysis$Urban_250m<=0.5, "0.25-0.5",
                                                     ifelse(gbs_analysis$Urban_250m>0.5 & gbs_analysis$Urban_250m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Arable_250m2 <- as.factor(ifelse(gbs_analysis$Arable_250m<=0.25, "0-0.25",
                                               ifelse(gbs_analysis$Arable_250m>0.25 & gbs_analysis$Arable_250m<=0.5, "0.25-0.5",
                                                      ifelse(gbs_analysis$Arable_250m>0.5 & gbs_analysis$Arable_250m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Woodland_250m2 <- as.factor(ifelse(gbs_analysis$Woodland_250m<=0.25, "0-0.25",
                                                 ifelse(gbs_analysis$Woodland_250m>0.25 & gbs_analysis$Woodland_250m<=0.5, "0.25-0.5",
                                                        ifelse(gbs_analysis$Woodland_250m>0.5 & gbs_analysis$Woodland_250m<=0.75, "0.5-0.75", "0.75-1"))))


gbs_analysis$Urban_500m2 <- as.factor(ifelse(gbs_analysis$Urban_500m<=0.25, "0-0.25",
                                              ifelse(gbs_analysis$Urban_500m>0.25 & gbs_analysis$Urban_500m<=0.5, "0.25-0.5",
                                                     ifelse(gbs_analysis$Urban_500m>0.5 & gbs_analysis$Urban_500m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Arable_500m2 <- as.factor(ifelse(gbs_analysis$Arable_500m<=0.25, "0-0.25",
                                               ifelse(gbs_analysis$Arable_500m>0.25 & gbs_analysis$Arable_500m<=0.5, "0.25-0.5",
                                                      ifelse(gbs_analysis$Arable_500m>0.5 & gbs_analysis$Arable_500m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Woodland_500m2 <- as.factor(ifelse(gbs_analysis$Woodland_500m<=0.25, "0-0.25",
                                                 ifelse(gbs_analysis$Woodland_500m>0.25 & gbs_analysis$Woodland_500m<=0.5, "0.25-0.5",
                                                        ifelse(gbs_analysis$Woodland_500m>0.5 & gbs_analysis$Woodland_500m<=0.75, "0.5-0.75", "0.75-1"))))

gbs_analysis2 <- gbs_analysis[,c(1,7,26:34)] # grid reference, presence of long grass and new landscape variables
gbs_analysis3 <- melt(setDT(gbs_analysis2), id.vars = c("grid_reference","garden_long_grass"), variable.name = "landscape")
dat_presence <- gbs_analysis3 %>%                    
  group_by(landscape, value, garden_long_grass) %>%          
  summarise(unique_sites = n_distinct(grid_reference)) 
dat_presence$garden_long_grass  <- factor(dat_presence$garden_long_grass)

land_lookup <- c(Urban_100m2="Urban 100m", Arable_100m2="Arable 100m", Woodland_100m2="Woodland 100m", 
                 Urban_250m2="Urban 250m", Arable_250m2="Arable 250m", Woodland_250m2="Woodland 250m",
                 Urban_500m2="Urban 500m", Arable_500m2="Arable 500m", Woodland_500m2="Woodland 500m")
dat_presence$landscape <- as.character(land_lookup[dat_presence$landscape])

dat_presence$garden_long_grass <- ifelse(dat_presence$garden_long_grass==0, "Absent", "Present") # 0 = no long grass, 1 = long grass

no_sites_presence <- ggplot(data=dat_presence, aes(x=value, y=unique_sites, fill=garden_long_grass)) +
  geom_bar(stat='identity') +
  facet_grid(.~landscape)+
  labs(x="Proportion of landscape in buffer", y="Number of sites")+
  guides(fill=guide_legend(title="Long grass"))+
  scale_y_continuous(breaks=seq(0,800,by=100))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size=24), legend.title=element_text(size=24), 
        legend.text=element_text(size=24))
no_sites_presence
#

ggplot(gbs_analysis, aes(x = Woodland_500m, y = log(total_abund))) +
  geom_point(size = .1) +
  stat_smooth(method = "lm", se = FALSE)+
  facet_wrap(~ garden_long_grass_area)

coplot(log(total_abund) ~ Arable_250m | garden_long_grass_area, data = gbs_analysis)


## Same again for area of long grass

# subset analysis to only sites with long grass present 
gbs_analysis_area <- gbs_analysis[gbs_analysis$garden_long_grass==1,]
length(unique(gbs_analysis_area$grid_reference)) # 284

quantile(gbs_analysis_area$garden_long_grass_area, 0.05)
quantile(gbs_analysis_area$garden_long_grass_area, 0.25)
quantile(gbs_analysis_area$garden_long_grass_area, 0.5)
quantile(gbs_analysis_area$garden_long_grass_area, 0.75)
quantile(gbs_analysis_area$garden_long_grass_area, 0.95)

gbs_analysis_area$garden_long_grass_area2 <- as.factor(ifelse(gbs_analysis_area$garden_long_grass_area<=1, "<2m2",
                                                          ifelse(gbs_analysis_area$garden_long_grass_area>1 & gbs_analysis_area$garden_long_grass_area<=3, "2-3m2",
                                                                 ifelse(gbs_analysis_area$garden_long_grass_area>3 & gbs_analysis_area$garden_long_grass_area<=8, "4-8m2",
                                                                        ifelse(gbs_analysis_area$garden_long_grass_area>8 & gbs_analysis_area$garden_long_grass_area<=20, "9-20m2",
                                                                               ifelse(gbs_analysis_area$garden_long_grass_area>20 & gbs_analysis_area$garden_long_grass_area<=100, "21-100m2", "101-400m2"))))))

gbs_analysis2 <- gbs_analysis_area[,c(1,26:35)] # grid reference, new area of long grass and new landscape variables
gbs_analysis3 <- melt(setDT(gbs_analysis2), id.vars = c("grid_reference","garden_long_grass_area2"), variable.name = "landscape")
dat_area <- gbs_analysis3 %>%                    
  group_by(landscape, value, garden_long_grass_area2) %>%          
  summarise(unique_sites = n_distinct(grid_reference)) 
dat_area$garden_long_grass_area2  <- factor(dat_area$garden_long_grass_area2,levels = c("<2m2", "2-3m2", "4-8m2", "9-20m2", "21-100m2", "101-400m2"))

land_lookup <- c(Urban_100m2="Urban 100m", Arable_100m2="Arable 100m", Woodland_100m2="Woodland 100m", 
                 Urban_250m2="Urban 250m", Arable_250m2="Arable 250m", Woodland_250m2="Woodland 250m",
                 Urban_500m2="Urban 500m", Arable_500m2="Arable 500m", Woodland_500m2="Woodland 500m")
dat_area$landscape <- as.character(land_lookup[dat_area$landscape])

my_leg <- c(expression(paste("<2m"^"2")),
           expression(paste("2-3m"^"2")),
           expression(paste("4-8m"^"2")), 
           expression(paste("9-20m"^"2")), 
           expression(paste("21-100m"^"2")), 
           expression(paste("101-400m"^"2")))

no_sites_area <- ggplot(data=dat_area, aes(x=value, y=unique_sites, fill=garden_long_grass_area2)) +
  geom_bar(stat='identity') +
  facet_grid(.~landscape)+
  labs(x="Proportion of landscape in buffer", y="Number of sites")+
  guides(fill=guide_legend(title="Area of long grass"))+
  scale_y_continuous(breaks=seq(0,300,by=50))+
  scale_fill_discrete(name="Area of long grass",
                      labels=c(my_leg[1], 
                               my_leg[2],
                               my_leg[3],
                               my_leg[4],
                               my_leg[5],
                               my_leg[6]))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size=24), legend.title=element_text(size=24), 
        legend.text=element_text(size=24), legend.text.align = 0)
no_sites_area


## repeat for ivy (all species)

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_ivy.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden presence of long grass and garden size
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
length(unique(gbs_analysis$grid_reference)) # 581 gardens 

# split area and landscape into 4 categories and count number of sites

gbs_analysis$Urban_100m2 <- as.factor(ifelse(gbs_analysis$Urban_100m<=0.25, "0-0.25",
                                             ifelse(gbs_analysis$Urban_100m>0.25 & gbs_analysis$Urban_100m<=0.5, "0.25-0.5",
                                                    ifelse(gbs_analysis$Urban_100m>0.5 & gbs_analysis$Urban_100m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Arable_100m2 <- as.factor(ifelse(gbs_analysis$Arable_100m<=0.25, "0-0.25",
                                              ifelse(gbs_analysis$Arable_100m>0.25 & gbs_analysis$Arable_100m<=0.5, "0.25-0.5",
                                                     ifelse(gbs_analysis$Arable_100m>0.5 & gbs_analysis$Arable_100m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Woodland_100m2 <- as.factor(ifelse(gbs_analysis$Woodland_100m<=0.25, "0-0.25",
                                                ifelse(gbs_analysis$Woodland_100m>0.25 & gbs_analysis$Woodland_100m<=0.5, "0.25-0.5",
                                                       ifelse(gbs_analysis$Woodland_100m>0.5 & gbs_analysis$Woodland_100m<=0.75, "0.5-0.75", "0.75-1"))))

gbs_analysis$Urban_250m2 <- as.factor(ifelse(gbs_analysis$Urban_250m<=0.25, "0-0.25",
                                             ifelse(gbs_analysis$Urban_250m>0.25 & gbs_analysis$Urban_250m<=0.5, "0.25-0.5",
                                                    ifelse(gbs_analysis$Urban_250m>0.5 & gbs_analysis$Urban_250m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Arable_250m2 <- as.factor(ifelse(gbs_analysis$Arable_250m<=0.25, "0-0.25",
                                              ifelse(gbs_analysis$Arable_250m>0.25 & gbs_analysis$Arable_250m<=0.5, "0.25-0.5",
                                                     ifelse(gbs_analysis$Arable_250m>0.5 & gbs_analysis$Arable_250m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Woodland_250m2 <- as.factor(ifelse(gbs_analysis$Woodland_250m<=0.25, "0-0.25",
                                                ifelse(gbs_analysis$Woodland_250m>0.25 & gbs_analysis$Woodland_250m<=0.5, "0.25-0.5",
                                                       ifelse(gbs_analysis$Woodland_250m>0.5 & gbs_analysis$Woodland_250m<=0.75, "0.5-0.75", "0.75-1"))))

gbs_analysis$Urban_500m2 <- as.factor(ifelse(gbs_analysis$Urban_500m<=0.25, "0-0.25",
                                             ifelse(gbs_analysis$Urban_500m>0.25 & gbs_analysis$Urban_500m<=0.5, "0.25-0.5",
                                                    ifelse(gbs_analysis$Urban_500m>0.5 & gbs_analysis$Urban_500m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Arable_500m2 <- as.factor(ifelse(gbs_analysis$Arable_500m<=0.25, "0-0.25",
                                              ifelse(gbs_analysis$Arable_500m>0.25 & gbs_analysis$Arable_500m<=0.5, "0.25-0.5",
                                                     ifelse(gbs_analysis$Arable_500m>0.5 & gbs_analysis$Arable_500m<=0.75, "0.5-0.75", "0.75-1"))))
gbs_analysis$Woodland_500m2 <- as.factor(ifelse(gbs_analysis$Woodland_500m<=0.25, "0-0.25",
                                                ifelse(gbs_analysis$Woodland_500m>0.25 & gbs_analysis$Woodland_500m<=0.5, "0.25-0.5",
                                                       ifelse(gbs_analysis$Woodland_500m>0.5 & gbs_analysis$Woodland_500m<=0.75, "0.5-0.75", "0.75-1"))))

gbs_analysis2 <- gbs_analysis[,c(1,6,26:34)] # grid reference, presence of ivy and new landscape variables
gbs_analysis3 <- melt(setDT(gbs_analysis2), id.vars = c("grid_reference","garden_ivy"), variable.name = "landscape")
dat_ivy <- gbs_analysis3 %>%                    
  group_by(landscape, value, garden_ivy) %>%          
  summarise(unique_sites = n_distinct(grid_reference)) 
dat_ivy$garden_ivy  <- factor(dat_ivy$garden_ivy)

land_lookup <- c(Urban_100m2="Urban 100m", Arable_100m2="Arable 100m", Woodland_100m2="Woodland 100m", 
                 Urban_250m2="Urban 250m", Arable_250m2="Arable 250m", Woodland_250m2="Woodland 250m",
                 Urban_500m2="Urban 500m", Arable_500m2="Arable 500m", Woodland_500m2="Woodland 500m")
dat_ivy$landscape <- as.character(land_lookup[dat_ivy$landscape])
dat_ivy$garden_ivy <- ifelse(dat_ivy$garden_ivy==0, "Absent", "Present") # 0 = no long grass, 1 = long grass

no_sites_ivy <- ggplot(data=dat_ivy, aes(x=value, y=unique_sites, fill=garden_ivy)) +
  geom_bar(stat='identity') +
  facet_grid(.~landscape)+
  labs(x="Proportion of landscape in buffer", y="Number of sites")+
  guides(fill=guide_legend(title="Flowering ivy"))+
  scale_y_continuous(breaks=seq(0,800,by=100))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size=24), legend.title=element_text(size=24), 
        legend.text=element_text(size=24))
no_sites_ivy


## Conclusion: run interactions with presence of long grass and presence of flowering ivy
## but not area of long grass - not enough data across gradients of land use 

# but for interactions with woodland, make sure we are not predicting beyond the data - clip the prediction line if necessary

# put the 3 graphs together in the supplementary material 

library(ggpubr)
plot <- ggarrange(no_sites_presence, no_sites_area, no_sites_ivy, labels=c("(a)", "(b)", "(c)"), nrow=3, ncol=1, font.label=list(color="black",size=24))
plot
ggsave(plot, file="Graphs/FigureS2.png", height=15, width=24)


#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################


###### 1. Presence of long grass #######


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

gbs_analysis$garden_long_grass <- ifelse(gbs_analysis$garden_long_grass==0, "Absent", "Present") # 0 = no long grass, 1 = long grass

## ABUNDANCE 100m ##

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_100m", "Arable_100m", "Grassland_100m", "Woodland_100m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_100m) # 0.28
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_100m) # -0.35
# remove grassland from analysis 

abund_100_mod <- lmer(log(total_abund) ~ garden_size + n_days + garden_long_grass*scale(Urban_100m) +
                           garden_long_grass*scale(Arable_100m) + scale(woodland_dist_m) + scale(Woodland_100m) +
                           scale(I(Woodland_100m^2)) + scale(grassland_dist_m) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis)
summary(abund_100_mod)

# backwards model selection
abund_100_mod_step <- lmerTest::step(abund_100_mod, reduce.random=FALSE, alpha.fixed=0.05)
abund_100_mod_final <- get_model(abund_100_mod_step)
summary(abund_100_mod_final)
# arable significant interaction

car::vif(abund_100_mod_final) 
r.squaredGLMM(abund_100_mod_final) # 14%
AIC(abund_100_mod_final) # 2770.12

testDispersion(abund_100_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_100_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_100_summary <- summary(abund_100_mod_final) #pulling out model averages
abund_100_summary<-as.data.frame(abund_100_summary$coefficients) #selecting full model coefficient averages (more conservative)
abund_100_summary$parameters <- row.names(abund_100_summary)
row.names(abund_100_summary) <- 1:nrow(abund_100_summary)
abund_100_summary$AIC <- AIC(abund_100_mod_final)
write.csv(abund_100_summary, file="Results/Abundance/Abundance_100m_presence.csv", row.names=FALSE)

# plot significant interactions

# Arable
abund_arable_100 <- ggpredict(abund_100_mod_final, terms=c("Arable_100m","garden_long_grass")) 
abund_arable_100_p <- ggplot() +
  geom_line(data=abund_arable_100, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=abund_arable_100, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Arable_100m, y=total_abund, colour=garden_long_grass), alpha=0.2)+
  labs(x="Proportion of arable \nin 100m buffer", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()
abund_arable_100_p
ggsave(abund_arable_100_p, file="Graphs/Interactions/Abundance_presence_arable_100m.png")
# abundance increases as arable increases when gardens have long grass
# but abundance declines as arable increases when gardens don't have long grass


## ABUNDANCE 250m ##

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_250m", "Arable_250m", "Grassland_250m", "Woodland_250m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.798
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_250m) # 0.30
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_250m) # -0.37
# remove grassland from analysis 

abund_250_mod <- lmer(log(total_abund) ~ garden_size + n_days + scale(Woodland_250m) + garden_long_grass*scale(Urban_250m) +
                        garden_long_grass*scale(Arable_250m) + scale(woodland_dist_m) + scale(grassland_dist_m) + (1|M_YEAR), 
                         na.action = "na.fail", data=gbs_analysis)
summary(abund_250_mod)

# backwards model selection
abund_250_mod_step <- lmerTest::step(abund_250_mod, reduce.random=FALSE, alpha.fixed=0.05)
abund_250_mod_final <- get_model(abund_250_mod_step)
summary(abund_250_mod_final)
# arable and urban significant interaction

car::vif(abund_250_mod_final) 
AIC(abund_250_mod_final) # 2723.979

testDispersion(abund_250_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_250_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_250_summary <- summary(abund_250_mod_final) #pulling out model averages
abund_250_summary<-as.data.frame(abund_250_summary$coefficients) 
abund_250_summary$parameters <- row.names(abund_250_summary)
row.names(abund_250_summary) <- 1:nrow(abund_250_summary)
abund_250_summary$AIC <- AIC(abund_250_mod_final)
write.csv(abund_250_summary, file="Results/Abundance/Abundance_250m_presence.csv", row.names=FALSE)

# plot significant interactions

# Arable
abund_arable_250 <- ggpredict(abund_250_mod_final, terms=c("Arable_250m","garden_long_grass")) 
abund_arable_250_p <- ggplot() +
  geom_line(data=abund_arable_250, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=abund_arable_250, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Arable_250m, y=total_abund, colour=garden_long_grass), alpha=0.1)+
  labs(x="Proportion of arable \nin 250m buffer", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()+
  theme(text = element_text(size=20), legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
abund_arable_250_p
ggsave(abund_arable_250_p, file="Graphs/Interactions/Abundance_presence_arable_250m.png", height=6, width=7)
# abundance increases as arable increases when gardens have long grass
# but abundance declines as arable increases when gardens don't have long grass
# same as 100m relationship
library(emmeans)
emtrends(abund_250_mod_final, ~ garden_long_grass, var="Arable_250m")

# Urban
abund_urban_250 <- ggpredict(abund_250_mod_final, terms=c("Urban_250m","garden_long_grass")) 
abund_urban_250_p <- ggplot() +
  geom_line(data=abund_urban_250, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=abund_urban_250, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Urban_250m, y=total_abund, colour=garden_long_grass), alpha=0.1)+
  labs(x="Proportion of urban \nin 250m buffer", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
abund_urban_250_p
ggsave(abund_urban_250_p, file="Graphs/Interactions/Abundance_presence_urban_250m.png")
# decline in abundance across urban gradient is less steep in gardens with long grass
# compared to gardens without long grass
emtrends(abund_250_mod_final, ~ garden_long_grass, var="Urban_250m")


## ABUNDANCE 500m ##

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_500m", "Arable_500m", "Grassland_500m", "Woodland_500m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.763
# threshold of 0.7 - need to remove with urban OR grassland
# urban and grassland distance are close at -0.68

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_500m) # 0.27
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_500m) # -0.35
# remove grassland from analysis 

abund_500_mod <- lmer(log(total_abund) ~ garden_size + n_days + scale(Woodland_500m) + garden_long_grass*scale(Urban_500m) +
                         garden_long_grass*scale(Arable_500m) + scale(woodland_dist_m) + scale(grassland_dist_m) + (1|M_YEAR), 
                       na.action = "na.fail", data=gbs_analysis)
summary(abund_500_mod)

# backwards model selection
abund_500_mod_step <- lmerTest::step(abund_500_mod, reduce.random=FALSE, alpha.fixed=0.05)
abund_500_mod_final <- get_model(abund_500_mod_step)
summary(abund_500_mod_final)
# arable and urban significant interaction (linear urban only - quadratic interaction removed)

car::vif(abund_500_mod_final) 
AIC(abund_500_mod_final) # 2724.31

testDispersion(abund_500_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_500_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_500_summary <- summary(abund_500_mod_final) #pulling out model coefficients
abund_500_summary<-as.data.frame(abund_500_summary$coefficients) 
abund_500_summary$parameters <- row.names(abund_500_summary)
row.names(abund_500_summary) <- 1:nrow(abund_500_summary)
abund_500_summary$AIC <- AIC(abund_500_mod_final)
write.csv(abund_500_summary, file="Results/Abundance/Abundance_500m_presence.csv", row.names=FALSE)

# plot significant interactions

# Arable
abund_arable_500 <- ggpredict(abund_500_mod_final, terms=c("Arable_500m [all]","garden_long_grass")) 
abund_arable_500_p <- ggplot() +
  geom_line(data=abund_arable_500, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=abund_arable_500, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Arable_500m, y=total_abund, colour=garden_long_grass), alpha=0.2)+
  labs(x="Proportion of arable \nin 500m buffer", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()
abund_arable_500_p
ggsave(abund_arable_500_p, file="Graphs/Interactions/Abundance_presence_arable_500m.png")
# abundance increases as arable increases when gardens have long grass
# but no change in abundance as arable increases when gardens don't have long grass
# same as 100m and 250m relationships

# Urban
abund_urban_500 <- ggpredict(abund_500_mod_final, terms=c("Urban_500m [all]","garden_long_grass")) 
abund_urban_500_p <- ggplot() +
  geom_line(data=abund_urban_500, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=abund_urban_500, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Urban_500m, y=total_abund, colour=garden_long_grass), alpha=0.2)+
  labs(x="Proportion of urban \nin 500m buffer", y="Total abundance")+
  scale_y_continuous(trans = 'log', breaks=c(50,150,300,450,600)) +
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()
abund_urban_500_p
ggsave(abund_urban_500_p, file="Graphs/Interactions/Abundance_presence_urban_500m.png")
# this plots the quadratic interaction, but it's the linear one that's significant...
# but when the urban quadratic term is removed from model, the urban linear interaction
# plot looks identical to the 250m one that will be in the manuscript


## RICHNESS 100m ##
# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_100m", "Arable_100m", "Grassland_100m", "Woodland_100m", 
                          "Urban_vs_Grassland_100m", "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_100m) # 0.29
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_100m) # -0.38
# remove grassland from analysis 

rich_100_mod <- lmer(site_rel_SR ~ garden_size + n_days + garden_long_grass*scale(Arable_100m) + scale(Woodland_100m) + 
                        garden_long_grass*scale(Urban_100m) + scale(woodland_dist_m) + scale(grassland_dist_m) + (1|M_YEAR), 
                        na.action = "na.fail", data=gbs_analysis)
summary(rich_100_mod)

# backwards model selection
rich_100_mod_step <- lmerTest::step(rich_100_mod, reduce.random=FALSE, alpha.fixed=0.05)
rich_100_mod_final <- get_model(rich_100_mod_step)
summary(rich_100_mod_final)
# urban significant interaction

car::vif(rich_100_mod_final) 
AIC(rich_100_mod_final) # -4503.56

testDispersion(rich_100_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_100_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_100_summary <- summary(rich_100_mod_final) #pulling out model coefficients
rich_100_summary<-as.data.frame(rich_100_summary$coefficients) 
rich_100_summary$parameters <- row.names(rich_100_summary)
row.names(rich_100_summary) <- 1:nrow(rich_100_summary)
rich_100_summary$AIC <- AIC(rich_100_mod_final)
write.csv(rich_100_summary, file="Results/Richness/Richness_100m_presence.csv", row.names=FALSE)

# plot significant interactions

# Urban
rich_urban_100 <- ggpredict(rich_100_mod_final, terms=c("Urban_100m","garden_long_grass")) 
rich_urban_100_p <- ggplot() +
  geom_line(data=rich_urban_100, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=rich_urban_100, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Urban_100m, y=site_rel_SR, colour=garden_long_grass), alpha=0.2)+
  labs(x="Proportion of urban \nin 100m buffer", y="Relative richness")+
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()
rich_urban_100_p
ggsave(rich_urban_100_p, file="Graphs/Interactions/Richness_presence_urban_100m.png")
# decline in richness across urban gradient is less steep in gardens with long grass
# compared to gardens without long grass


## RICHNESS 250m ##
# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_250m", "Arable_250m", "Grassland_250m", "Woodland_250m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.798
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_250m) # 0.29
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_250m) # -0.41
# remove grassland from analysis 

rich_250_mod <- lmer(site_rel_SR ~ garden_size + n_days + scale(Woodland_250m) + garden_long_grass*scale(Urban_250m) +
                        garden_long_grass*scale(Arable_250m) + scale(woodland_dist_m) + scale(grassland_dist_m) + (1|M_YEAR), 
                      data=gbs_analysis)
summary(rich_250_mod)

# backwards model selection
rich_250_mod_step <- lmerTest::step(rich_250_mod, reduce.random=FALSE, alpha.fixed=0.05)
rich_250_mod_final <- get_model(rich_250_mod_step)
summary(rich_250_mod_final)
# urban and woodland significant interaction

car::vif(rich_250_mod_final) 
AIC(rich_250_mod_final) # -4527.63

testDispersion(rich_250_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_250_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_250_summary <- summary(rich_250_mod_final) #pulling out model coefficients
rich_250_summary<-as.data.frame(rich_250_summary$coefficients) 
rich_250_summary$parameters <- row.names(rich_250_summary)
row.names(rich_250_summary) <- 1:nrow(rich_250_summary)
rich_250_summary$AIC <- AIC(rich_250_mod_final)
write.csv(rich_250_summary, file="Results/Richness/Richness_250m_presence.csv", row.names=FALSE)

# plot significant interactions

# Urban
rich_urban_250 <- ggpredict(rich_250_mod_final, terms=c("Urban_250m","garden_long_grass")) 
rich_urban_250_p <- ggplot() +
  geom_line(data=rich_urban_250, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=rich_urban_250, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Urban_250m, y=site_rel_SR, colour=garden_long_grass), alpha=0.1)+
  labs(x="Proportion of urban \nin 250m buffer", y="Relative richness")+
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()+
  theme(text = element_text(size=20), legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
rich_urban_250_p
ggsave(rich_urban_250_p, file="Graphs/Interactions/Richness_presence_urban_250m.png", height=6, width=7)
# decline in richness across urban gradient is in gardens without long grass
# whereas richness relatively stable across urban gradient in gardens with long grass
emtrends(rich_250_mod_final, ~ garden_long_grass, var="Urban_250m")

## RICHNESS 500m ##
# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_500m", "Arable_500m", "Grassland_500m", "Woodland_500m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.763
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_500m) # 0.28
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_500m) # -0.41
# remove grassland from analysis 

rich_500_mod <- lmer(site_rel_SR ~ garden_size + n_days + scale(Woodland_500m) + garden_long_grass*scale(Urban_500m) +
                        garden_long_grass*scale(Arable_500m) + scale(woodland_dist_m) + scale(grassland_dist_m) + (1|M_YEAR), 
                      na.action = "na.fail", data=gbs_analysis)
summary(rich_500_mod)

# backwards model selection
rich_500_mod_step <- lmerTest::step(rich_500_mod, reduce.random=FALSE, alpha.fixed=0.05)
rich_500_mod_final <- get_model(rich_500_mod_step)
summary(rich_500_mod_final)
# urban and woodland significant interaction

car::vif(rich_500_mod_final) 
AIC(rich_500_mod_final) # -4521.777

testDispersion(rich_500_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_500_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_500_summary <- summary(rich_500_mod_final) #pulling out model coefficients
rich_500_summary<-as.data.frame(rich_500_summary$coefficients) 
rich_500_summary$parameters <- row.names(rich_500_summary)
row.names(rich_500_summary) <- 1:nrow(rich_500_summary)
rich_500_summary$AIC <- AIC(rich_500_mod_final)
write.csv(rich_500_summary, file="Results/Richness/Richness_500m_presence.csv", row.names=FALSE)

# plot significant interactions

# Urban
rich_urban_500 <- ggpredict(rich_500_mod_final, terms=c("Urban_500m","garden_long_grass")) 
rich_urban_500_p <- ggplot() +
  geom_line(data=rich_urban_500, aes(x=x, y=predicted, colour=group), lwd=1) +
  geom_ribbon(data=rich_urban_500, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=.5)+
  geom_point(data=gbs_analysis, aes(x=Urban_500m, y=site_rel_SR, colour=garden_long_grass), alpha=0.2)+
  labs(x="Proportion of urban \nin 500m buffer", y="Relative richness")+
  guides(fill=guide_legend(title="Long grass"), colour=guide_legend(title="Long grass"))+
  scale_fill_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  scale_colour_manual(values=c("Present"="deepskyblue3", "Absent"="coral2"))+
  theme_classic()
rich_urban_500_p
ggsave(rich_urban_500_p, file="Graphs/Interactions/Richness_presence_urban_500m.png")
# decline in richness across urban gradient is in gardens without long grass
# whereas richness relatively stable across urban gradient in gardens with long grass


## only put graphs in manuscript of lowest AIC models
# 250m for abundance (500m is equivalent though.. but technically slightly higher)
# 250m for richness
library(ggpubr)
plots <- ggarrange(abund_arable_250_p, abund_urban_250_p, rich_urban_250_p, labels=c("(a)", "(b)", "(c)"), nrow=1, ncol=3,
                   common.legend = TRUE, legend="right")
plots
ggsave(plots, file="Graphs/Figure5.png", height=6, width=15)





#################################################################################################
#################################################################################################


###### 2. Presence of flowering ivy #######

# load data
gbs_analysis <- read.csv("Data/GBS_analysis_final_ivy.csv", header=TRUE)

# remove extreme value
gbs_analysis <- gbs_analysis[!gbs_analysis$total_abund>=2000,]

## remove NAs from garden presence of long grass and garden size
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
length(unique(gbs_analysis$grid_reference)) # 581 gardens 

gbs_analysis$garden_ivy <- ifelse(gbs_analysis$garden_ivy==0, "Absent", "Present") # 0 = no long grass, 1 = long grass

## ABUNDANCE 100m ##

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_100m", "Arable_100m", "Grassland_100m", "Woodland_100m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_100m) # 0.28
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_100m) # -0.35
# remove grassland from analysis 

abund_100_mod <- lmer(log(total_abund) ~ garden_size + n_days + scale(Woodland_100m) + garden_ivy*scale(Urban_100m) +
                        garden_ivy*scale(Arable_100m) + scale(woodland_dist_m) + scale(I(Woodland_100m^2)) +
                        scale(grassland_dist_m) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis)
summary(abund_100_mod)

# backwards model selection
abund_100_mod_step <- lmerTest::step(abund_100_mod, reduce.random=FALSE, alpha.fixed=0.05)
abund_100_mod_final <- get_model(abund_100_mod_step)
summary(abund_100_mod_final)
# woodland and urban significant interaction

car::vif(abund_100_mod_final) 
r.squaredGLMM(abund_100_mod_final) # 4%
AIC(abund_100_mod_final) # 3220.11

testDispersion(abund_100_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_100_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_100_summary <- summary(abund_100_mod_final) #pulling out model averages
abund_100_summary<-as.data.frame(abund_100_summary$coefficients) #selecting full model coefficient averages (more conservative)
abund_100_summary$parameters <- row.names(abund_100_summary)
row.names(abund_100_summary) <- 1:nrow(abund_100_summary)
abund_100_summary$AIC <- AIC(abund_100_mod_final)
write.csv(abund_100_summary, file="Results/Abundance/Abundance_100m_ivy.csv", row.names=FALSE)

## ABUNDANCE 250m ##

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_250m", "Arable_250m", "Grassland_250m", "Woodland_250m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_250m) # 0.28
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_250m) # -0.35
# remove grassland from analysis 

abund_250_mod <- lmer(log(total_abund) ~ garden_size + n_days + scale(Woodland_250m) + garden_ivy*scale(Urban_250m) +
                        garden_ivy*scale(Arable_250m) + scale(woodland_dist_m) + 
                        scale(grassland_dist_m) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis)
summary(abund_250_mod)

# backwards model selection
abund_250_mod_step <- lmerTest::step(abund_250_mod, reduce.random=FALSE, alpha.fixed=0.05)
abund_250_mod_final <- get_model(abund_250_mod_step)
summary(abund_250_mod_final)
# woodland significant interaction

car::vif(abund_250_mod_final) 
r.squaredGLMM(abund_250_mod_final) # 4%
AIC(abund_250_mod_final) # 3227.87

testDispersion(abund_250_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_250_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_250_summary <- summary(abund_250_mod_final) #pulling out model averages
abund_250_summary<-as.data.frame(abund_250_summary$coefficients) #selecting full model coefficient averages (more conservative)
abund_250_summary$parameters <- row.names(abund_250_summary)
row.names(abund_250_summary) <- 1:nrow(abund_250_summary)
abund_250_summary$AIC <- AIC(abund_250_mod_final)
write.csv(abund_250_summary, file="Results/Abundance/Abundance_250m_ivy.csv", row.names=FALSE)

## ABUNDANCE 500m ##

# get response distribution
hist(gbs_analysis$total_abund)
gbs_analysis$total_abund_log <- log(gbs_analysis$total_abund)
hist(gbs_analysis$total_abund_log)
# use log transformed

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_500m", "Arable_500m", "Grassland_500m", "Woodland_500m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$total_abund, gbs_analysis$Grassland_500m) # 0.28
cor.test(gbs_analysis$total_abund, gbs_analysis$Urban_500m) # -0.35
# remove grassland from analysis 

abund_500_mod <- lmer(log(total_abund) ~ garden_size + n_days + scale(Woodland_500m) + garden_ivy*scale(Urban_500m) +
                        garden_ivy*scale(Arable_500m) + scale(woodland_dist_m) + scale(grassland_dist_m) + (1|M_YEAR), 
                        na.action = "na.fail", data=gbs_analysis)
summary(abund_500_mod)

# backwards model selection
abund_500_mod_step <- lmerTest::step(abund_500_mod, reduce.random=FALSE, alpha.fixed=0.05)
abund_500_mod_final <- get_model(abund_500_mod_step)
summary(abund_500_mod_final)
# urban (quadratic) and woodland significant interaction

car::vif(abund_500_mod_final) 
r.squaredGLMM(abund_500_mod_final) # 4%
AIC(abund_500_mod_final) # 3224.17 (v. similar to 250m)

testDispersion(abund_500_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = abund_500_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
abund_500_summary <- summary(abund_500_mod_final) #pulling out model averages
abund_500_summary<-as.data.frame(abund_500_summary$coefficients) #selecting full model coefficient averages (more conservative)
abund_500_summary$parameters <- row.names(abund_500_summary)
row.names(abund_500_summary) <- 1:nrow(abund_500_summary)
abund_500_summary$AIC <- AIC(abund_500_mod_final)
write.csv(abund_500_summary, file="Results/Abundance/Abundance_500m_ivy.csv", row.names=FALSE)


## RICHNESS 100m ##

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_100m", "Arable_100m", "Grassland_100m", "Woodland_100m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_100m) # 0.28
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_100m) # -0.35
# remove grassland from analysis 

rich_100_mod <- lmer(site_rel_SR ~ garden_size + n_days + scale(Woodland_100m) + garden_ivy*scale(Urban_100m) +
                        garden_ivy*scale(Arable_100m) + scale(woodland_dist_m) + 
                        scale(grassland_dist_m) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis)
summary(rich_100_mod)

# backwards model selection
rich_100_mod_step <- lmerTest::step(rich_100_mod, reduce.random=FALSE, alpha.fixed=0.05)
rich_100_mod_final <- get_model(rich_100_mod_step)
summary(rich_100_mod_final)
# urban significant interaction

car::vif(rich_100_mod_final) 
r.squaredGLMM(rich_100_mod_final) # 4%
AIC(rich_100_mod_final) # -3091.51

testDispersion(rich_100_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_100_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_100_summary <- summary(rich_100_mod_final) #pulling out model averages
rich_100_summary<-as.data.frame(rich_100_summary$coefficients) #selecting full model coefficient averages (more conservative)
rich_100_summary$parameters <- row.names(rich_100_summary)
row.names(rich_100_summary) <- 1:nrow(rich_100_summary)
rich_100_summary$AIC <- AIC(rich_100_mod_final)
write.csv(rich_100_summary, file="Results/Richness/Richness_100m_ivy.csv", row.names=FALSE)

ggpredict(rich_100_mod_final, terms="garden_ivy") %>% plot() # present slightly lower SR
ggpredict(rich_100_mod_final, terms="Urban_100m") %>% plot() # slightly lower SR at high urban levels


## RICHNESS 250m ##

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_250m", "Arable_250m", "Grassland_250m", "Woodland_250m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_250m) # 0.28
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_250m) # -0.35
# remove grassland from analysis 

rich_250_mod <- lmer(site_rel_SR ~ garden_size + n_days + scale(Woodland_250m) + garden_ivy*scale(Urban_250m) +
                       garden_ivy*scale(Arable_250m) + scale(woodland_dist_m) + 
                       scale(grassland_dist_m) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis)
summary(rich_250_mod)

# backwards model selection
rich_250_mod_step <- lmerTest::step(rich_250_mod, reduce.random=FALSE, alpha.fixed=0.05)
rich_250_mod_final <- get_model(rich_250_mod_step)
summary(rich_250_mod_final)
# woodland significant interaction

car::vif(rich_250_mod_final) 
r.squaredGLMM(rich_250_mod_final) # 4%
AIC(rich_250_mod_final) # -3084.69

testDispersion(rich_250_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_250_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_250_summary <- summary(rich_250_mod_final) #pulling out model averages
rich_250_summary<-as.data.frame(rich_250_summary$coefficients) #selecting full model coefficient averages (more conservative)
rich_250_summary$parameters <- row.names(rich_250_summary)
row.names(rich_250_summary) <- 1:nrow(rich_250_summary)
rich_250_summary$AIC <- AIC(rich_250_mod_final)
write.csv(rich_250_summary, file="Results/Richness/Richness_250m_ivy.csv", row.names=FALSE)


## RICHNESS 500m ##

# Check correlation between all variables
round(cor(gbs_analysis[,c("Urban_500m", "Arable_500m", "Grassland_500m", "Woodland_500m", 
                          "woodland_dist_m", "grassland_dist_m")]),3)
# urban and grassland -0.792
# threshold of 0.7 - need to remove with urban OR grassland

cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Grassland_500m) # 0.28
cor.test(gbs_analysis$site_rel_SR, gbs_analysis$Urban_500m) # -0.35
# remove grassland from analysis 

rich_500_mod <- lmer(site_rel_SR ~ garden_size + n_days + scale(Woodland_500m) + garden_ivy*scale(Urban_500m) +
                       garden_ivy*scale(Arable_500m) + scale(woodland_dist_m) + 
                       scale(grassland_dist_m) + (1|M_YEAR), na.action = "na.fail", data=gbs_analysis)
summary(rich_500_mod)

# backwards model selection
rich_500_mod_step <- lmerTest::step(rich_500_mod, reduce.random=FALSE, alpha.fixed=0.05)
rich_500_mod_final <- get_model(rich_500_mod_step)
summary(rich_500_mod_final)
# no significant interactions

car::vif(rich_500_mod_final) 
r.squaredGLMM(rich_500_mod_final) # 4%
AIC(rich_500_mod_final) # -3113.56

testDispersion(rich_500_mod_final) ## looks good
simulationOutput <- simulateResiduals(fittedModel = rich_500_mod_final, plot = F)
plot(simulationOutput) # no assumptions violated 

# save model output
rich_500_summary <- summary(rich_500_mod_final) #pulling out model averages
rich_500_summary<-as.data.frame(rich_500_summary$coefficients) #selecting full model coefficient averages (more conservative)
rich_500_summary$parameters <- row.names(rich_500_summary)
row.names(rich_500_summary) <- 1:nrow(rich_500_summary)
rich_500_summary$AIC <- AIC(rich_500_mod_final)
write.csv(rich_500_summary, file="Results/Richness/Richness_500m_ivy.csv", row.names=FALSE)
