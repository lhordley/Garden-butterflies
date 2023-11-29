##########################
#### user: Lisbeth Hordley
#### date: September 2022
#### info: Compare UKBMS and GBS collated indices
options(scipen = 100)

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

## read in collated indices 
gbs_indices <- read.csv("Data/Abundance data/Annual_collated_index_daily.csv", header=TRUE)
ukbms_indices <- read.csv("Data/Abundance data/GB_GAI_collated_indices_1976-2021.csv", header=TRUE)

ukbms_indices$data <- "UKBMS"
gbs_indices$data <- "GBS"

species <- unique(gbs_indices$sp)
ukbms_indices <- ukbms_indices[ukbms_indices$COMMON_NAME %in% species,]
ukbms_indices <- ukbms_indices[ukbms_indices$YEAR>=2016,]

ukbms_indices <- ukbms_indices[,c("COMMON_NAME","YEAR","TRMOBS","data")]
ukbms_indices <- ukbms_indices %>% group_by(COMMON_NAME) %>% mutate(TRMOBS=2+TRMOBS-mean(TRMOBS)) %>% ungroup()
colnames(ukbms_indices) <- c("sp", "M_YEAR", "LCI", "data")
gbs_indices <- gbs_indices[,c("sp","M_YEAR","LCI","data")]

all_indices <- rbind(ukbms_indices, gbs_indices)

# create graph comparing the two schemes for each species
species <- unique(gbs_indices$sp)
for(i in species){print(i)
  
  temp_plot <- ggplot(all_indices[all_indices$sp==i,], aes(x=M_YEAR, y=LCI, colour=data))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=2, linetype='dotted', col = 'black')+
    ggtitle(i)+
    labs(x="Year", y=expression('log '['(10)']*' Collated Index'))+
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  ggsave(temp_plot, file=paste0("Graphs/Species collated indices/Collated_index_comparison_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}

###### Quantitatively compare UKBMS and GBS collated indices ######

# Compare change in annual collated indices between years between UKBMS and GBS
annual_gr <- all_indices %>% group_by(sp, data) %>% mutate(AGR=LCI-lag(LCI))
annual_gr <- na.omit(annual_gr)
annual_gr$LCI <- NULL
# annual_gr <- annual_gr[!annual_gr$sp=="Painted Lady",]

annual_gr2<- spread(annual_gr, key = data, value = AGR) # change data from wide to long

annual_gr3 <- annual_gr2 %>% mutate(differences = UKBMS - GBS)

# Plot graphs 
annual_gr2$M_YEAR <- as.factor(annual_gr2$M_YEAR)
AGR_plot <- ggplot(annual_gr2, aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS change in annual \ncollated index", y="GBS change in annual \ncollated index")+
  scale_y_continuous(limits=c(-1.8,1.4))+
  scale_x_continuous(limits=c(-1.8,1.4))+
  geom_abline(intercept = 0, slope = 1, size = 1, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
AGR_plot
# save file
ggsave(AGR_plot, file="Graphs/UKBMS_GBS_AGR_comparison_23spp.png")

annual_gr3 <- annual_gr3[!annual_gr3$sp=="Painted Lady",]

AGR_plot2 <- ggplot(annual_gr3, aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS change in annual \ncollated index", y="GBS change in annual \ncollated index")+
  scale_y_continuous(limits=c(-0.7,0.7))+
  scale_x_continuous(limits=c(-0.7,0.7))+
  geom_abline(intercept = 0, slope = 1, size = 1, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
AGR_plot2
# save file
ggsave(AGR_plot2, file="Graphs/UKBMS_GBS_AGR_comparison_22spp.png")

library(ggpubr)
agr_plots <- ggarrange(AGR_plot, AGR_plot2, labels = c("(a)", "(b)"),
                       ncol = 2, nrow = 1)
agr_plots
ggsave(agr_plots, file="Graphs/AGR_comparison_final.png", height=4, width=9)

### Paired t-test across all years ### 
res <- t.test(annual_gr3$UKBMS, annual_gr3$GBS, paired = TRUE)
res # non-significant for all species (still non-significant when PL removed)

ggboxplot(annual_gr, x = "data", y = "AGR", 
          color = "data", palette = c("#00AFBB", "#E7B800"),
          order = c("GBS", "UKBMS")) # no difference overall


### Paired t-test for each consecutive year comparison ###
t_test_final <- NULL
years <- unique(annual_gr3$M_YEAR)

for(i in years){print(i)
  annual_gr4 <- annual_gr3[annual_gr3$M_YEAR==i,]
   
  t_test <- t.test(annual_gr4$UKBMS, annual_gr4$GBS, paired = TRUE) 
  t_test_temp <- data.frame(year=i,mean_diff=t_test$estimate,se=t_test$stderr,
                            t_value=t_test$statistic, p_value=t_test$p.value)
  t_test_final <- rbind(t_test_temp, t_test_final)
}
# 2017 still significant, 2021 just non-significant

# Plot graphs for each consecutive year comparison

AGR_plot1 <- ggplot(annual_gr2[annual_gr2$Year=="2016-2017",], aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS change in \nannual collated index", y="GBS change in \nannual collated index")+
  scale_y_continuous(limits=c(-0.3,0.6))+
  scale_x_continuous(limits=c(-0.3,0.6))+
  geom_abline(intercept = 0, slope = 1, size = 0.5, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))+
  facet_grid(~Year)
AGR_plot1

AGR_plot2 <- ggplot(annual_gr2[annual_gr2$Year=="2017-2018",], aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS change in \nannual collated index", y="GBS change in \nannual collated index")+
  scale_y_continuous(limits=c(-0.7,0.6))+
  scale_x_continuous(limits=c(-0.7,0.6))+
  geom_abline(intercept = 0, slope = 1, size = 0.5, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))+
  facet_grid(~Year)
AGR_plot2

AGR_plot3 <- ggplot(annual_gr2[annual_gr2$Year=="2018-2019",], aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS change in \nannual collated index", y="GBS change in \nannual collated index")+
  scale_y_continuous(limits=c(-0.5,1.4))+
  scale_x_continuous(limits=c(-0.5,1.4))+
  geom_abline(intercept = 0, slope = 1, size = 0.5, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))+
  facet_grid(~Year)
AGR_plot3

AGR_plot4 <- ggplot(annual_gr2[annual_gr2$Year=="2019-2020",], aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS change in \nannual collated index", y="GBS change in \nannual collated index")+
  scale_y_continuous(limits=c(-1.8,0.5))+
  scale_x_continuous(limits=c(-1.8,0.5))+
  geom_abline(intercept = 0, slope = 1, size = 0.5, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))+
  facet_grid(~Year)
AGR_plot4

AGR_plot5 <- ggplot(annual_gr2[annual_gr2$Year=="2020-2021",], aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS change in \nannual collated index", y="GBS change in \nannual collated index")+
  scale_y_continuous(limits=c(-0.4,0.7))+
  scale_x_continuous(limits=c(-0.4,0.7))+
  geom_abline(intercept = 0, slope = 1, size = 0.5, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(text = element_text(size=16), legend.title=element_text(size=16), 
        legend.text=element_text(size=16))+
  facet_grid(~Year)
AGR_plot5

AGR_plots <- ggarrange(AGR_plot1, AGR_plot2, AGR_plot3, AGR_plot4, AGR_plot5, labels=c("(a)", "(b)", "(c)", "(d)", "(e)"),
                                                                                       nrow=1, ncol=5)
AGR_plots
ggsave(AGR_plots, file="Graphs/UKBMS_GBS_AGR_comparison_yearly.png", height=4, width=17)

