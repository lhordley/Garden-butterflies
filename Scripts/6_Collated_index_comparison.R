##########################
#### user: Lisbeth Hordley
#### date: September 2022
#### info: Compare UKBMS and GBS collated indices

library(tidyr)

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

# Quantitatively compare UKBMS and GBS collated indices
all_indices2<- spread(all_indices, key = data, value = LCI)
# are the data normally distributed?
# Shapiro-Wilk normality test for UKBMS
shapiro.test(all_indices2$UKBMS) # => p < 0.001
# Shapiro-Wilk normality test for GBS
shapiro.test(all_indices2$GBS) # => p < 0.001
# Not normally distributed

ggqqplot(all_indices2$GBS, ylab = "GBS")

cor.test(all_indices2$UKBMS, all_indices2$GBS, 
         method = "spearman")
# r=0.87, p<0.001

# Plot data
collated_index_plot <- ggplot(all_indices2, aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS Collated Index", y="GBS Collated Index")+
  scale_y_continuous(limits=c(1.3,3.2))+
  scale_x_continuous(limits=c(1.3,3.2))+
  geom_abline(intercept = 0, slope = 1, size = 1, linetype="dashed", colour="grey") +
  annotate("text", x=1.5, y=3.1, label="r=0.87, p<0.001")+
  theme_classic()
collated_index_plot
# save file
ggsave(collated_index_plot, file="Graphs/UKBMS_GBS_collated_index_comparison.png")

# Also compare Annual Growth Rates (AGRs) between UKBMS and GBS
# This is just the change in collated index from the previous year
# Choose this method because comparing collated indices means comparing the exact values, but because the methodology and
# survey time differ between UKBMS and GBS, this isn't a fair comparison
# Whereas comparing AGRs is a fairer comparison as we're comparing the rate of change rather than an exact yearly value

annual_gr <- all_indices %>% group_by(sp, data) %>% mutate(AGR=LCI-lag(LCI))
annual_gr <- na.omit(annual_gr)
annual_gr$LCI <- NULL

# Quantitatively compare UKBMS and GBS AGRs
annual_gr2<- spread(annual_gr, key = data, value = AGR)
# are the data normally distributed?
# Shapiro-Wilk normality test for UKBMS
shapiro.test(annual_gr2$UKBMS) # => p < 0.001
# Shapiro-Wilk normality test for GBS
shapiro.test(annual_gr2$GBS) # => p < 0.001
# Not normally distributed

ggqqplot(annual_gr2$UKBMS, ylab = "UKBMS")

cor.test(annual_gr2$UKBMS, annual_gr2$GBS, 
         method = "spearman")
# r=0.85, p<0.001

# Plot data
AGR_plot <- ggplot(annual_gr2, aes(x=UKBMS, y=GBS))+
  geom_point(size=3, alpha=0.2)+
  labs(x="UKBMS Annual Growth Rate", y="GBS Annual Growth Rate")+
  scale_y_continuous(limits=c(-1.8,1.4))+
  scale_x_continuous(limits=c(-1.8,1.4))+
  geom_abline(intercept = 0, slope = 1, size = 1, linetype="dashed", colour="grey") +
  annotate("text", x=-1.4, y=1.2, label="r=0.85, p<0.001")+
  theme_classic()
AGR_plot
# save file
ggsave(AGR_plot, file="Graphs/UKBMS_GBS_AGR_comparison.png")
