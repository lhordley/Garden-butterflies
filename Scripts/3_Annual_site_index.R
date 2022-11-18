##########################
#### user: Lisbeth Hordley
#### date: July 2022
#### info: Calculate site abundance indices

#devtools::install_github("RetoSchmucki/rbms")
library(rbms)
library(data.table)
library(dplyr)
library(ggplot2)

# load cleaned data
gbs_data <- read.csv("Data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
# DAILY ANCHOR ARE CORRECT FILES

# Get data into right format for rbms
# create two files: gbs_count and gbs_visit
gbs_count <- gbs_data[,c("grid_reference","date","common_name","day","month","year","quantity")]
gbs_count <- gbs_count %>% transform(., SITE_ID=match(grid_reference, unique(grid_reference))) # assign a unique code to each site
site_match <- unique(gbs_count[,c("grid_reference", "SITE_ID")])
gbs_count$grid_reference <- NULL
colnames(gbs_count) <- c("DATE", "SPECIES", "DAY", "MONTH", "YEAR", "COUNT", "SITE_ID")
gbs_count <- gbs_count[,c(7,1:6)]

gbs_visit <- unique(gbs_count[,c("SITE_ID", "DATE")])

#################################
# 1. From counts to flight curve
#################################

# initialise a time-series with day-week-month-year information
ts_date <- rbms::ts_dwmy_table(InitYear = 2016, LastYear = 2021, WeekDay1 = 'monday')
# each day each year from 2016 to 2021

# add monitoring season to the new time-series, providing startmonth and endmonth arguments
# note: can also refine further by adding in startday and endday
# also define resolution of time series (weekly or daily, timeunit='w' or 'd')
# anchor argument adds zeros before and after monitoring season

ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 1, EndMonth = 12, StartDay = 1, EndDay = 31, 
                                   CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 2, AnchorLag = -2,
                                   TimeUnit = 'd') # ts date and ts season same length

# Now add site visits to the time series
ts_season_visit <- rbms::ts_monit_site(ts_season, gbs_count)

# create loop to calculate and plot flight curve, and calculate sindex and plot annual indices for each species, each year
species <- unique(gbs_data$common_name) # 31 species
sindex_final <- NULL
pheno_final <- NULL
annual_indices <- NULL
#options(warn=2)
for(i in species) {
print(i)
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, gbs_count, sp = i)
  ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 100, MinVisit = 4, MinOccur = 2, 
                                        MinNbrSite = 5, MaxTrial = 3, GamFamily = 'poisson', SpeedGam = FALSE, 
                                        CompltSeason = TRUE, SelectYear = NULL, TimeUnit = 'd')
  
  # extract pheno object: contains the shape of annual flight curves, standardised to sum 1
  pheno <- ts_flight_curve$pheno
  pheno$sp <- i
  pheno_final <- rbind(pheno, pheno_final, fill=TRUE)
  
  # extract site index for each site, year and species
  # impute_count() function uses the count data generated from ts_season_count() function and the flight curves
  # it looks for the phenology available to estimate and input missing values
  # imputation are made on a daily basis

  impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='d')

  # this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year,
  # and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season

  # From the imputed count, the site index can be calculated for each site, or with a filter that will only keep sites that have
  # been monitored at least a certain proportion of the flight curve

  sindex_temp <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.10) # only keep sites that have been monitored for 10% of flight curve
  sindex_final <- rbind(sindex_temp, sindex_final) # this is the annual abundance indices for each site

  # extract annual collated index and plot to compare against UKBMS trends
  co_index <- collated_index(data = sindex_temp, s_sp = i, sindex_value = "SINDEX", glm_weights = TRUE, rm_zero = TRUE)
  co_index <- co_index$col_index
  # NSITE_OBS is number of sites that the species has been recorded at - assume all these are used to calculate collated index
  # NSITE is the number of sites (total) for that year
  # transform index to a log(10) scale
  co_index_b <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
  co_index_logInd <- co_index_b[BOOTi == 0, .(M_YEAR, COL_INDEX)][, log(COL_INDEX)/log(10), by = M_YEAR][, mean_logInd := mean(V1)]
  data.table::setnames(co_index_logInd, "V1", "logInd"); setkey(co_index_logInd, M_YEAR); setkey(co_index_b, M_YEAR)
  co_index_b <- merge(co_index_b, co_index_logInd, all.x = TRUE)
  b1 <- data.table(M_YEAR = co_index_b$M_YEAR, LCI = 2 + co_index_b$logInd - co_index_b$mean_logInd)
  b1 <- merge(b1, co_index, by="M_YEAR")
  b1$sp <- i
  annual_indices <- rbind(b1, annual_indices) # this is the annual collated index across all sites

}
# Can't estimate flight curves for Small Blue - left with 30 species
options(op)

# This produces 3 files:
# 1. Pheno_final = each species' flight curves based on GBS data with anchors on 1st + 2nd Jan and 30th + 31st Dec
# 2. Sindex_final = annual abundance index for each site across all 30 species based on raw counts and imputed values from flight curve
# 3. Annual_indices = annual collated index across all sites and species

length(unique(sindex_final$SPECIES)) # 30 species
length(unique(annual_indices$sp)) # 30 species
length(unique(pheno_final$sp)) # 31 species

sindex_final <- merge(sindex_final, site_match, by="SITE_ID", all.x=TRUE) # add grid references back in
# save each file
write.csv(sindex_final, file="Data/Site_index_GBS_daily.csv", row.names=FALSE)
write.csv(annual_indices, file="Data/Annual_collated_index_daily.csv", row.names=FALSE)
saveRDS(pheno_final, file="Data/Flight_curves_daily.rds")


################ Put this into new script! 
annual_indices <- read.csv("Data/Annual_collated_index_daily.csv", header=TRUE)
# read in UKBMS collated indices and check similarity
ukbms_indices <- read.csv("Data/GB_GAI_collated_indices_1976-2021.csv", header=TRUE)

ukbms_indices$data <- "UKBMS"
annual_indices$data <- "GBS"
species <- "Brimstone"

ukbms_indices <- ukbms_indices[ukbms_indices$YEAR>=2016,]

ukbms_indices <- ukbms_indices[,c("COMMON_NAME","YEAR","TRMOBS","data")]
ukbms_indices2 <- ukbms_indices %>% group_by(COMMON_NAME) %>% summarise(TRMOBS=TRMOBS-mean(TRMOBS)+2)
ukbms_indices2$YEAR <- ukbms_indices$YEAR
ukbms_indices2$data <- "UKBMS"
colnames(ukbms_indices2) <- c("sp", "LCI", "M_YEAR", "data")
ukbms_indices <- ukbms_indices[which(ukbms_indices$sp %in% species), ]
annual_indices <- annual_indices[,c("sp","M_YEAR","LCI","data")]

all_indices <- rbind(ukbms_indices2, annual_indices)


for(i in species){print(i)
  
  temp_plot <- ggplot(all_indices[all_indices$sp==i,], aes(x=M_YEAR, y=LCI, colour=data))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=2, linetype='dotted', col = 'black')+
    ggtitle(i)+
    labs(x="Year", y=expression('log '['(10)']*' Collated Index'))+
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  #ggsave(temp_plot, file=paste0("../Graphs/Species collated indices/Collated_index_comparison_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}



## Put flight curves and comparison collated index plot in same pdf 

library(ggplot2)
library(gridExtra)

annual_indices <- read.csv("../Data/Annual_collated_index_daily.csv", header=TRUE)
# read in UKBMS collated indices and check similarity
ukbms_indices <- read.csv("../Data/GB_GAI_collated_indices_1976-2021.csv", header=TRUE)

ukbms_indices$data <- "UKBMS"
annual_indices$data <- "GBS"
species <- unique(annual_indices$sp)

ukbms_indices <- ukbms_indices[ukbms_indices$YEAR>=2016,]

ukbms_indices <- ukbms_indices[,c("COMMON_NAME","YEAR","TRMOBS","data")]
colnames(ukbms_indices) <- c("sp", "M_YEAR", "LCI", "data")
ukbms_indices <- ukbms_indices[which(ukbms_indices$sp %in% species), ]
annual_indices <- annual_indices[,c("sp","M_YEAR","LCI","data")]

all_indices <- rbind(ukbms_indices, annual_indices)
all_indices <- all_indices[!all_indices$sp=="Small Blue",]


pheno <- readRDS("../Data/Flight_curves_daily.rds")
pheno <- pheno[!pheno$SPECIES=="Small Blue",] # can't estimate flight curve
sp<-unique(pheno$SPECIES)

# also remove Dark Green Fritillary - can only estimate filght curve in 2020
# White-letter Hairstreak
# Brown Hairstreak

pdf('../Graphs/Flight_curves_collated_indices.pdf', height = 8, width = 6, onefile=TRUE)
for (i in sp) {
  print(i)
  flight_curve <- ggplot(pheno[pheno$SPECIES==i,], aes(x=trimDAYNO, y=NM, colour=M_YEAR))+
    geom_line()+
    labs(x="Monitoring Day", y="Relative Abundance", colour="")+
    ggtitle(i)+
    scale_x_continuous(breaks=seq(0,366, by=50)) +
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  
  collated_index <- ggplot(all_indices[all_indices$sp==i,], aes(x=M_YEAR, y=LCI, colour=data))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=2, linetype='dotted', col = 'black')+
    labs(x="Year", y=expression('log '['(10)']*' Collated Index'))+
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  
  grid.arrange(flight_curve, collated_index)
}
dev.off()





