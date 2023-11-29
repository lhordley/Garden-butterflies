##########################
#### user: Lisbeth Hordley
#### date: July 2022
#### info: Calculate site abundance indices

#devtools::install_github("RetoSchmucki/rbms")
library(rbms)
library(data.table)
library(dplyr)
library(ggplot2)

####################### Abundance calculation for all sites and species ###################################

# load cleaned data
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE) # main dataset
# gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered_top_rec.csv", header=TRUE) # top 50% recorded sites

# gbs_sites <- read.csv("Data/GBS_sites_GEnS_CTM.csv", header=TRUE)
# 
# # filter to only include sites within one bioclimatic zone (cool, temperate and moist)
# gbs_data <- gbs_data[gbs_data$grid_reference %in% gbs_sites$grid_reference,]
# length(unique(gbs_data$grid_reference)) # 785

# Get data into right format for rbms
# create two files: gbs_count and gbs_visit
gbs_count <- gbs_data[,c("grid_reference","date","common_name","day","month","year","quantity")]
gbs_count <- gbs_count %>% transform(., SITE_ID=match(grid_reference, unique(grid_reference))) # assign a unique code to each site
site_match <- unique(gbs_count[,c("grid_reference", "SITE_ID")])
gbs_count$grid_reference <- NULL
colnames(gbs_count) <- c("DATE", "SPECIES", "DAY", "MONTH", "YEAR", "COUNT", "SITE_ID")
gbs_count <- gbs_count[,c(7,1:6)]

gbs_count <- unique(gbs_count)

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
ts_season_visit <- rbms::ts_monit_site(ts_season, gbs_visit)

# create loop to calculate and plot flight curve, and calculate sindex and plot annual indices for each species, each year
species <- unique(gbs_data$common_name) # 31 species
pheno_final <- NULL
#options(warn=2)
start_time <- Sys.time() 
for(i in species) {
print(i)
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, gbs_count, sp = i)
  ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 300, MinVisit = 4, MinOccur = 2, 
                                        MinNbrSite = 5, MaxTrial = 3, GamFamily = 'poisson', SpeedGam = FALSE, 
                                        CompltSeason = TRUE, SelectYear = NULL, TimeUnit = 'd')
  
  # extract pheno object: contains the shape of annual flight curves, standardised to sum 1
  pheno <- ts_flight_curve$pheno
  pheno$sp <- i
  pheno_final <- rbind(pheno, pheno_final, fill=TRUE)
}
# Can't estimate flight curves for Small Blue - left with 30 species

# This produces pheno_final = each species' flight curves based on GBS data with anchors on 1st + 2nd Jan and 30th + 31st Dec
length(unique(pheno_final$sp)) # 31 species (Small blue is just NAs or zeros)
end_time <- Sys.time() # takes ~ 7 hours to calculate flight curves

# save flight curves for all 31 species 
saveRDS(pheno_final, file="Data/Abundance data/Flight_curves_GBS_daily_anchor.rds")

###
pheno_final <- readRDS("Data/Abundance data/Flight_curves_GBS_daily_anchor.rds")
pheno_final2 <- readRDS("Data/Abundance data/Flight_curves_GBS_daily_anchor_GEnS.rds")

### Some comparisons of flight curves with and without climate zone 
library(gridExtra)
pdf('Graphs/Species flight curves/Flight_curve_comparison.pdf',width = 14)
for (i in unique(pheno_final$SPECIES)) {
  print(i)
  flight_curve1 <- ggplot(na.omit(pheno_final[pheno_final$SPECIES==i,]), aes(x=trimDAYNO, y=NM, colour=M_YEAR))+
    geom_line()+
    labs(x="Monitoring Day", y="Relative Abundance", colour="")+
    ggtitle(i, ": No environmental zone")+
    #scale_x_continuous(breaks=seq(0,366, by=50)) +
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  
  flight_curve2 <- ggplot(na.omit(pheno_final2[pheno_final2$SPECIES==i,]), aes(x=trimDAYNO, y=NM, colour=M_YEAR))+
    geom_line()+
    labs(x="Monitoring Day", y="Relative Abundance", colour="")+
    ggtitle(i, ": Environmental zone")+
    #scale_x_continuous(breaks=seq(0,366, by=50)) +
    theme_classic()+ 
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  grid.arrange(flight_curve1, flight_curve2, nrow = 1, ncol=2)
  # plot(flight_curve1)
  # plot(flight_curve2)
  #ggsave(flight_curve, file=paste0("Graphs/Species flight curves/Flight_curve_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}
dev.off()

pheno_1 <- pheno[,c("SPECIES", "SITE_ID", "YEAR", "MONTH", "DAY", "NM")]
pheno_2 <- pheno2[,c("SPECIES", "SITE_ID", "YEAR", "MONTH", "DAY", "NM")]
colnames(pheno_1)[6] <- "NM_nozone"
colnames(pheno_2)[6] <- "NM_climatezone"

pheno_final <- merge(pheno_1, pheno_2, by=c("SPECIES", "SITE_ID", "YEAR", "MONTH", "DAY"))
pheno_final$NM_diff <- pheno_final$NM_nozone - pheno_final$NM_climatezone

hist(pheno_final$NM_diff) # most around zero
pheno_final$YEAR <- as.factor(pheno_final$YEAR)
# plot for each species
pdf('Graphs/Species flight curves/Flight_curve_difference.pdf',width = 14)
for (i in unique(pheno_final$SPECIES)) {
  print(i)
  flight_curve1 <- ggplot(na.omit(pheno_final[pheno_final$SPECIES==i,]), aes(x=NM_nozone, y=NM_climatezone, colour=YEAR))+
    geom_line()+
    labs(x="Relative Abundnace [original]", y="Relative Abundance [climate zone]", colour="")+
    ggtitle(i)+
    #scale_x_continuous(breaks=seq(0,366, by=50)) +
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  plot(flight_curve1)
  Sys.sleep(2)
}
dev.off()



###################################################
# 2. From flight curve to site and collated indices
###################################################

pheno <- readRDS("Data/Abundance data/Flight_curves_GBS_daily_anchor.rds")
pheno2 <- readRDS("Data/Abundance data/Flight_curves_GBS_daily_anchor_GEnS.rds")
# Make sure to create ts_season_visit before running this section

## Calculate site indices

species <- unique(gbs_count$SPECIES) # 31 species
sindex_final <- NULL
annual_indices <- NULL
#options(warn=2)
for(i in species) {
  print(i)
  
# extract site index for each site, year and species
# impute_count() function uses the count data generated from ts_season_count() function and the flight curves
# it looks for the phenology available to estimate and input missing values
# imputation are made on a daily basis
ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, gbs_count, sp = i)
tryCatch({
  impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='d')
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year,
# and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season

# From the imputed count, the site index can be calculated for each site, or with a filter that will only keep sites that have
# been monitored at least a certain proportion of the flight curve

sindex_temp <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.10) # only keep sites that have been monitored for 10% of flight curve
sindex_temp$n_observed <- length(na.omit(impt_counts$COUNT[impt_counts$COUNT>0]))
sindex_final <- rbind(sindex_temp, sindex_final) # this is the annual abundance indices for each site
# Species with errors: Dark Green Fritillary, Wall, White-letter Hairstreak, Clouded Yellow, Purple Hairstreak
# Small Blue, Green Hairstreak and Essex Skipper
# These species get removed below anyway because they do not have a flight curve for each year

}
sindex_final <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)


## Calculate annual abundance indices

species <- unique(sindex_final$SPECIES) # 31 species
annual_indices <- NULL
#options(warn=2)
for(i in species) {
  print(i)
  
sindex_temp <- sindex_final[sindex_final$SPECIES==i,]
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

# This produces 2 files:
# 1. Sindex_final = annual abundance index for each site across all 30 species based on raw counts and imputed values from flight curve
# 2. Annual_indices = annual collated index across all sites and species

sindex_dup <- sindex_final[duplicated(sindex_final) | duplicated(sindex_final, fromLast=TRUE), ]
# not sure why there are duplicated entries...
sindex_final <- unique(sindex_final)

length(unique(sindex_final$SPECIES)) # 30 species
length(unique(annual_indices$sp)) # 30 species

write.csv(annual_indices, file="Data/Abundance data/Annual_collated_index_daily_GEnS.csv", row.names=FALSE)

annual_indices <- read.csv("Data/Abundance data/Annual_collated_index_daily_GEnS.csv", header=TRUE)

# Which species don't have a flight curve estimated for each year? 
# Can still estimate an abundance index for these species, but it uses a flight curve from a different year
# Plot the flight curves to see which species to exclude:

for (i in unique(pheno_final$SPECIES)) {
  print(i)
  flight_curve <- ggplot(na.omit(pheno_final[pheno_final$SPECIES==i,]), aes(x=trimDAYNO, y=NM, colour=M_YEAR))+
    geom_line()+
    labs(x="Monitoring Day", y="Relative Abundance", colour="")+
    ggtitle(i)+
    #scale_x_continuous(breaks=seq(0,366, by=50)) +
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  
  ggsave(flight_curve, file=paste0("Graphs/Species flight curves/Flight_curve_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}

# Also create a way of doing this using pheno_final? Then filter based on that
# filter by length of M_YEAR - species need to have at least 6 years (2016-2021)
flight_curve_years <- pheno %>% na.omit() %>% group_by(SPECIES) %>% summarise(n_years=n_distinct(M_YEAR))
# 7 species do not have enough data to estimate a flight curve for each year
# Remove these species from sindex_final and annual_indices
sp_exclude <- subset(flight_curve_years, n_years<6, select=SPECIES)
sindex_final <- sindex_final[!sindex_final$SPECIES %in% sp_exclude$SPECIES,]
length(unique(sindex_final$SPECIES)) # 23 species
# add grid references back in
sindex_final <- merge(sindex_final, site_match, by="SITE_ID", all.x=TRUE) 
# remove zeros where a species was not observed at a site
sindex_final <- sindex_final[!sindex_final$SINDEX==0,]
# save site index file = abundance values for 23 species at 823 sites from 2016-2021
write.csv(sindex_final, file="Data/Abundance data/Site_index_GBS_daily.csv", row.names=FALSE)

# Exclude species from annual indices
annual_indices <- annual_indices[!annual_indices$sp %in% sp_exclude$SPECIES,]
length(unique(annual_indices$sp)) # 23
write.csv(annual_indices, file="Data/Abundance data/Annual_collated_index_daily.csv", row.names=FALSE)

## Compare site indices for each species with and without climate zone
sindex_final <- read.csv("Data/Abundance data/Site_index_GBS_daily.csv", header=TRUE)
sindex_final2 <- read.csv("Data/Abundance data/Site_index_GBS_daily_GEnS.csv", header=TRUE)

length(unique(sindex_final$SITE_ID))
sindex_final2 <- unique(sindex_final2)

sindex_final$TOTAL_NM <- NULL
sindex_final2$TOTAL_NM <- NULL
sindex_final$SITE_ID <- NULL
sindex_final2$SITE_ID <- NULL

colnames(sindex_final)[3] <- "SINDEX_nozone"
colnames(sindex_final2)[3] <- "SINDEX_climatezone"

sindex_final3 <- merge(sindex_final, sindex_final2, by=c("SPECIES", "M_YEAR", "grid_reference"))
sindex_final3$M_YEAR <- as.factor(sindex_final3$M_YEAR)
pdf('Graphs/Site_index_comparison.pdf',width = 14)
for (i in unique(sindex_final3$SPECIES)) {
  print(i)
comparison <- ggplot(na.omit(sindex_final3[sindex_final3$SPECIES==i,]), aes(x=SINDEX_nozone, y=SINDEX_climatezone, colour=M_YEAR))+
  geom_point()+
  labs(x="Site Index [original]", y="Site Index [climate zone]", colour="")+
  ggtitle(i)+
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype="dashed", colour="grey") +
  theme_classic()+
  theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
plot(comparison)
Sys.sleep(2)
}
dev.off()





############################# Abundance calculation for ivy (autumn) #####################################

# load cleaned data
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered_ivy.csv", header=TRUE)
# DAILY ANCHOR ARE CORRECT FILES

# gbs_sites <- read.csv("Data/GBS_sites_GEnS_CTM.csv", header=TRUE) # 785 sites
# # filter to only include sites within one bioclimatic zone (cool, temperate and moist)
# gbs_data <- gbs_data[gbs_data$grid_reference %in% gbs_sites$grid_reference,]
# length(unique(gbs_data$grid_reference)) # 592

# Get data into right format for rbms
# create two files: gbs_count and gbs_visit
gbs_count <- gbs_data[,c("grid_reference","date","common_name","day","month","year","quantity")]
gbs_count <- gbs_count %>% transform(., SITE_ID=match(grid_reference, unique(grid_reference))) # assign a unique code to each site
site_match <- unique(gbs_count[,c("grid_reference", "SITE_ID")])
gbs_count$grid_reference <- NULL
colnames(gbs_count) <- c("DATE", "SPECIES", "DAY", "MONTH", "YEAR", "COUNT", "SITE_ID")
gbs_count <- gbs_count[,c(7,1:6)]

gbs_visit <- unique(gbs_count[,c("SITE_ID", "DATE")])

# initialise a time-series with day-week-month-year information
ts_date <- rbms::ts_dwmy_table(InitYear = 2016, LastYear = 2021, WeekDay1 = 'monday')
# each day each year from 2016 to 2021

# add monitoring season to the new time-series, providing startmonth and endmonth arguments
# note: can also refine further by adding in startday and endday
# also define resolution of time series (weekly or daily, timeunit='w' or 'd')
# anchor argument adds zeros before and after monitoring season

# Try starting Holly Blue sindex from 19th June to capture second flight period - day 170
ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 9, EndMonth = 11, StartDay = 1, EndDay = 30, 
                                   CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 2, AnchorLag = -2,
                                   TimeUnit = 'd') # ts date and ts season same length

# Now add in site and count columns with count values as zeros where M_SEASON == 1 (i.e for september-november each year) and NA 
# for all other days
ts_season_visit <- rbms::ts_monit_site(ts_season, gbs_count)

pheno <- readRDS("Data/Abundance data/Flight_curves_GBS_daily_anchor.rds")

species <- unique(pheno$SPECIES) # 31 species
species2 <- unique(gbs_data$common_name)
sp_final <- intersect(species, species2)
sindex_final <- NULL
#options(warn=2)
# 1:8,10:14,16:17,20:23,25:26,29,31
for(i in sp_final) {print(i) # 23 species with records in Sept-Nov
  
  if(i=="Small Blue"){next}
  
  # extract site index for each site, year and species
  # impute_count() function uses the count data generated from ts_season_count() function and the flight curves
  # it looks for the phenology available to estimate and input missing values
  # imputation are made on a daily basis
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, gbs_count, sp = i) # adds in observed counts for species i
  tryCatch({
    impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='d')
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # this catches errors where there is no flight curve during sept-nov
  # no imputed counts when m_season = 0 (i.e. outside of sept-nov)
  # this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year,
  # and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season
  
  # From the imputed count, the site index can be calculated for each site, or with a filter that will only keep sites that have
  # been monitored at least a certain proportion of the flight curve
  
  sindex_temp <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.1) # not using minFC here - as that removes sites that look like
  # they've only monitored a small proportion of flight curve, but they have monitored more, it's just because we're restricting to
  # months sept-nov
  # to get around this - limit sites to those used in main GBS analysis which were restricted to minFC of 10% for whole year
  sindex_final <- rbind(sindex_temp, sindex_final) # this is the annual abundance indices for each site
  
}

length(unique(sindex_final$SPECIES)) # 15 species

# filter by length of M_YEAR - species need to have at least 6 years (2016-2021)
flight_curve_years <- pheno %>% na.omit() %>% group_by(SPECIES) %>% summarise(n_years=n_distinct(M_YEAR))
# 7 species do not have enough data to estimate a flight curve for each year
# Remove these species from sindex_final and annual_indices
sp_exclude <- subset(flight_curve_years, n_years<6, select=SPECIES)
sindex_final <- sindex_final[!sindex_final$SPECIES %in% sp_exclude$SPECIES,]
length(unique(sindex_final$SPECIES)) # still 15 species
# add grid references back in
sindex_final <- merge(sindex_final, site_match, by="SITE_ID", all.x=TRUE) 
# remove zeros where a species was not observed at a site
sindex_final <- sindex_final[!sindex_final$SINDEX==0,]
length(unique(sindex_final$grid_reference)) # 746
# save site index file = abundance values for 23 species at 746 sites from 2016-2021
write.csv(sindex_final, file="Data/Abundance data/Site_index_GBS_daily_ivy.csv", row.names=FALSE)



############################# Abundance calculation for ivy (summer) #####################################
# Holly Blue (2nd gen) only

# load cleaned data
gbs_data <- read.csv("Data/Raw GBS data/GBS_2016_2021_cleaned_filtered.csv", header=TRUE)
# DAILY ANCHOR ARE CORRECT FILES

length(unique(gbs_data$grid_reference))

# Get data into right format for rbms
# create two files: gbs_count and gbs_visit
gbs_count <- gbs_data[,c("grid_reference","date","common_name","day","month","year","quantity")]
gbs_count <- gbs_count %>% transform(., SITE_ID=match(grid_reference, unique(grid_reference))) # assign a unique code to each site
site_match <- unique(gbs_count[,c("grid_reference", "SITE_ID")])
gbs_count$grid_reference <- NULL
colnames(gbs_count) <- c("DATE", "SPECIES", "DAY", "MONTH", "YEAR", "COUNT", "SITE_ID")
gbs_count <- gbs_count[,c(7,1:6)]

gbs_visit <- unique(gbs_count[,c("SITE_ID", "DATE")])

# initialise a time-series with day-week-month-year information
ts_date <- rbms::ts_dwmy_table(InitYear = 2016, LastYear = 2021, WeekDay1 = 'monday')
# each day each year from 2016 to 2021

# add monitoring season to the new time-series, providing startmonth and endmonth arguments
# note: can also refine further by adding in startday and endday
# also define resolution of time series (weekly or daily, timeunit='w' or 'd')
# anchor argument adds zeros before and after monitoring season

# Try starting Holly Blue sindex from 29th June to capture second flight period - day 180 (see code below to decide this)
ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 6, EndMonth = 12, StartDay = 29, EndDay = 31, 
                                   CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 2, AnchorLag = -2,
                                   TimeUnit = 'd') # ts date and ts season same length

# Now add in site and count columns with count values as zeros where M_SEASON == 1 and NA 
# for all other days
ts_season_visit <- rbms::ts_monit_site(ts_season, gbs_count)

pheno <- readRDS("Data/Abundance data/Flight_curves_GBS_daily_anchor.rds")
# find when HB is at lowest abundance between the two flight curves to quantitatively determine when to start
# calculating abundance for 2nd flight period
pheno_hb <- pheno[pheno$SPECIES=="Holly Blue",]
pheno_hb <- pheno_hb[pheno_hb$trimDAYNO>150 & pheno_hb$trimDAYNO<230,]
min_abund <- pheno_hb %>% group_by(M_YEAR) %>% filter(NM==min(NM)) %>% ungroup()
min_abund <- unique(min_abund[,c(4,11,17)])
min_abund <- min_abund %>% group_by(M_YEAR) %>% filter(trimDAYNO==min(trimDAYNO)) %>% ungroup()
mean(min_abund$trimDAYNO) # 179.5 (round up to day 180) - start filter on 29th June

flight_curve <- ggplot(na.omit(pheno_hb), aes(x=trimDAYNO, y=NM, colour=M_YEAR))+
  geom_line()+
  labs(x="Monitoring Day", y="Relative Abundance", colour="")+
  scale_x_continuous(breaks=seq(0,366, by=10)) +
  theme_classic()+
  theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
flight_curve

i <- "Holly Blue"
sindex_final <- NULL
#options(warn=2)

# extract site index for each site, year and species
# impute_count() function uses the count data generated from ts_season_count() function and the flight curves
# it looks for the phenology available to estimate and input missing values
# imputation are made on a daily basis
ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, gbs_count, sp = "Holly Blue") # adds in observed counts for species i
impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='d')
# no imputed counts when m_season = 0
# this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year,
# and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season

# From the imputed count, the site index can be calculated for each site, or with a filter that will only keep sites that have
# been monitored at least a certain proportion of the flight curve
sindex_temp <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.1) # not using minFC here - as that removes sites that look like
# they've only monitored a small proportion of flight curve, but they have monitored more, it's just because we're restricting to
# months sept-nov
# this removes site/year combinations where TOTAL_NM<0.1 - now 806 sites (some of these will be 0 - where HB was not recorded)
sindex_final <- rbind(sindex_temp, sindex_final) # this is the annual abundance indices for each site

# filter by length of M_YEAR - species need to have at least 6 years (2016-2021)
flight_curve_years <- pheno %>% na.omit() %>% group_by(SPECIES) %>% summarise(n_years=n_distinct(M_YEAR))
# 7 species do not have enough data to estimate a flight curve for each year
# Remove these species from sindex_final and annual_indices
sp_exclude <- subset(flight_curve_years, n_years<6, select=SPECIES)
sindex_final <- sindex_final[!sindex_final$SPECIES %in% sp_exclude$SPECIES,]
# add grid references back in
sindex_final <- merge(sindex_final, site_match, by="SITE_ID", all.x=TRUE) 
# remove zeros where a species was not observed at a site
sindex_final <- sindex_final[!sindex_final$SINDEX==0,]
length(unique(sindex_final$grid_reference)) # 645 sites
write.csv(sindex_final, file="Data/Abundance data/Site_index_GBS_daily_ivy_HollyBlue.csv", row.names=FALSE)

