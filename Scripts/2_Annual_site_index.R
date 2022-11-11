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
                                   CompltSeason = TRUE, Anchor = FALSE, TimeUnit = 'd') # ts date and ts season same length

# Now add site visits to the time series
ts_season_visit <- rbms::ts_monit_site(ts_season, gbs_count)

# create loop to calculate and plot flight curve, and calculate sindex and plot annual indices for each species, each year
species <- unique(gbs_data$common_name)
sindex_final <- NULL
pheno_final <- NULL
annual_indices <- NULL
#options(warn=2)
for(i in species) {
print(i)
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, gbs_count, sp = i)
  ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 100, MinVisit = 4, MinOccur = 2, 
                                        MinNbrSite = 5, MaxTrial = 3, GamFamily = 'poisson', SpeedGam = FALSE, 
                                        CompltSeason = TRUE, SelectYear = 2019, TimeUnit = 'd')
  
  # extract pheno object: contains the shape of annual flight curves, standardised to sum 1
  pheno <- ts_flight_curve$pheno
  pheno$sp <- i
  pheno_final <- rbind(pheno, pheno_final, fill=TRUE)
  
  # temp <- ggplot(pheno, aes(x=trimDAYNO, y=NM, colour=M_YEAR))+
  #   geom_line()+
  #   labs(x="Monitoring Day", y="Relative Abundance", colour="")+
  #   ggtitle(i)+
  #   scale_x_continuous(breaks=seq(0,366, by=50)) +
  #   theme_classic()+
  #   theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  # ggsave(temp, file=paste0("../Graphs/Species flight curves/Flight_curve_", i,".png"), width = 20, height = 15, units = "cm")
  # Sys.sleep(2)
  
  # extract site index for each site, year and species
  # impute_count() function uses the count data generated from ts_season_count() function and the flight curves
  # it looks for the phenology available to estimate and input missing values
  # imputation can be made on a weekly or daily basis

  impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='d')

  # this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year,
  # and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season

  # From the imputed count, the site index can be calculated for each site, or with a filter that will only keep sites that have
  # been monitored at least a certain proportion of the flight curve

  sindex_temp <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.10)
  sindex_final <- rbind(sindex_temp, sindex_final)

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
  annual_indices <- rbind(b1, annual_indices)
  
  # annual_index <- ggplot(b1, aes(x=M_YEAR, y=LCI)) +
  #   geom_point()+
  #   geom_line()+
  #   geom_hline(yintercept=2, linetype='dotted', col = 'black')+
  #   ggtitle(i)+
  #   labs(x="Year", y=expression('log '['(10)']*' Collated Index'))+
  #   theme(title = element_text(size = 12))+
  #   theme_classic()
  # # mean total butterfly count expected across gardens each year
  # ggsave(annual_index, file=paste0("../Graphs/Species collated indices/Collated_index_", i,".png"), width = 30, height = 25, units = "cm")
  # Sys.sleep(2)


}
# Can't estimate flight curves for Small Blue - left with 30 species
options(op)

length(unique(sindex_final$SPECIES)) # 30 species
length(unique(annual_indices$sp)) # 30 species
length(unique(pheno_final$sp)) # 31 species

sindex_final <- merge(sindex_final, site_match, by="SITE_ID", all.x=TRUE)
write.csv(sindex_final, file="../Data/Site_index_GBS_daily.csv", row.names=FALSE)
write.csv(annual_indices, file="../Data/Annual_collated_index_daily.csv", row.names=FALSE)
saveRDS(pheno_final, file="../Data/Flight_curves_daily.rds")

# plot anchor vs no anchor flight curves
pheno_anchor <- read.csv("../Data/Flight_curves_daily_anchor.csv", header=TRUE)
pheno_noanchor <- read.csv("../Data/Flight_curves_daily_no_anchor.csv", header=TRUE)

pheno_anchor$anchor <- "anchor"
pheno_noanchor$anchor <- "no_anchor"

pheno <- rbind(pheno_anchor, pheno_noanchor)
pheno$M_YEAR <- as.factor(pheno$M_YEAR)

sp <- unique(pheno$SPECIES)

for(i in sp){
  print(i)

flight <- ggplot(pheno[pheno$SPECIES==i,], aes(x=trimDAYNO, y=NM, colour=M_YEAR))+
  geom_line()+
  labs(x="Monitoring Day", y="Relative Abundance", colour="")+
  ggtitle(i)+
  scale_x_continuous(breaks=seq(0,366, by=50)) +
  facet_grid(rows = vars(anchor))+
  theme_classic()+
  theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  
  ggsave(flight, file=paste0("../Graphs/Species flight curves/Anchor comparison/Flight_curve_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}

sindex <- read.csv("Data/Site_index_GBS.csv", header=TRUE)
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


for(i in species){print(i)
  
  temp_plot <- ggplot(all_indices[all_indices$sp==i,], aes(x=M_YEAR, y=LCI, colour=data))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=2, linetype='dotted', col = 'black')+
    ggtitle(i)+
    labs(x="Year", y=expression('log '['(10)']*' Collated Index'))+
    theme_classic()+
    theme(title = element_text(size = 14), legend.text=element_text(size=12),axis.text=element_text(size=12))
  ggsave(temp_plot, file=paste0("../Graphs/Species collated indices/Collated_index_comparison_", i,".png"), width = 20, height = 15, units = "cm")
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













































## add the line of the first year
yr <- unique(pheno[order(M_YEAR), as.numeric(as.character(M_YEAR))])

if("trimWEEKNO" %in% names(pheno)){
  plot(pheno[M_YEAR == yr[1], trimWEEKNO], pheno[M_YEAR == yr[1], NM], type = 'l', ylim = c(0, max(pheno[, NM])), xlab = 'Monitoring Week', ylab = 'Relative Abundance')
} else {
  plot(pheno[M_YEAR == yr[1], trimDAYNO], pheno[M_YEAR == yr[1], NM], type = 'l', ylim = c(0, max(pheno[, NM])), xlab = 'Monitoring Day', ylab = 'Relative Abundance')
  
}
## add individual curves for additional years
if(length(yr) > 1) {
  i <- 2
  for(y in yr[-1]){
    if("trimWEEKNO" %in% names(pheno)){
      points(pheno[M_YEAR == y , trimWEEKNO], pheno[M_YEAR == y, NM], type = 'l', col = i)
    } else {
      points(pheno[M_YEAR == y, trimDAYNO], pheno[M_YEAR == y, NM], type = 'l', col = i)
    }
    i <- i + 1
  }
}

## add legend
legend('topright', legend = c(yr), col = c(seq_along(c(yr))), lty = 1, bty = 'n')


############################################################


# 2. Inputting missing counts and producing annual collated indices


# impute_count() function uses the count data generated from ts_season_count() function and the flight curves
# it looks for the phenology available to estimate and input missing values
# imputation can be made on a weekly or daily basis

impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, ts_flight_curve=pheno, YearLimit= NULL, TimeUnit='d')

# this produces a data.table with original counts and imputed counts over the monitoring season, total count per site and year, 
# and total proportion of flight curve covered by the visits and SINDEX = sum of both observed and imputed counts over sampling season

# From the imputed count, the site index can be calculated for each site, or with a filter that will only keep sites that have 
# been monitored at least a certain proportion of the flight curve

sindex <- rbms::site_index(butterfly_count = impt_counts, MinFC = 0.10)

# Can also calculate an annual collated index by fitting GLMs where sites and years are modelled as factors
# proportion of flight curve sampled used as a weight for the GLM
# also remove all sites where the species was not observed (rm_zero=TRUE)

co_index <- collated_index(data = sindex, s_sp = "Pieris rapae", sindex_value = "SINDEX", glm_weights = TRUE, rm_zero = TRUE)

# transform index to a log(10) scale
co_index <- co_index$col_index
co_index_b <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
co_index_logInd <- co_index_b[BOOTi == 0, .(M_YEAR, COL_INDEX)][, log(COL_INDEX)/log(10), by = M_YEAR][, mean_logInd := mean(V1)]

## merge the mean log index with the full bootstrap dataset
data.table::setnames(co_index_logInd, "V1", "logInd"); setkey(co_index_logInd, M_YEAR); setkey(co_index_b, M_YEAR)
co_index_b <- merge(co_index_b, co_index_logInd, all.x = TRUE)

ggplot(co_index, aes(x=M_YEAR, y=COL_INDEX)) +
  geom_point()+
  geom_line()+
  labs(x="Year", y="Collated Index")+
  theme_classic()
# mean total butterfly count expected across gardens each year

# now can plot the log scaled indices
col_pal <- c("cyan4", "orange", "orangered2")

## compute the metric used for the graph of the Collated Log-Index centred around 2 (observed, bootstrap sample, credible confidence interval, linear trend)
b1 <- data.table(M_YEAR = co_index_b$M_YEAR, LCI = 2 + co_index_b$logInd - co_index_b$mean_logInd)
b2 <- data.table(M_YEAR = co_index_b[BOOTi == 0, M_YEAR], LCI= 2 + co_index_b[BOOTi == 0, logInd] - co_index_b[BOOTi == 0, mean_logInd])

lm_mod <- try(lm(LCI ~ M_YEAR, data = b2), silent=TRUE)

plot(b1, col = adjustcolor( "cyan4", alpha.f = 0.2),
     xlab = "year", ylab = expression('log '['(10)']*' Collated Index'),
     xaxt="n", type = 'n')

axis(1, at = b2$M_YEAR)

points(b2[!is.na(LCI),], type = 'l', lty=2, col="grey70")
points(b2, type = 'l', lwd=1.3, col = col_pal[1])
points(b2[!is.na(LCI),], col= col_pal[1], pch=19)
abline(h=2, lty=2)

if (class(lm_mod)[1] != "try-error"){
  points(b2$M_YEAR, as.numeric(predict(lm_mod, newdata = b2, type = "response")),
         type = 'l', col='maroon4', lwd = 1.5, lty = 1)
}

# Bootstrap confidence interval - n sites indices are randomly resampled (with replacement) k times to produce a 
# distribution of annual collated indices
set.seed(218795)
bootsample <- rbms::boot_sample(sindex, boot_n = 200)

co_index <- list()

## for progression bar, uncomment the following
## pb <- txtProgressBar(min = 0, max = dim(bootsample$boot_ind)[1], initial = 0, char = "*",  style = 3)

for(i in c(0,seq_len(dim(bootsample$boot_ind)[1]))){
  
  co_index[[i+1]] <- rbms::collated_index(data = sindex, s_sp = 2, sindex_value = "SINDEX", bootID=i, boot_ind= bootsample, glm_weights=TRUE, rm_zero=TRUE)
  
  ## for progression bar, uncomment the following
  ## setTxtProgressBar(pb, i)
  
}

## collate and append all the result in a data.table object
co_index <- rbindlist(lapply(co_index, FUN = "[[","col_index"))

# Annual log indices as well as their average are computed from the original sample
co_index_b <- co_index[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
co_index_logInd <- co_index_b[BOOTi == 0, .(M_YEAR, COL_INDEX)][, log(COL_INDEX)/log(10), by = M_YEAR][, mean_logInd := mean(V1)]

## merge the mean log index with the full bootstrap dataset
data.table::setnames(co_index_logInd, "V1", "logInd"); setkey(co_index_logInd, M_YEAR); setkey(co_index_b, M_YEAR)
co_index_b <- merge(co_index_b, co_index_logInd, all.x = TRUE)

data.table::setkey(co_index_b, BOOTi, M_YEAR)
co_index_b[ , boot_logInd := log(COL_INDEX)/log(10)]

# Now derive 95% confidence intervals
b1 <- data.table(M_YEAR = co_index_b$M_YEAR, LCI = 2 + co_index_b$boot_logInd - co_index_b$mean_logInd)
b2 <- data.table(M_YEAR = co_index_b[BOOTi == 0, M_YEAR], LCI= 2 + co_index_b[BOOTi == 0, logInd] - co_index_b[BOOTi == 0, mean_logInd])
b5 <- b1[co_index_b$BOOTi != 0, quantile(LCI, 0.025, na.rm = TRUE), by = M_YEAR]
b6 <- b1[co_index_b$BOOTi != 0, quantile(LCI, 0.975, na.rm = TRUE), by = M_YEAR]
lm_mod <- try(lm(LCI ~ M_YEAR, data = b2), silent=TRUE)

## define graph axis limits and color scheme
yl <- c(floor(min(b5$V1, na.rm=TRUE)), ceiling(max(b6$V1, na.rm=TRUE)))

col_pal <- c("cyan4", "orange", "orangered2")

## draw the plot for the selected species
plot(b1, ylim = yl, col = adjustcolor( "cyan4", alpha.f = 0.2),
     xlab = "year", ylab = expression('log '['(10)']*' Collated Index'),
     xaxt="n", type = 'n')

axis(1, at = b2$M_YEAR)

segments(x0 = as.numeric(unlist(b5[,1])), y0 = as.numeric(unlist(b5[,2])),
         x1 = as.numeric(unlist(b6[,1])), y1 = as.numeric(unlist(b6[,2])),
         col = col_pal[2], lwd = 2)
points(b2[!is.na(LCI),], type = 'l', lty=2, col="grey70")
points(b2, type = 'l', lwd=1.3, col = col_pal[1])
points(b2[!is.na(LCI),], col= col_pal[1], pch=19)
abline(h=2, lty=2)

if (class(lm_mod)[1] != "try-error"){
  points(b2$M_YEAR, as.numeric(predict(lm_mod, newdata = b2, type = "response")),
         type = 'l', col='maroon4', lwd = 1.5, lty = 1)
}


# plot mean abundance for sp=2 for each year
count_summary <- ts_season_count %>% group_by(WEEK, YEAR) %>%
  summarise(mean = mean(COUNT, na.rm=TRUE),
            lower = quantile(COUNT, 0.05, na.rm = TRUE),
            upper = quantile(COUNT, 0.95, na.rm = TRUE))

count_summary <- na.omit(count_summary)
count_summary$YEAR <- as.factor(count_summary$YEAR)
ggplot(count_summary, aes(WEEK, mean, color=YEAR))+
  geom_point(size = 2)+
  #geom_errorbar(aes(ymin = lower, ymax = upper))+
  geom_line()+
  ylab("Count")+
  xlab("Occasion (week)")+
  theme_classic()+
  theme(text = element_text(size = 18))
  #scale_x_continuous(breaks = seq(0,25,5)) 

# looks fairly similar to other plot







