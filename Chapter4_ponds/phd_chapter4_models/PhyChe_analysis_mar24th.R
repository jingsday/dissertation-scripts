rm(list=ls())
setwd('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Analysis')

#####################################################################################################

Pond <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/D2.4_PondCharacteristics_20240607.csv', header=T, sep = ",")
id <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/D2.4_PondID_20240124.csv', header=T, sep = ",")
PC <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/D2.4_PhysicoChemistry_20240905.csv', header=T, sep = ";")
LU5 <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/D2.4_LandUse_5m_20240727.csv', header=T, sep = ";")
LU100 <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/D2.4_LandUse_100m_20240724.csv', header=T, sep = ";")
LU <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/LU_current.csv', header=T, sep = ",")
LC <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/LandCover_selectedradii.csv', header=T, sep = ",")
A <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/Animals.csv', header=T, sep = ",")
Hetog100 <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/landscape_heterogeneity_100_allponds.csv', header=T, sep = ",")
Hetog5 <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/landscape_heterogeneity_5_allponds.csv', header=T, sep = ",")
HA <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/D2.4_HumanActivity_20240205.csv', header=T, sep = ",")
HP <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/Hydroperiod_length_data_sampling_resampling.csv', header=T, sep = ",")
ecels <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/ECELS GEA_2024-03-21.csv', header=T, sep = ",")
climate <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/Current_climate.csv', header=T, sep = ",")
fish <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/D2.2_Fish_20230711.csv', header=T, sep = ",")
#temp <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/temp.csv', header=T, sep = ",")
#prec <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/precipitation.csv', header=T, sep = ",")
#clim <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/climate_predictors_PONDERFUL.csv', header=T, sep = ",")

#notes fitxer PC
# Measurement date=NA: algunes ponds no tenen cap data, eliminar
# Algunes ponds dades gasos i dades de resampling (més d'una data per any)
# "<" límit de detecció com a valor
# Chla Dinamarca = CHLa_RFU (fluorescència)
# dades de nutrients pel mostreig gasos, biodiversitat i resampling. A l'hora de fer el merge amb biodiversitat i altres temes, fer-ho per data més propera

# Pond variables
id <- id[,c("Pond_ID","PondCode","X","Y","Pondscape","Country","Strat_Survey","Resampling")]
id = id[id$Strat_Survey==1,]
Pond <- Pond[,c("PondCode","AssessmentDate","Nat_res","Area","Depth","Pond_dries","Rel_waterlev","PVI")]
Pond$year <-substr(Pond$AssessmentDate, 1, 4)

# Physicochemical variables: TN, TP, O2....
PC[] <- lapply(PC, gsub, pattern='<', replacement='')
PC$NOTES <- NULL
PC$Notes..different.date.O2.spls.taken <- NULL
PC <-PC[!is.na(PC$MeasurementDate),]
PC[,-c(1:3)] <- sapply(PC[,-c(1:3)],as.numeric)
PC$year=as.integer(substr(PC[,'MeasurementDate'], 1, 4))
PC[!is.na(PC$TN) & PC$TN >= 10,]$TN <- NA
PC[!is.na(PC$CHLa_Spectro) & PC$CHLa_Spectro >= 800,]$CHLa_Spectro <- NA

#PC table seasons
PCs=PC[!is.na(PC$TN),c(1:3,12,13,38)]
# Count the number of unique seasons for each pond
library(dplyr)
pond_season_count <- PCs %>%
  distinct(PondCode, season_col) %>%  # Remove duplicate seasons for the same pond
  group_by(PondCode) %>%
  summarise(seasons_sampled = n_distinct(season_col)) %>%
  ungroup()
# Summarize the number of ponds sampled in 1, 2, or 3 seasons
season_summary <- pond_season_count %>%
  group_by(seasons_sampled) %>%
  summarise(number_of_ponds = n(),
            ponds = list(PondCode))
# Filter the ponds that have been sampled in all three seasons
ponds_sampled_in_all_three <- pond_season_count %>%
  filter(seasons_sampled == 3) %>%
  pull(PondCode)  # Extract the pond IDs

#######select the first date for each PC
library(dplyr)
PC$MeasurementDate <-as.Date(PC$MeasurementDate,"%Y-%m-%d")

# Keep only the row with the first measurement date where TN does not have NA
#PCunique <- PC %>%  group_by(PondCode) %>%
# filter(!is.na(TN)) %>%   filter(MeasurementDate == min(MeasurementDate))

# Mean of the spring measurements date where TN does not have NA
PCunique <- PC[,-2] %>% filter(season_col == "Spring") %>% group_by(PondCode,season_col) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Annual mean
PCannual <- PC[,-c(2,3)] %>% group_by(PondCode) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

# Mean of the hydroperiod
HPunique <- HP %>%  group_by(PondCode) %>%
  summarise(Hydeoperiod_length = mean(Hydeoperiod_length, na.rm = TRUE), .groups = "drop")

# Mean of the Pond
Pond$PVI = as.numeric(Pond$PVI)
Pondunique <- Pond %>%  group_by(PondCode) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

#fish
fish$FishS <- rowSums(fish[, c(5:55)])
fish$Fish =0
fish[fish$FishS>0,]$Fish <- 1
fish2=fish[,c("PondCode","Fish")]
fish22 <- fish2 %>%  group_by(PondCode) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

############################

#Climate variables
# temp$monthT = rowMeans(temp[,9:39], na.rm = T)
# meanT = aggregate(temp$monthT, list(temp$Pond_ID,temp$year), FUN=mean) 
# names(meanT)<-c("Pond_ID","year","T1")
# meanT$T1 = meanT$T1 - 273.15
# 
# prec$monthP = rowSums(prec[,9:39], na.rm = T)
# Prec = aggregate(prec$monthP, list(prec$Pond_ID,prec$year), FUN=sum) 
# names(Prec)<-c("Pond_ID","year","P1")
# 
# prec$monthP = rowMeans(prec[,9:39], na.rm = T)
# meanP = aggregate(prec$monthP, list(prec$Pond_ID,prec$year), FUN=mean) 
# names(meanP)<-c("Pond_ID","year","P1mean")
# 
# clim = clim[,c(4,6:19)]
# names(clim)[4]<-"T40"
# names(clim)[5]<-"P40"

PC1=merge(id,PCannual,by=c("PondCode"))
PC11=merge(PC1,HPunique,by=c("PondCode"), all.x=T)
#PC2=merge(PC11,meanT,by=c("Pond_ID","year"), all.x=T)
#PC3=merge(PC2,meanP,by=c("Pond_ID","year"), all.x=T)
#PC33=merge(PC3,Prec,by=c("Pond_ID","year"), all.x=T)
#PC333=merge(PC33,clim,by=c("PondCode"), all.x=T) #I do the merge with PondCode because there is some weird error with pondID
PC3=merge(PC11,climate[,c(2,7:14)],by=c("PondCode"), all.x=T) 

#Pond characteristics
PC4=merge(PC3,Pondunique,by=c("PondCode"), all.x=T)

#fish
PC4f=merge(PC4,fish22,by=c("PondCode"), all.x=T)

#Livestock
A$LivestockQ=1
A[A$livestock == "Moderate",]$LivestockQ<-2
A[A$livestock == "High",]$LivestockQ<-3

PC5=merge(PC4f,A[,c(1,2,5)],by=c("Pond_ID"))

# Land use 5m and 100m (field data)
LU5$Natural_5 = LU5$Moorland + LU5$Rank_veg + LU5$Woodland + LU5$Marsh + LU5$Other
LU5$year <-substr(LU5$AssessmentDate, 1, 4)
LU5unique <- LU5 %>%  group_by(PondCode) %>%
  summarise(Natural_5 = mean(Natural_5, na.rm = TRUE), .groups = "drop")



#LU100$Cultivated_land_100  = LU100$Arable + LU100$Imp_grassland

LU100$Arable_100 = LU100$Arable
LU100$Imp_grassland_100 = LU100$Imp_grassland
LU100$Pasture_Open_nature_100 = LU100$Moorland + LU100$Rank_veg + LU100$Rock + LU100$Unimp_grassland + LU100$Semimp_grassland
LU100$Forest_100 = LU100$Woodland
LU100$Other_100 = LU100$Other + LU100$Marsh + LU100$Streams + LU100$Urban

LU100$year <-substr(LU100$AssessmentDate, 1, 4)

LU100unique <- LU100 %>%  group_by(PondCode) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))


PC5LU=merge(PC5,LU5unique,by=c("PondCode"), all.x=T)
PC6=merge(PC5LU,LU100unique[,c(1,16:20)],by=c("PondCode"), all.x=T)

# Land use 3000m
names(LU)<-c("Pond_ID","PondCode","X","Y","Pondscape","Country","urban_3000","crop_3000","pasture_3000","forestry_3000","undefined_3000","natural_3000")
PC6d=merge(PC6,LU[,c(2,7:12)],by=c("PondCode"), all.x=T)  #I do the merge with PondCode because there is some weird error with pondID

# Land cover
PC6e=merge(PC6d,LC[,c(1,6:20)],by=c("Pond_ID"), all.x=T) 

# ECELS
PC7=merge(PC6e,ecels[,c(2,10)],by=c("Pond_ID"), all.x=T)


#De quins països no tenim dades de nutrients?
library(dplyr)
summary_data <- PC7 %>%
  group_by(Country) %>%
  summarize(across(c(17, 18, 26:33), ~sum(!is.na(.)), .names = "count_{.col}"))

summary(PC7[PC7$Country == "Germany",]$Pond_dries)

#######################################################################################################
#### Table 1
#######################################################################################################
Table1 = PC7[,c(9:13,15,17:19,25:30,54:55,56,44,45:52,74:78,60,89)]
Table1[,28:32]=Table1[,28:32]*100
Table1=Table1[Table1$Depth>0,]

colnames(Table1) <-c("pH","O2 con","O2 sat","Turbidity","Sludge","Conductivity","TN","TP","DOC","Alkalinity","Ca","Mg","Na","K","TSS","Area","Depth",
                     "Temporality","Hydeoperiod length","Annual Mean Temperature","Temperature Seasonality","Max Temperature of Warmest Month","Min Temperature of Coldest Month","Temperature Annual Range",
                     "Annual Precipitation","Precipitation Seasonality","Precipitation of Driest Quarter","Aquatic","Cropland","Forest","Pastures.and.open.nature","Urban","Livestock intensity","ECELS") 
            


calculate_summary <- function(data) {
  summary_table <- data.frame(
    Variable = character(),
    `Mean (SD)` = character(),
    `Median (Q1/Q3)` = character(),
    Minimum = numeric(),
    Maximum = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var_name in names(data)) {
    variable <- data[[var_name]]
    mean_sd <- sprintf("%.2f (%.2f)", mean(variable, na.rm = TRUE), sd(variable, na.rm = TRUE))
    median_q1q3 <- sprintf("%.2f (%.2f/%.2f)", median(variable, na.rm = TRUE),
                           quantile(variable, 0.25, na.rm = TRUE),
                           quantile(variable, 0.75, na.rm = TRUE))
    min_value <- min(variable, na.rm = TRUE)
    max_value <- max(variable, na.rm = TRUE)
    
    summary_table <- rbind(summary_table, data.frame(
      Variable = var_name,
      `Mean (SD)` = mean_sd,
      `Median (Q1/Q3)` = median_q1q3,
      Minimum = min_value,
      Maximum = max_value,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_table)
}

# Generate the summary table
summary_Table1 <- calculate_summary(Table1)
summary_Table1$Minimum <- round(summary_Table1$Minimum, 2)
summary_Table1$Maximum <- round(summary_Table1$Maximum, 2)

# Create a Word document and add the table
library(officer)
doc_Table1 <- read_docx() %>%
  body_add_table(value = summary_Table1, style = "Normal") %>%
  body_add_par("Summary Statistics Table", style = "heading 1")

# Save the document to a file
print(doc_Table1, target = "summary_Table1.docx")
write.csv(summary_Table1, file = "summary_Table1.csv", row.names = FALSE)


#######################################################################################################
#### Table 2
#######################################################################################################

#SEPARATE PONDS WITHOUT INFLUENCE OF LAND USE 0% CROPLAND URBAN

Table2 = PC7[,c(9:13,15,17:19,25:30,54:55,56,44,45:52,74:78,60)]
Table2[,28:31]=Table2[,28:31]*100
Table2=Table2[Table2$Depth>0,]
Table2=Table2[Table2$Cropland<5 & Table2$Urban<5 & Table2$Animals_cont<2,]


colnames(Table2) <-c("pH","O2 con","O2 sat","Turbidity","Sludge","Conductivity","TN","TP","DOC","Alkalinity","Ca","Mg","Na","K","TSS","Area","Depth",
                     "Temporality","Hydeoperiod length","Annual Mean Temperature","Temperature Seasonality","Max Temperature of Warmest Month","Min Temperature of Coldest Month","Temperature Annual Range",
                     "Annual Precipitation","Precipitation Seasonality","Precipitation of Driest Quarter","Aquatic","Cropland","Forest","Pastures.and.open.nature","Urban","Livestock intensity") 



calculate_summary <- function(data) {
  summary_table <- data.frame(
    Variable = character(),
    `Mean (SD)` = character(),
    `Median (Q1/Q3)` = character(),
    Minimum = numeric(),
    Maximum = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var_name in names(data)) {
    variable <- data[[var_name]]
    mean_sd <- sprintf("%.2f (%.2f)", mean(variable, na.rm = TRUE), sd(variable, na.rm = TRUE))
    median_q1q3 <- sprintf("%.2f (%.2f/%.2f)", median(variable, na.rm = TRUE),
                           quantile(variable, 0.25, na.rm = TRUE),
                           quantile(variable, 0.75, na.rm = TRUE))
    min_value <- min(variable, na.rm = TRUE)
    max_value <- max(variable, na.rm = TRUE)
    
    summary_table <- rbind(summary_table, data.frame(
      Variable = var_name,
      `Mean (SD)` = mean_sd,
      `Median (Q1/Q3)` = median_q1q3,
      Minimum = min_value,
      Maximum = max_value,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_table)
}

# Generate the summary table
summary_Table2 <- calculate_summary(Table2)
summary_Table2$Minimum <- round(summary_Table2$Minimum, 2)
summary_Table2$Maximum <- round(summary_Table2$Maximum, 2)

# Create a Word document and add the table
library(officer)
doc_Table2 <- read_docx() %>%
  body_add_table(value = summary_Table2, style = "Normal") %>%
  body_add_par("Summary Statistics Table", style = "heading 1")

# Save the document to a file
print(doc_Table2, target = "summary_Table2.docx")
write.csv(summary_Table2, file = "summary_Table2.csv", row.names = FALSE)

#######################################################################################################
#### Table 3
#######################################################################################################

#SEPARATE PONDS WITHOUT INFLUENCE OF LAND USE >50% CROPLAND URBAN >50

Table3 = PC7[,c(9:13,15,17:19,25:30,54:55,56,44,45:52,74:78,60)]
Table3[,28:31]=Table3[,28:31]*100
Table3=Table3[Table3$Depth>0,]
Table3=Table3[Table3$Cropland>50 | Table3$Animals_cont>2,]

colnames(Table3) <-c("pH","O2 con","O2 sat","Turbidity","Sludge","Conductivity","TN","TP","DOC","Alkalinity","Ca","Mg","Na","K","TSS","Area","Depth",
                     "Temporality","Hydeoperiod length","Annual Mean Temperature","Temperature Seasonality","Max Temperature of Warmest Month","Min Temperature of Coldest Month","Temperature Annual Range",
                     "Annual Precipitation","Precipitation Seasonality","Precipitation of Driest Quarter","Aquatic","Cropland","Forest","Pastures.and.open.nature","Urban","Livestock intensity") 



calculate_summary <- function(data) {
  summary_table <- data.frame(
    Variable = character(),
    `Mean (SD)` = character(),
    `Median (Q1/Q3)` = character(),
    Minimum = numeric(),
    Maximum = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var_name in names(data)) {
    variable <- data[[var_name]]
    mean_sd <- sprintf("%.2f (%.2f)", mean(variable, na.rm = TRUE), sd(variable, na.rm = TRUE))
    median_q1q3 <- sprintf("%.2f (%.2f/%.2f)", median(variable, na.rm = TRUE),
                           quantile(variable, 0.25, na.rm = TRUE),
                           quantile(variable, 0.75, na.rm = TRUE))
    min_value <- min(variable, na.rm = TRUE)
    max_value <- max(variable, na.rm = TRUE)
    
    summary_table <- rbind(summary_table, data.frame(
      Variable = var_name,
      `Mean (SD)` = mean_sd,
      `Median (Q1/Q3)` = median_q1q3,
      Minimum = min_value,
      Maximum = max_value,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_table)
}

# Generate the summary table
summary_Table3 <- calculate_summary(Table3)
summary_Table3$Minimum <- round(summary_Table3$Minimum, 2)
summary_Table3$Maximum <- round(summary_Table3$Maximum, 2)

# Create a Word document and add the table
library(officer)
doc_Table3 <- read_docx() %>%
  body_add_table(value = summary_Table3, style = "Normal") %>%
  body_add_par("Summary Statistics Table", style = "heading 1")

# Save the document to a file
print(doc_Table3, target = "summary_Table3.docx")
write.csv(summary_Table3, file = "summary_Table3.csv", row.names = FALSE)



#################################################################################################################################################
# Data transformation to achieve normality 
#################################################################################################################################################

# this was done many months ago so there might be some variables missing the normalisation
# so far, most transformations do not improve the normality of the variable. Important to check the normality of residuals in the models..

#normalitzem
library(bestNormalize)
shapiro.test(PC7$bio1)
bio1n <- bestNormalize(PC7$bio1)
plot(bio1n, leg_loc = "bottomright")
PC7$bio1.t <- predict(bio1n, newdata = PC7$bio1)
hist(PC7$bio1.t)
shapiro.test(PC7$bio1.t)

shapiro.test(PC7$bio4)
bio4n <- bestNormalize(PC7$bio4)
plot(bio4n, leg_loc = "bottomright")
PC7$bio4.t <- predict(bio4n, newdata = PC7$bio4)
hist(PC7$bio4.t)
shapiro.test(PC7$bio4.t)

shapiro.test(PC7$bio5)
bio5n <- bestNormalize(PC7$bio5)
plot(bio5n, leg_loc = "bottomright")
PC7$bio5.t <- predict(bio5n, newdata = PC7$bio5)
hist(PC7$bio5.t)
shapiro.test(PC7$bio5.t)

shapiro.test(PC7$bio6)
bio6n <- bestNormalize(PC7$bio6)
plot(bio6n, leg_loc = "bottomright")
PC7$bio6.t <- predict(bio6n, newdata = PC7$bio6)
hist(PC7$bio6.t)
shapiro.test(PC7$bio6.t)

shapiro.test(PC7$bio7)
bio7n <- bestNormalize(PC7$bio7)
plot(bio7n, leg_loc = "bottomright")
PC7$bio7.t <- predict(bio7n, newdata = PC7$bio7)
hist(PC7$bio7.t)
shapiro.test(PC7$bio7.t)

shapiro.test(PC7$bio12)
bio12n <- bestNormalize(PC7$bio12)
plot(bio12n, leg_loc = "bottomright")
PC7$bio12.t <- predict(bio12n, newdata = PC7$bio12)
hist(PC7$bio12.t)
shapiro.test(PC7$bio12.t)

shapiro.test(PC7$bio15)
bio15n <- bestNormalize(PC7$bio15)
plot(bio15n, leg_loc = "bottomright")
PC7$bio15.t <- predict(bio15n, newdata = PC7$bio15)
hist(PC7$bio15.t)
shapiro.test(PC7$bio15.t)

shapiro.test(PC7$bio17)
bio17n <- bestNormalize(PC7$bio17)
plot(bio17n, leg_loc = "bottomright")
PC7$bio17.t <- predict(bio17n, newdata = PC7$bio17)
hist(PC7$bio17.t)
shapiro.test(PC7$bio17.t)

shapiro.test(PC7$Area)
Arean <- bestNormalize(PC7$Area)
plot(Arean, leg_loc = "bottomright")
PC7$Area.t <- predict(Arean, newdata = PC7$Area)
hist(PC7$Area.t)
shapiro.test(PC7$Area.t)

shapiro.test(PC7$Depth)
Depthn <- bestNormalize(PC7$Depth)
plot(Depthn, leg_loc = "bottomright")
PC7$Depth.t <- predict(Depthn, newdata = PC7$Depth)
hist(PC7$Depth.t)
shapiro.test(PC7$Depth.t)

shapiro.test(PC7$Animals_cont)
Animals_contn <- bestNormalize(PC7$Animals_cont)
plot(Animals_contn, leg_loc = "bottomright")
PC7$Animals_cont.t <- predict(Animals_contn, newdata = PC7$Animals_cont)
hist(PC7$Animals_cont.t)
shapiro.test(PC7$Animals_cont.t)

shapiro.test(PC7$Hydeoperiod_length)
Hydeoperiod_lengthn <- bestNormalize(PC7$Hydeoperiod_length)
plot(Hydeoperiod_lengthn, leg_loc = "bottomright")
PC7$Hydeoperiod_length.t <- predict(Hydeoperiod_lengthn, newdata = PC7$Hydeoperiod_length)
hist(PC7$Hydeoperiod_length.t)
shapiro.test(PC7$Hydeoperiod_length.t)

shapiro.test(PC7$ECELS)
ECELSn <- bestNormalize(PC7$ECELS)
plot(ECELSn, leg_loc = "bottomright")
PC7$ECELS.t <- predict(ECELSn, newdata = PC7$ECELS)
hist(PC7$ECELS.t)
shapiro.test(PC7$ECELS.t)

ECELSn <- bestNormalize(PC7$ECELS)
PC7$ECELS.t <- predict(ECELSn, newdata = PC7$ECELS)

# Transform with orderNorm all landuse categories
library(dplyr)
cols_to_normalize <- 62:88
normalize_results <- lapply(PC7[, cols_to_normalize], bestNormalize)



hist(PC7$Cropland_500.t)
shapiro.test(PC7$Cropland_500.t)
hist(PC7$Arable_100.t)



############### Other way to transformar land use data
# quantile normalization
###############
# transformar_columna <- function(x) {
#   perc <- rank(x, na.last = "keep") / (sum(!is.na(x)) + 1)  # Calcular los percentiles manteniendo los NA
#   result <- ifelse(!is.na(x), qnorm(perc), NA)  # Transformar los percentiles a la distribución normal, manteniendo los NA
#   return(result)
# }
# PC7_T <- as.data.frame(apply(PC7[,73:90], 2, transformar_columna))
# colnames(PC7_T) <- paste0(colnames(PC7_T), ".t")
# 
# library(dplyr)
# normalityVar1<-PC7_T %>%
#   summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
#                              p.value = shapiro.test(.)$p.value))
# PC72=cbind(PC7,PC7_T)




#PCA Landuse 5 
# s'ha de refer, PC7[,73:80] no actualitzat, s'han afegit variables
# library(vegan)
# pca5 <- rda(PC7[,73:80], scale=T)
# biplot(pca5, display = c("sites", "species"),
#        type = c("text","points"))
# 
# pca5.t <- rda(na.omit(PC7[,105:112]), scale=T)
# biplot(pca5.t, display = c("sites", "species"),
#        type = c("text","points"))
# #PCA Landuse 100
# library(vegan)
# pca100 <- rda(PC7[,82:89], scale=T)
# biplot(pca100, display = c("sites", "species"),
#        type = c("text","points"))
# 
# pca100.t <- rda(PC7[,114:121], scale=T)
# biplot(pca100.t, display = c("sites", "species"),
#        type = c("text","points"))






#correlations
# 
library(psych) 
png("correlation_PC.t.png", height=15, width=15, units = "in", res = 300)
pairs.panels(PC7[,c(94:106)], scale=T)
dev.off()

# ALWAYS: just use the climatic variables: bio1, bio7 and bio12 (forget about all other climatic variables)
# General: check correlation forest vs. bio12. If >0.6, remove forest
# Temperate: remove forest and bio7
# Mediterranean: remove forest, bio12, bio7
# Continental: remove pastures, aquatic500, area; and check correlation cropland and bio1, cropland and bio1 and animals_cont and bio1. If>0.6 keep just 1: animals_cont?
# Subtropical: remove forest, bio12, bio7
############################################################################################
# PCA nutrients
############################################################################################
library(ggplot2)
library(ggfortify)
library(factoextra)
nutr=na.omit(PC7[,c("PondCode","Country","TN","TP","Cond","pH","O2_con")])

n_pca <- prcomp(nutr[,-c(1:2)], scale=TRUE)
nutrPC <- as.data.frame(n_pca$x)

nutr=cbind(nutr,nutrPC[,1:2])
names(nutr)[8:9]<-c("PC_PC1","PC_PC2")
PC7 = merge(PC7,nutr[,c(1,8:9)],by="PondCode", all.x=T)


write.csv(PC7,"C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/PhyChe_XY.csv")

png("PC_PCA.png", height=5, width=7, units = "in", res = 300)
fviz_pca_biplot(n_pca, repel = TRUE, col.var = "deepskyblue", geom="point")
dev.off()

png("PC_PCA_percvar.png", height=5, width=5, units = "in", res = 300)
fviz_eig(n_pca, addlabels = TRUE)
dev.off()
############################################################################################
# RDA nutrients
############################################################################################

library(vegan)
PC8=na.omit(PC7[,c("Country","TN","TP","Cond","pH","O2_con","T1.t","P1.t","Nat_res","Area.t","Depth.t","Pond_dries.t","Hetog5.t","Hetog100.t","LU_PCA.t","Animals_cont.t","LivestockQ.t","ECELS.t")])

en=PC8[,c("T1.t","P1.t","Nat_res","Area.t","Depth.t","Pond_dries.t","Hetog5.t","Hetog100.t","LU_PCA.t","Animals_cont.t","LivestockQ.t","ECELS.t")]
n=PC8[,c("TN","TP","Cond","pH","O2_con")]
dgen <- en+3
decorana(na.omit(n)) #fem RDA: if the length of the first axis in my data is less then 2 = RDA, more then 2 = CCA. Not a foolproof, but quite robust as using unimodal method with linear data is not a big problem, but using linear method for unimodal data could be

rda1 <- rda(n ~ T1.t+P1.t+Area.t+Depth.t+Pond_dries.t+Hetog5.t+Hetog100.t+LU_PCA.t+Animals_cont.t, data=PC8, scale=T)
summary(rda1)
anova(rda1)
anova(rda1, by="terms", permu=200)#term significance
RsquareAdj(rda1)
vif.cca(rda1) #High VIF values (typically greater than 10) suggest multicollinearity and may indicate that the corresponding variables are highly correlated with each other. It's generally recommended to address multicollinearity issues by either removing one of the correlated variables or by using dimensionality reduction techniques before fitting the CCA model.
(anova(rda1, by="terms", permu=200))$Variance
#All values are below 10, and most are below 5, which indicates that multicollinearity among these predictors shouldn’t be a problem for the model. 



png("rda_nutr.png", height=7, width=7, units = "in", res = 300)
library(RColorBrewer)
col8 <- brewer.pal(8, "Set1")

plot(rda1, type="n", scaling=3,xlim = c(-2, 3))
text(rda1, display="species", cex=1, col="gray32", scaling=3)           
points(rda1, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=col8[as.factor(PC8$Country)]) 
text(rda1, scaling=3, display="bp", col="black", cex=1)                           
legend("topleft", legend=levels(as.factor(PC8$Country)), bty="n", col="gray32", pch=21, cex=1, pt.bg=col8)
dev.off()




############################################################################################
# RDA 2 with ECELS
############################################################################################
PC8=na.omit(PC7[,c("Country","TN","TP","Cond","pH","O2_con","T1.t","P1.t","Nat_res","Area.t","Depth.t","Pond_dries.t","Hetog5.t","Hetog100.t","LU_PCA.t","Animals_cont.t","LivestockQ.t","ECELS.t")])

# !!!!!!! s'ha de refer, noves variables no actualitzades!!!!!!!!!!

en=PC8[,c("T1.t","P1.t","Nat_res","Area.t","Depth.t","Pond_dries.t","Hetog5.t","Hetog100.t","LU_PCA.t","Animals_cont.t","LivestockQ.t","ECELS.t")]
n=PC8[,c("TN","TP","Cond","pH","O2_con")]
dgen <- en+3
decorana(na.omit(n)) #fem RDA: if the length of the first axis in my data is less then 2 = RDA, more then 2 = CCA. Not a foolproof, but quite robust as using unimodal method with linear data is not a big problem, but using linear method for unimodal data could be

rda2 <- rda(n ~ T1.t+P1.t+Area.t+Depth.t+Pond_dries.t+Hetog5.t+Hetog100.t+LU_PCA.t+Animals_cont.t+ECELS.t, data=PC8, scale=T)
summary(rda2)
anova(rda2)
anova(rda2, by="terms", permu=200)#term significance
RsquareAdj(rda2)
vif.cca(rda2) #High VIF values (typically greater than 10) suggest multicollinearity and may indicate that the corresponding variables are highly correlated with each other. It's generally recommended to address multicollinearity issues by either removing one of the correlated variables or by using dimensionality reduction techniques before fitting the CCA model.
(anova(rda2, by="terms", permu=200))$Variance
#All values are below 10, and most are below 5, which indicates that multicollinearity among these predictors shouldn’t be a problem for the model. 



png("rda2_nutr_with_ECELS.png", height=7, width=7, units = "in", res = 300)
library(RColorBrewer)
col8 <- brewer.pal(8, "Set1")

plot(rda2, type="n", scaling=3,xlim = c(-2, 3))
text(rda2, display="species", cex=1, col="gray32", scaling=3)           
points(rda2, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=col8[as.factor(PC8$Country)]) 
text(rda2, scaling=3, display="bp", col="black", cex=1)                           
legend("topleft", legend=levels(as.factor(PC8$Country)), bty="n", col="gray32", pch=21, cex=1, pt.bg=col8)
dev.off()

############################################################################################
# RDA 3 pigments
############################################################################################
PC9=na.omit(PC7[,c("Country","TN","TP","Cond","pH","O2_con","CHLa","PC","T1.t","P1.t","Nat_res","Area.t","Depth.t","Pond_dries.t","Hetog5.t","Hetog100.t","LU_PCA.t","Animals_cont.t","LivestockQ.t","ECELS.t")])

# !!!!!!! s'ha de refer, noves variables no actualitzades!!!!!!!!!!


en7=PC9[,c("T1.t","P1.t","Nat_res","Area.t","Depth.t","Pond_dries.t","Hetog5.t","Hetog100.t","LU_PCA.t","Animals_cont.t","LivestockQ.t","ECELS.t")]
n7=PC9[,c("TN","TP","Cond","pH","O2_con","CHLa","PC")]
dgen7 <- en7+3
decorana(na.omit(n7)) #fem RDA: if the length of the first axis in my data is less then 2 = RDA, more then 2 = CCA. Not a foolproof, but quite robust as using unimodal method with linear data is not a big problem, but using linear method for unimodal data could be

rda3 <- rda(n7 ~ T1.t+P1.t+Area.t+Depth.t+Pond_dries.t+Hetog5.t+Hetog100.t+LU_PCA.t+Animals_cont.t+ECELS.t, data=PC9, scale=T)
summary(rda3)
anova(rda3)
anova(rda3, by="terms", permu=200)#term significance
RsquareAdj(rda3)
vif.cca(rda3) #High VIF values (typically greater than 10) suggest multicollinearity and may indicate that the corresponding variables are highly correlated with each other. It's generally recommended to address multicollinearity issues by either removing one of the correlated variables or by using dimensionality reduction techniques before fitting the CCA model.
(anova(rda3, by="terms", permu=200))$Variance
#All values are below 10, and most are below 5, which indicates that multicollinearity among these predictors shouldn’t be a problem for the model. 



png("rda3_pigments.png", height=7, width=7, units = "in", res = 300)
library(RColorBrewer)
col8 <- brewer.pal(8, "Set1")

plot(rda3, type="n", scaling=3,xlim = c(-2, 3))
text(rda3, display="species", cex=1, col="gray32", scaling=3)           
points(rda3, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=col8[as.factor(PC9$Country)]) 
text(rda3, scaling=3, display="bp", col="black", cex=1)                           
legend("topleft", legend=levels(as.factor(PC9$Country)), bty="n", col="gray32", pch=21, cex=1, pt.bg=col8)
dev.off()

####################################################################################
### NMDS land use
####################################################################################
library(vegan)
L100<-LU100[,-c(1:3,16:17)]
nmds <- metaMDS(L100, distance = "bray", k = 2)
stressplot(nmds)
plot(nmds, type = "n") # Create empty plot
text(nmds, display = "sites") # Add sample labels to the plot
points(nmds, col = "red") # Plot samples
nmds$stress #Lower stress values (<0.3) indicate better preservation of the original distances.
anova(nmds)
stress_value / nmds$GOF$Ms[1]

png("NMDS_ALL_BrayCurtis.png", height=7, width=7, units = "in", res = 300)
plot(nmds)
ordiplot(nmds,type="n")
orditorp(nmds,display="species",col="black",air=0.1)
orditorp(nmds,display="sites",air=0.01,cex=0.1)
dev.off()





####################################################################################
### LME
####################################################################################
# in the article https://esajournals.onlinelibrary.wiley.com/doi/10.1890/ES12-00116.1 there is a very clear explanation of the method

# 1. repassar correlacions entre variables. Amb aquests models correlacions entre variables >0.6 no funcionen
# 2. decidir si començar amb variables independents transformades o fer-ho després només si cal


#lme() because you can't include a correlation structure using lmer()
library(nlme)
library(arm)
library(car)
library(MASS)
library(MuMIn) #rsquared for MLM 


#qqnorm( All, ~ resid(., type = "n") )
#plot(Variogram( All, form = ~ latitude + longitude ))
#ML is necessary because comparisons using REML are not valid when the fixed effects change
#REML is generally considered to give better estimates for the random effects,
#though, so the usual advice is to fit your best model using REML for your final inference and reporting.

#this is just an example, all LU categories need to be added and decide which are the best predictor climatic variables
#full model
# ¿¿¿¿LU arable is not added because it a percentage and it will be explained by the intercept. But I'm not sure it is
# correct when there are other predictor variables that are not part of the %???
# Also, I need to check how to deal with this issue in the random effects


#correlation = corExp(form = ~ xlong + ylat, nugget=T) --> the spatial autocorrelation. If Uruguay is present, I wouldn't use it

# select best model or 3 best models
stepAIC(m1)

anova(m1,m1b)

# with the best model, in this case m1:
summary(m1) #Table 1, Jackson et al 2012 (fixed+random)
ranef(m1) #Table 3, Jackson et al 2012 (fixed+random)
var(m1$residuals)
intervals(m1)
r.squaredGLMM(m1)


################################################################################################################################################
# analyze data: lmer (MLM)
# compute random effects pvalues #Table 1, Jackson et al 2012
################################################################################################################################################

##### T1
m1 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
          random = list(~1|PondCode, ~1|Country,~ 0+Aable_100|Country,~ 0+Imp_grassland_100|Country,~ 0+Open_nature_100|Country,~ 0+Forest_100|Country,~ 0+Pasture_100|Country,~ 0+Urban_100|Country,~ 0+Other_100|Country),
          correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)

#Without Arable_100
m1_1 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
            random = list(~1|PondCode, ~1|Country,~ 0+Imp_grassland_100|Country,~ 0+Open_nature_100|Country,~ 0+Forest_100|Country,~ 0+Pasture_100|Country,~ 0+Urban_100|Country,~ 0+Other_100|Country),
            correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)

#Without Imp_grassland_100
m1_2 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
            random = list(~1|PondCode, ~1|Country,~ 0+Aable_100|Country,~ 0+Open_nature_100|Country,~ 0+Forest_100|Country,~ 0+Pasture_100|Country,~ 0+Urban_100|Country,~ 0+Other_100|Country),
            correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)

#Without Open_nature_100
m1_3 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
            random = list(~1|PondCode, ~1|Country,~ 0+Aable_100|Country,~ 0+Imp_grassland_100|Country,~ 0+Forest_100|Country,~ 0+Pasture_100|Country,~ 0+Urban_100|Country,~ 0+Other_100|Country),
            correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)

#Without Forest_100
m1_4 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
            random = list(~1|PondCode, ~1|Country,~ 0+Aable_100|Country,~ 0+Imp_grassland_100|Country,~ 0+Open_nature_100|Country,~ 0+Pasture_100|Country,~ 0+Urban_100|Country,~ 0+Other_100|Country),
            correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)

#Without Pasture_100
m1_5 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
            random = list(~1|PondCode, ~1|Country,~ 0+Aable_100|Country,~ 0+Imp_grassland_100|Country,~ 0+Open_nature_100|Country,~ 0+Forest_100|Country,~ 0+Urban_100|Country,~ 0+Other_100|Country),
            correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)

#Without Urban_100
m1_6 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
            random = list(~1|PondCode, ~1|Country,~ 0+Aable_100|Country,~ 0+Imp_grassland_100|Country,~ 0+Open_nature_100|Country,~ 0+Forest_100|Country,~ 0+Pasture_100|Country,~ 0+Other_100|Country),
            correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)

#Without Other_100
m1_7 <- lme(fixed = TP ~ Imp_grassland_100 +Open_nature_100 + Forest_100 + Pasture_100 + Urban_100 + Other_100 + Forest_5  + Animals_cont + Area + Depth + T1 + P1, 
            random = list(~1|PondCode, ~1|Country,~ 0+Aable_100|Country,~ 0+Imp_grassland_100|Country,~ 0+Open_nature_100|Country,~ 0+Forest_100|Country,~ 0+Pasture_100|Country,~ 0+Urban_100|Country),
            correlation = corExp(form = ~ xlong + ylat, nugget=T),method = "ML",data=PC7)


# if deviance in LMER is the same as (-2)*logLik, this is correct:
LLm1 <- c((-2)*logLik(m1)[1], (-2)*logLik(m1_1)[1], (-2)*logLik(m1_2)[1], (-2)*logLik(m1_3)[1], (-2)*logLik(m1_4)[1], (-2)*logLik(m1_5)[1], (-2)*logLik(m1_6)[1], (-2)*logLik(m1_7)[1])

mlm.pvalsm1 <- c(1-pchisq(LLm1[2] - LLm1[1],1), 1-pchisq(LLm1[3] - LLm1[1],1), 1-pchisq(LLm1[4] - LLm1[1],1), 1-pchisq(LLm1[5] - LLm1[1],1), 1-pchisq(LLm1[6] - LLm1[1],1), 1-pchisq(LLm1[7] - LLm1[1],1), 1-pchisq(LLm1[8] - LLm1[1],1))

mlm.pvalsm1 #Table1


# ens quedem amb random slope significatives, només





###########################################################
### Regression tree analysis
###########################################################
## Implementation of regression tree (ref. Gareth et al, 2021. An introduction to statistical learning)
## Classification and Regression Trees (CART) with rpart

library(rpart) #for fitting decision trees
library(rpart.plot) #for plotting decision trees

#db <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/full_table_first_measured_PCA_included.csv', header=T, sep = ",")

#### ALL VARIABLES

## Subset data and change columns' name
PC7_sel100<-PC7[,c("PC_PC1","Forest_5","Arable_100","Animals_cont","Depth","Area","Hydeoperiod_length","T1","P1")]
PC7_sel500<-PC7[,c("PC_PC1","Forest_5","Cropland_500","Animals_cont","Depth","Area","Hydeoperiod_length","T1","P1")]
PC7_sel3000<-PC7[,c("PC_PC1","Forest_5","crop_3000","Animals_cont","Depth","Area","Hydeoperiod_length","T1","P1")]
PC7_sel5000<-PC7[,c("PC_PC1","Forest_5","Cropland_5000","Animals_cont","Depth","Area","Hydeoperiod_length","T1","P1")]

PC7_sel=PC7_sel100[!is.na(PC7_sel100$PC_PC1),]
## Data split (train and test, 85% and 15%, respectively)
set.seed(1234) # for reproducibility
sample_set <- sample(nrow(PC7_sel), round(nrow(PC7_sel)*.85), replace = FALSE)
PC7_train <- PC7_sel[sample_set, ] #train dataset (80% whole PC7)
PC7_test <- PC7_sel[-sample_set, ] #test dataset (20% whole PC7)


## Fitting the regression tree with all data
mod1 <-
  rpart(
    (PC_PC1) ~ .,
    method = "anova",
    data = PC7_sel
  )

rpart.plot(mod1) #plot the first regression tree on the train dataset


## Fitting the regression tree with train data
mod1 <-
  rpart(
    (PC_PC1) ~ .,
    method = "anova",
    data = PC7_train
  )

rpart.plot(mod1) #plot the first regression tree on the train dataset

## Use the model on the test dataset and compute mse
pred.tree = predict(mod1, PC7_test)
mse <- mean(((pred.tree - PC7_test$PC_PC1)^2))
mse


## Try to prune the tree
printcp(mod1) #lowest xerror cp=0.252021
bestcp <- mod1$cptable[which.min(mod1$cptable[,"xerror"]),"CP"]
bestcp #confirm the same cp

pruned.tree <- prune(mod1, cp = bestcp)

## Use the model (pruned tree) on the test dataset and compute mse
pred.prune = predict(pruned.tree, PC7_test)
mse <- mean((pred.prune - PC7_test$PC_PC1)^2)
mse #mse is higher, so it is better the tree not pruned

rpart.plot(pruned.tree)


#### LAND USE SELECTION
#db <- read.csv('C:/Users/650475/OneDrive - Universitat de Vic/UVic/Projectes en curs/PONDERFUL/2 BASIC RESEARCH/PhyChe/Data/Unified data - Intranet/full_table_first_measured_PCA_included.csv', header=T, sep = ",")
#PC7 <-merge(db,PC7[,c("Pond_ID","Pastures.and.open.nature_5000","Urban_5000")],by="Pond_ID")

lu_sel100<-PC7[,c("PC_PC1","Arable_100","Imp_grassland_100","Open_nature_100","Pasture_100","Urban_100","Aquatic_100","Other_100")]
lu_sel500<-PC7[,c("PC_PC1","Aquatic_500","Cropland_500","Forest_500","Pastures.and.open.nature_500","Urban_500")]
lu_sel3000<-PC7[,c("PC_PC1","urban_3000","crop_3000","pasture_3000","forestry_3000","undefined_3000","natural_3000")]
lu_sel5000<-PC7[,c("PC_PC1","Aquatic_5000","Cropland_5000","Forest_5000","Pastures.and.open.nature_5000","Urban_5000")]

lu_sel=lu_sel100[!is.na(lu_sel100$PC_PC1),]
## Data split (train and test, 85% and 15%, respectively)
set.seed(1234) # for reproducibility
sample_set <- sample(nrow(lu_sel), round(nrow(lu_sel)*.85), replace = FALSE)
lu_train <- lu_sel[sample_set, ] #train dataset (80% whole lu)
lu_test <- lu_sel[-sample_set, ] #test dataset (20% whole lu)


## Fitting the regression tree
mod1 <-
  rpart(
    (PC_PC1) ~ .,
    method = "anova",
    data = lu_sel
  )

rpart.plot(mod1) #plot the first regression tree on the train dataset

## Use the model on the test dataset and compute mse
pred.tree = predict(mod1, lu_test)
mse <- mean(((pred.tree - lu_test$PC_PC1)^2))
mse


## Try to prune the tree
printcp(mod1) #lowest xerror cp=0.252021
bestcp <- mod1$cptable[which.min(mod1$cptable[,"xerror"]),"CP"]
bestcp #confirm the same cp

pruned.tree <- prune(mod1, cp = bestcp)

## Use the model (pruned tree) on the test dataset and compute mse
pred.prune = predict(pruned.tree, lu_test)
mse <- mean((pred.prune - lu_test$PC_PC1)^2)
mse #mse is higher, so it is better the tree not pruned





######################################
#### OUTLIERS
######################################

o <- PC7 %>%
  dplyr::select(PondCode,Pond_ID,Country, Pondscape,TN,TP, pH, O2_con, Turb,ECELS,Arable_100,Urban_100,Cropland_500,urban_3000,crop_3000,pasture_3000,forestry_3000,undefined_3000,natural_3000, Cropland_5000,Forest_5, 
                Animals_cont, Area, Depth, T1, P1,Hydeoperiod_length)

o10 = o[o$TN>=10,]
write.csv(o10,"Outliers_TN10.csv")


opc10 = PC[PC$TN>=10,]
write.csv(opc10,"Outliers_TN10_intraannual.csv")


opc10b=PC[PC$PondCode %in% c("FB1_IN1","QUI_ex5","QUI_18"), ]
write.csv(opc10b,"Outliers_TN10_3_all.csv")






############################################################################################
# Land use BAR PLOT TN TP
############################################################################################
summary(PC7$TN)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.3200  0.8652  1.3310  1.8090  2.2850  7.9700       1 

summary(PC7$TP)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.3200  0.8652  1.3310  1.8090  2.2850  7.9700       1 

TN_80p <- PC7 %>%
  group_by(Country) %>%
  summarize(Top_20_TN = quantile(TN, 0.80, na.rm = TRUE))  # 0.80 for the top 20%
DTN_80p <- PC7 %>%
  inner_join(TN_80p, by = "Country") %>%
  filter(TN >= Top_20_TN)

TN_20p <- PC7 %>%
  group_by(Country) %>%
  summarize(Bottom_20_TN = quantile(TN, 0.20, na.rm = TRUE))  # 0.20 for the bottom 20%
DTN_20p <- PC7 %>%
  inner_join(TN_20p, by = "Country") %>%
  filter(TN <= Bottom_20_TN)

TP_80p <- PC7 %>%
  group_by(Country) %>%
  summarize(Top_20_TP = quantile(TP, 0.80, na.rm = TRUE))  # 0.80 for the top 20%
DTP_80p <- PC7 %>%
  inner_join(TP_80p, by = "Country") %>%
  filter(TP >= Top_20_TP)

TP_20p <- PC7 %>%
  group_by(Country) %>%
  summarize(Bottom_20_TP = quantile(TP, 0.20, na.rm = TRUE))  # 0.20 for the bottom 20%
DTP_20p <- PC7 %>%
  inner_join(TP_20p, by = "Country") %>%
  filter(TP <= Bottom_20_TP)

library(dplyr)
# Calculate the mean of land use categories per country
TN80_LU <- DTN_80p %>%
  group_by(Country) %>%
  summarize(
    Mean_Aquatic_500 = mean(Aquatic_500, na.rm = TRUE),
    Mean_Forest_500 = mean(Forest_500, na.rm = TRUE),
    Mean_Pastures.and.open.nature_500 = mean(Pastures.and.open.nature_500, na.rm = TRUE),
    Mean_Cropland_500 = mean(Cropland_500, na.rm = TRUE),
    Mean_Urban_500 = mean(Urban_500, na.rm = TRUE)
  )
TN80_LU$Percentile = "Top 20th Percentile, TN"

TN20_LU <- DTN_20p %>%
  group_by(Country) %>%
  summarize(
    Mean_Aquatic_500 = mean(Aquatic_500, na.rm = TRUE),
    Mean_Forest_500 = mean(Forest_500, na.rm = TRUE),
    Mean_Pastures.and.open.nature_500 = mean(Pastures.and.open.nature_500, na.rm = TRUE),
    Mean_Cropland_500 = mean(Cropland_500, na.rm = TRUE),
    Mean_Urban_500 = mean(Urban_500, na.rm = TRUE)
  )
TN20_LU$Percentile = "Bottom 20th Percentile, TN"

TP80_LU <- DTP_80p %>%
  group_by(Country) %>%
  summarize(
    Mean_Aquatic_500 = mean(Aquatic_500, na.rm = TRUE),
    Mean_Forest_500 = mean(Forest_500, na.rm = TRUE),
    Mean_Pastures.and.open.nature_500 = mean(Pastures.and.open.nature_500, na.rm = TRUE),
    Mean_Cropland_500 = mean(Cropland_500, na.rm = TRUE),
    Mean_Urban_500 = mean(Urban_500, na.rm = TRUE)
  )
TP80_LU$Percentile = "Top 20th Percentile, TP"

TP20_LU <- DTP_20p %>%
  group_by(Country) %>%
  summarize(
    Mean_Aquatic_500 = mean(Aquatic_500, na.rm = TRUE),
    Mean_Forest_500 = mean(Forest_500, na.rm = TRUE),
    Mean_Pastures.and.open.nature_500 = mean(Pastures.and.open.nature_500, na.rm = TRUE),
    Mean_Cropland_500 = mean(Cropland_500, na.rm = TRUE),
    Mean_Urban_500 = mean(Urban_500, na.rm = TRUE)
  )
TP20_LU$Percentile = "Bottom 20th Percentile, TP"

TN80_A <- DTN_80p %>%
  group_by(Country) %>%
  summarize(Mean_Animals_cont = mean(Animals_cont, na.rm = TRUE))
TN80_A$Percentile = "Top 20th Percentile, TN"

TN20_A <- DTN_20p %>%
  group_by(Country) %>%
  summarize(Mean_Animals_cont = mean(Animals_cont, na.rm = TRUE))
TN20_A$Percentile = "Bottom 20th Percentile, TN"

TP80_A <- DTP_80p %>%
  group_by(Country) %>%
  summarize(Mean_Animals_cont = mean(Animals_cont, na.rm = TRUE))
TP80_A$Percentile = "Top 20th Percentile, TP"

TP20_A <- DTP_20p %>%
  group_by(Country) %>%
  summarize(Mean_Animals_cont = mean(Animals_cont, na.rm = TRUE))
TP20_A$Percentile = "Bottom 20th Percentile, TP"


TNTP_LU = rbind(TN80_LU,TN20_LU,TP80_LU,TP20_LU)
TNTP_LU$Percentile <- factor(TNTP_LU$Percentile, levels=c('Bottom 20th Percentile, TN', 'Top 20th Percentile, TN', 'Bottom 20th Percentile, TP', 'Top 20th Percentile, TP'))

TNTP_A = rbind(TN80_A,TN20_A,TP80_A,TP20_A)
TNTP_A$Percentile <- factor(TNTP_A$Percentile, levels=c('Bottom 20th Percentile, TN', 'Top 20th Percentile, TN', 'Bottom 20th Percentile, TP', 'Top 20th Percentile, TP'))
TNTP_A$Country <- factor(TNTP_A$Country, levels=c("Uruguay","Turkey","Spain","Switzerland","Belgium","Germany","UK","Denmark"))

# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Reshape the data from wide to long format
TNTP_LU_long <- TNTP_LU  %>%
  pivot_longer(cols = starts_with("Mean"),  # Select all columns starting with "Mean"
               names_to = "LandUseCategory", 
               values_to = "MeanPercentage") %>%
  mutate(LandUseCategory = recode(LandUseCategory,
                                  "Mean_Aquatic_500" = "Aquatic",
                                  "Mean_Forest_500" = "Forest",
                                  "Mean_Pastures.and.open.nature_500" = "Pastures and Open Nature",
                                  "Mean_Cropland_500" = "Cropland",
                                  "Mean_Urban_500" = "Urban"))
TNTP_LU_long$LandUseCategory <- factor(TNTP_LU_long$LandUseCategory, levels=c('Urban', 'Cropland', 'Pastures and Open Nature', 'Forest','Aquatic'))
TNTP_LU_long$Country <- factor(TNTP_LU_long$Country, levels=c("Uruguay","Turkey","Spain","Switzerland","Belgium","Germany","UK","Denmark"))

# Define the colors for each land use category
land_use_colors <- c(
  'Urban' = '#7f7f7f',                # Gray for Urban
  'Cropland' = '#ff7f0e',              # Orange for Cropland
  'Pastures and Open Nature' = '#bcbd22',  # Light green/yellow for Pastures
  'Forest' = '#2ca02c',                # Green for Forest
  'Aquatic' = '#1f77b4'                # Blue for Aquatic
)

# Create the bar plot LU
png("TNTP_LU.png", height=7, width=10, units = "in", res = 600)
ggplot(TNTP_LU_long, aes(x = Country, y = MeanPercentage, fill = LandUseCategory)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  coord_flip() +
  scale_fill_manual(values = land_use_colors) +
  facet_wrap(~Percentile) +
  theme_minimal(base_size = 15) +  # Clean theme with larger font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  theme(aspect.ratio = 1)  # Ensures a square plot (1:1 aspect ratio)
dev.off()



# Create the bar plot Livestock
png("TNTP_A.png", height=7, width=4, units = "in", res = 600)
ggplot(TNTP_A, aes(x = Country, y = Mean_Animals_cont)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  facet_wrap(~Percentile) +
  theme_minimal(base_size = 15) +  # Clean theme with larger font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) 
dev.off()




############################################################################################
# Land cover BAR PLOT ALONE, not related with TN and TP
############################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

land_cover_colors <- c(
  'Urban' = '#7f7f7f',                     # Gray for Urban
  'Cropland' = '#ff7f0e',                  # Orange for Cropland
  'Pastures and Open Nature' = '#bcbd22',  # Light green/yellow for Pastures
  'Forest' = '#2ca02c',                    # Green for Forest
  'Aquatic' = '#1f77b4'                    # Blue for Aquatic
)

#boxplot
data_long <- PC7 %>%
  select(Country, Aquatic_500, Cropland_500, Forest_500, Pastures.and.open.nature_500, Urban_500) %>%
  pivot_longer(cols = -Country, names_to = "Land_cover_Category", values_to = "Value")

# Map land cover types to your defined colors
data_long <- data_long %>%
  mutate(Land_cover_Category = recode(Land_cover_Category,
                                Aquatic_500 = 'Aquatic',
                                Cropland_500 = 'Cropland',
                                Forest_500 = 'Forest',
                                Pastures.and.open.nature_500 = 'Pastures and Open Nature',
                                Urban_500 = 'Urban'))


png("boxplot_land_cover_by_country.png", height = 7, width = 10, units = "in", res = 600)
ggplot(data_long, aes(x = Land_cover_Category, y = Value, fill = Land_cover_Category)) +
  geom_boxplot(outlier.size = 2, outlier.colour = "red") +  # Customize outliers
  labs(title = "Distribution of Land cover Types by Country",
       x = "Land cover Type",
       y = "Value") +
  facet_wrap(~ Country, ncol=1) +  # Create separate plots for each country
  scale_fill_manual(values = land_cover_colors) +  # Apply custom colors
  theme_minimal(base_size = 15) +  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) 
dev.off()

#barplot
# Reshape the data to long format and calculate means (or other summaries)
data_long <- PC7 %>%
  select(Country, Aquatic_500, Cropland_500, Forest_500, Pastures.and.open.nature_500, Urban_500) %>%
  pivot_longer(cols = -Country, names_to = "Land_cover_Category", values_to = "Value") %>%
  mutate(Land_cover_Category = recode(Land_cover_Category,
                                Aquatic_500 = 'Aquatic',
                                Cropland_500 = 'Cropland',
                                Forest_500 = 'Forest',
                                Pastures.and.open.nature_500 = 'Pastures and Open Nature',
                                Urban_500 = 'Urban')) %>%
  group_by(Country, Land_cover_Category) %>%
  summarize(Mean_Value = mean(Value, na.rm = TRUE), .groups = 'drop')  # Calculate mean values

data_long$Land_cover_Category <- factor(data_long$Land_cover_Category, levels=c('Urban', 'Cropland', 'Pastures and Open Nature', 'Forest','Aquatic'))
data_long$Country <- factor(data_long$Country, levels=c("Uruguay","Turkey","Spain","Switzerland","Belgium","Germany","UK","Denmark"))


# Create the bar plot
png("horizontal_stacked_barplot_land_cover.png", height = 7, width = 10, units = "in", res = 600)
ggplot(data_long, aes(x = Country, y = Mean_Value, fill = Land_cover_Category)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  labs(title = "Mean Distribution of Land cover Types by Country",
       x = "Country",
       y = "Mean Value",
       fill = "Land Use Type") +
  scale_fill_manual(values = land_cover_colors) +  # Apply custom colors
  coord_flip() +  # Flip coordinates for horizontal bars
  theme_minimal(base_size = 15) +  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold x-axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )
dev.off()


############################################################################################
#Livestock
############################################################################################
A_mean <- PC7 %>%
  group_by(Country) %>%
  summarize(Mean_Value = mean(Animals_cont, na.rm = TRUE), .groups = 'drop')  # Calculate mean values
A_mean$Country <- factor(A_mean$Country, levels=c("Uruguay","Turkey","Spain","Switzerland","Belgium","Germany","UK","Denmark"))


# Create the bar plot Livestock
png("Livestock_byCountry.png", height=7, width=4, units = "in", res = 600)
ggplot(A_mean, aes(x = Country, y = Mean_Value)) +
  geom_bar(stat = "identity") +  
  coord_flip()  +
  theme_minimal(base_size = 15) +  # Clean theme with larger font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) 
dev.off()

############################################################################################
# Climate
############################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

#boxplot
data_long <- PC7 %>%
  select(Country, bio1, bio5, bio7, bio12) %>%
  pivot_longer(cols = -Country, names_to = "Category", values_to = "Value")
data_long$Country <- factor(data_long$Country, levels=c("Uruguay","Turkey","Spain","Switzerland","Belgium","Germany","UK","Denmark"))

data_long$Category = as.factor(data_long$Category)
data_long$Category <- factor(data_long$Category, 
                             levels = c("bio1", "bio5", "bio7", "bio12"),
                             labels = c("Annual Mean Temperature", 
                                        "Maximum Temperature of Warmest Month", 
                                        "Temperature Annual Range", 
                                        "Annual Precipitation"))

levels(data_long$Country)[levels(data_long$Country) == "Turkey"] <- "Türkiye"


png("boxplot_climate_by_country.png", height = 12, width = 10, units = "in", res = 600)
ggplot(data_long, aes(x = Country, y = Value, fill = Country)) +
  geom_boxplot(outlier.size = 2, outlier.colour = "red") +  # Customize outliers
  facet_wrap(~ Category, ncol=1,scales = "free") +  # Create separate plots for each country
  theme_minimal(base_size = 15) +  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) 
dev.off()

############################################################################################
# Depth and Area Box PLOT ALONE, not related with TN and TP
############################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

#boxplot
data_long <- PC7 %>%
  select(Country, Depth, Area) %>%
  pivot_longer(cols = -Country, names_to = "Category", values_to = "Value")
data_long$Country <- factor(data_long$Country, levels=c("Uruguay","Turkey","Spain","Switzerland","Belgium","Germany","UK","Denmark"))

levels(data_long$Country)[levels(data_long$Country) == "Turkey"] <- "Türkiye"

#data_long[data_long$Value >20000,] --> 55792

data_long= data_long[data_long$Value <20000,]

png("boxplot_depth_area_by_country.png", height = 10, width = 10, units = "in", res = 600)
ggplot(data_long, aes(x = Country, y = Value, fill = Country)) +
  geom_boxplot(outlier.size = 2, outlier.colour = "red") +  # Customize outliers
  facet_wrap(~ Category, ncol=1,scales = "free") +  # Create separate plots for each country
  theme_minimal(base_size = 15) +  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) 
dev.off()







df$y_trans <- ifelse(df$y > 100, df$y - 80, df$y)

ggplot(df, aes(x, y_trans)) +
  geom_line() +
  facet_wrap(~ category) +
  scale_y_continuous(
    breaks = c(0, 20, 50, 120), 
    labels = c(0, 20, 50, "120+")
  ) +

############################################################################################
# TN and TP Box PLOT ALONE
############################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)

#boxplot
data_long <- PC7 %>%
  select(Country, TN, TP) %>%
  pivot_longer(cols = -Country, names_to = "Category", values_to = "Value")
data_long$Country <- factor(data_long$Country, levels=c("Uruguay","Turkey","Spain","Switzerland","Belgium","Germany","UK","Denmark"))

levels(data_long$Country)[levels(data_long$Country) == "Turkey"] <- "Türkiye"


png("boxplot_TN_TP_by_country.png", height = 7, width = 10, units = "in", res = 600)
ggplot(data_long, aes(x = Country, y = Value, fill = Country)) +
  geom_boxplot(outlier.size = 2, outlier.colour = "red") +  # Customize outliers
  facet_wrap(~ Category, ncol=1,scales = "free") +  # Create separate plots for each country
  theme_minimal(base_size = 15) +  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(y = "mg/L")
dev.off()




############################################################################################
# AREA DEPTH
############################################################################################
library(ggplot2)

plot1 <- ggplot(PC7, aes(x = Depth.t, y = Area.t)) +
  geom_point(color = "black", size = 2) +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # Regression line
  labs(
    x = "Pond Depth (cm)", 
    y = expression(Pond~Area~(m^2))
  ) +  # Create separate plots for each country
  theme_minimal(base_size = 15) +  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"))

############################################################################################
# DO DEPTH
############################################################################################
library(ggplot2)
library(patchwork) 

library(bestNormalize)
shapiro.test(PC7$O2_con)
O2_conn <- bestNormalize(PC7$O2_con)
plot(O2_conn, leg_loc = "bottomright")
PC7$O2_con.t <- predict(O2_conn, newdata = PC7$bio1)
hist(PC7$O2_con.t)
shapiro.test(PC7$O2_con.t)


plot2 <- ggplot(PC7, aes(x = Depth.t, y = O2_con)) +
  geom_point(color = "black", size = 2)  +
  labs(
    x = "Pond Depth (cm)", 
    y = "Dissolved oxygen concentration  (mg/L)"
  ) +  # Create separate plots for each country
  theme_minimal(base_size = 15) +  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Bold and centered title
    axis.text.x = element_text(face = "bold"),  # Bold axis text
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"))

# Combine plots into two columns
combined_plot <- plot1 + plot2 + plot_layout(ncol = 2)

# Save to a PNG file
ggsave("Scatter_plot_Depth_Area_DO.png", plot = combined_plot, width = 10, height = 5, dpi = 300)


