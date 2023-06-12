# SAE : HIV Example - Malawi (women aged 15-29)
# Original Author: Taylor Okonek

################################################################
#########   load libraries
################################################################
#### Libraries ####

rm(list = ls())
library(rgdal)
library(SUMMER)
library(INLA)
library(haven)
library(rgeos)
library(ggplot2)
library(maptools)
library(geoR)
library(gstat)
library(raster)
library(verification)
library(scoringRules)
library(dplyr)


library(tidyverse)
library(readstata13)
library(survey)
library(spdep)
library(gridExtra)
library(fastmatch)
library(ggridges)
library(GGally)
library(knitr)
library(doMC)
library(viridis)
library(fields)
library(ggpubr)

################################################################
#########   set parameters
################################################################

## set directory
country<-'Malawi'

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-4)], collapse = "/")

# path to AfricaAdmin2Estimates folder 
main_dir<-'E:/Dropbox/AfricaAdmin2Estimates/Data/'
frame.year<-  sampleFrame <- 2008

country.abbrev<-'mwi' #for population rasters
gadm.abbrev = "MWI"

data_dir<-paste0(home_dir,'/Data/')
res_dir<-paste0(home_dir,'/Results/')
pop_dir<-paste0(main_dir,'/Population/')
cluster_dir<-paste0(main_dir,'countryDataFolders/',country)
code_dir <- paste0(home_dir,'/Scripts/')

# create directories to store prepared data
setwd(paste0(res_dir,country,sep=''))

if(!dir.exists(paths = paste0('HIV'))){
  dir.create(path = paste0('HIV'))
}

setwd(paste(data_dir,country,sep=''))

if(!dir.exists(paths = paste0('HIV'))){
  dir.create(path = paste0('HIV'))
}


################################################################
#########   load helper functions
################################################################

#setwd(paste0(code_dir,'/',country,'/Stunting'))
#source('helper_functions_stunt.R')


################################################################
#########   load polygons
################################################################

#### Load polygon files ####
setwd(paste(data_dir,country,sep=''))

poly.path <- paste0("shapeFiles_gadm")

poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_")
poly.layer.adm1 <- paste('new', gadm.abbrev,
                         '1', sep = "_")
poly.layer.adm2 <- paste('new', gadm.abbrev,
                         '2', sep = "_")
poly.layer.adm3 <- paste('gadm36', gadm.abbrev,
                         '3', sep = "_")

poly.adm0 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm0)) 
# use encoding to read special characters
poly.adm1 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm1))

if(sum(grepl(paste('new', gadm.abbrev,
                   '2', sep = "_"), list.files(poly.path))) != 0){
  poly.adm2 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                       layer = as.character(poly.layer.adm2))}

if(exists("poly.adm2")){
  proj4string(poly.adm0) <- proj4string(poly.adm1)  <- proj4string(poly.adm2)
}else{
  proj4string(poly.adm0) <- proj4string(poly.adm1)
} 

poly.adm3 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm3)) 

load(paste0('shapeFiles_gadm/', country, '_Amat.rda'))
load(paste0('shapeFiles_gadm/', country, '_Amat_Names.rda'))


################################################################
#########   prepare admin3 
################################################################

if(FALSE){
setwd(paste(data_dir,country,sep=''))

poly.adm3 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm3)) 

plot(poly.adm3)
  
admin3.mat <- poly2nb(SpatialPolygons(poly.adm3@polygons))
admin3.mat <- nb2mat(admin3.mat, zero.policy = TRUE)
colnames(admin3.mat) <- rownames(admin3.mat) <- paste0("admin3_", 1:dim(admin3.mat)[1])
admin3.names <- data.frame(GADM = poly.adm3@data$NAME_2,
                           Internal = rownames(admin3.mat))


save(admin1.mat, admin2.mat, admin3.mat,file = paste0(poly.path,'/', country, '_Amat.rda'))
save(admin1.names, admin2.names, admin3.names, file = paste0(poly.path, '/', country, '_Amat_Names.rda'))
}
################################################################
#########  load urban fraction for female 15-29 population 
################################################################

## Extract populations and urban/rural

setwd(paste0(res_dir,country,'/UR/Fractions/',sep=''))
load('uncrc_adm2_thresh_f15_29_natl_frac.rda')
load('uncrc_adm2_thresh_f15_29_adm1_frac.rda')
load('uncrc_adm2_thresh_f15_29_adm2_frac.rda')

f15_29.adm2.frac[f15_29.adm2.frac$admin.name=='Likoma',]$Frac<- 1352/10414


################################################################
######### prepare spatial
################################################################


Admin1Map <- poly.adm1
Admin2Map <- poly.adm2
Admin3Map <- poly.adm3


# Read DHS coordinates
setwd(paste(data_dir,country,sep=''))

GPS.path <- paste0("dhsFlat/MWGE7AFL")

corList = readOGR(dsn = GPS.path,
                  layer = "MWGE7AFL")

idxKeep = corList$LATNUM != 0
corList = corList[idxKeep,]


# Get graph
nameVec_adm2 = Admin2Map$NAME_1
for(i in 1:length(nameVec_adm2)){
  nameVec_adm2[i] = paste(nameVec_adm2[i], ":", Admin2Map$NAME_2[i], sep = "")
}

Admin1Graph = getAmat(Admin1Map, Admin1Map$NAME_1)
Admin2Graph = getAmat(Admin2Map, names = nameVec_adm2)

# admin-3 graph
nameVec_adm3 = Admin3Map$NAME_1
for(i in 1:length(nameVec_adm3)){
  nameVec_adm3[i] = paste(nameVec_adm3[i], ":", Admin3Map$NAME_2[i], sep = "")
}

Admin3Graph = getAmat(Admin3Map, names = nameVec_adm3)

# Make INLA graph object
Admin3GraphINLA = inla.read.graph(Admin3Graph)

################################################################
######### helper functions
################################################################

# useful functions
make_region <- function(x) {
  str_split(x," ") %>% unlist %>% nth(1)
}
make_region <- Vectorize(make_region)
make_urban <- function(x) {
  str_split(x," ") %>% unlist %>% nth(3)
}
make_urban <- Vectorize(make_urban)
split_sex <- function(x) {
  x %>% str_sub(7,8)
}


################################################################
######### load HIV data
################################################################

setwd(paste0(data_dir,country,'/dhsStata/',sep=''))

hiv_df <- readstata13::read.dta13("MWAR7ADT/MWAR7AFL.DTA", 
                                  generate.factors = TRUE)
#hiv_df2 <- read_dta("MWAR7ADT/MWAR7AFL.DTA")

# changing variable names so they're more intuitive
hiv_df <- hiv_df %>% mutate(clustid = hivclust,
                            id = hivnumb) %>%
  dplyr::select(-hivclust, -hivnumb)

################################################################
######### prepare household
################################################################

setwd(paste0(data_dir,country,'/dhsStata/',sep=''))

# read in household recode data - we need this for sex apparently
households_df <- readstata13::read.dta13("MWHR7HDT/MWHR7HFL.DTA",
                                         generate.factors = TRUE)

#households_df <- read_dta("MWHR7HDT/MWHR7HFL.DTA")

sex_columns <- which((colnames(households_df) %>% str_detect("hv104")) == TRUE)
households_temp <- households_df %>% dplyr::select(hv001, hv002, hv003, hv022) 

households_temp2 <- households_df[,c(3,4,sex_columns)] %>% 
  gather(hivline, sex, hv104_01:hv104_20) 

#mutate(hivline = sapply(hivline, function(x) {x %>% str_split("_") %>% unlist %>% nth(2)})) %>%
split_sex <- Vectorize(split_sex)
households_temp2$hivline <- split_sex(households_temp2$hivline)
households_temp2$hivline <- households_temp2$hivline %>% as.numeric() 

households_temp <- right_join(households_temp, households_temp2) 
households_temp <- households_temp %>% dplyr::select(-hv003)

colnames(households_temp) <- c("clustid","id","strata","hivline","sex")



################################################################
######### prepare individual
################################################################

setwd(paste0(data_dir,country,'/dhsStata/',sep=''))

#indiv_df <- read_dta("MWIR7HDT/MWIR7HFL.DTA")
indiv_df <- readstata13::read.dta13("MWIR7HDT/MWIR7HFL.DTA",
                                    generate.factors = TRUE)


# v022, v023 contain admin 2 stuff, can merge in by v001 (clustid) 
# if wealth index is something we want stratified by urban and rural, we should use v190a, v191a instead!!
indiv_temp <- indiv_df %>% 
  dplyr::select(v001, v002, v003, v012, v106, v133, v150, v190, v191, v151) %>% 
  distinct()
colnames(indiv_temp) <- c("clustid","id","hivline","age",
                          "education","education_year","head_of_household","wealth_index",
                          "wealth_index_factor_score","head_of_household_sex")



setwd(paste0(data_dir,country,'/HIV',sep=''))
#load(file = 'HIV_mwi.RData')

################################################################
######### prepare HIV data
################################################################

# remove some rows we know we can remove so that the merge is easier
households_temp <- households_temp %>% filter(hivline %in% (hiv_df$hivline %>% unique())) 


# Merge strata and admin 2 levels into hiv_df
hiv_df <- left_join(hiv_df, households_temp, by = c('clustid'='clustid', 'id'='id', 'hivline'='hivline'))
hiv_df <- left_join(hiv_df, indiv_temp, by = c('clustid'='clustid', 'id'='id', 'hivline'='hivline'))

# create hiv indicator variable
hiv_df <- hiv_df %>% 
  mutate(hiv_ind = ifelse(hiv03 == "hiv  positive" | hiv03 == "hiv2 positive hiv1 & hiv2 positive", 1,
                          ifelse(hiv03 == "hiv negative", 0, NA)))

hiv_df$clusterIdx <- hiv_df$clustid



################################################################
######### Add latitude longitude, Admin1, Admin2
################################################################

smallGeo = data.frame(clusterIdx = corList$DHSCLUST, urban = corList$URBAN_RURA,
                      lon = corList$LONGNUM, lat = corList$LATNUM,
                      admin1 = corList$ADM1NAME)

# add cluster xy information
hiv_df = merge(hiv_df, smallGeo,by= "clusterIdx")

# Use coordinates to place in admin2
clusters = unique(hiv_df[,c("clusterIdx", "lon.x", "lat.x")])
gridP = data.frame(Longitude = as.vector(clusters$lon),
                   Latitude = as.vector(clusters$lat))
coordinates(gridP) = ~ Longitude + Latitude
proj4string(gridP) = proj4string(Admin2Map)


admin1 = over(gridP, Admin1Map)
admin2 = over(gridP, Admin2Map)




# Replace NAs with nearest area
#proj4string(gridP) = proj4string(Admin2Map)

idx = which(is.na(admin2$NAME_2))
n = length(idx)
nearestArea = list()
distToNearest = c()

for(i in 1:n){
  gDists = gDistance(gridP[idx[i],], Admin1Map, byid = TRUE)
  tmpArea = Admin1Map[which.min(gDists),]
  distToNearest[i] = min(gDists)
  admin1$NAME_1[idx[i]] = tmpArea$NAME_1
}

for(i in 1:n){
  gDists = gDistance(gridP[idx[i],], Admin2Map, byid = TRUE)
  tmpArea = Admin2Map[which.min(gDists),]
  distToNearest[i] = min(gDists)
  admin2$NAME_2[idx[i]] = tmpArea$NAME_2
  admin2$NAME_1[idx[i]] = tmpArea$NAME_1
}

# Replace NAs with nearest area for admin-3
proj4string(gridP) = proj4string(Admin3Map)
admin3 = over(gridP, Admin3Map)

idx = which(is.na(admin3$NAME_2))
n = length(idx)
nearestArea = list()
distToNearest = c()

for(i in 1:n){
  gDists = gDistance(gridP[idx[i],], Admin3Map, byid = TRUE)
  tmpArea = Admin3Map[which.min(gDists),]
  distToNearest[i] = min(gDists)
  admin3$NAME_2[idx[i]] = tmpArea$NAME_2
  admin3$NAME_1[idx[i]] = tmpArea$NAME_1
}


# Fix names
admin2$NAME_2_fix <- NA
for(i in 1:dim(admin2)[1]){
  admin2$NAME_2_fix[i] = paste(admin2$NAME_1[i], ":", admin2$NAME_2[i], sep = "")
}

admin3$NAME_2_fix <- NA
for(i in 1:dim(admin3)[1]){
  admin3$NAME_2_fix[i] = paste(admin3$NAME_1[i], ":", admin3$NAME_2[i], sep = "")
}


clusters = list(clusterIdx = clusters$clusterIdx, 
                admin2 = admin2$NAME_2_fix,
                admin1.char = admin1$NAME_1,
                admin2.char = admin2$NAME_2,
                admin3 = admin3$NAME_2_fix,
                admin3.char = admin3$NAME_2)

# Add geographic information to dataset
hiv_df = merge(hiv_df, clusters, by = "clusterIdx")

hiv_df$admin2Idx = 0
for(i in 1:length(hiv_df$admin2Idx)){
  #print(i)
  hiv_df$admin2Idx[i] = which(nameVec_adm2 == hiv_df$admin2[i])
}

hiv_df$admin3Idx = 0
for(i in 1:length(hiv_df$admin3Idx)){
  #print(i)
  hiv_df$admin3Idx[i] = which(nameVec_adm3 == hiv_df$admin3[i])
}


# Make factor
hiv_df$admin2Fac = factor(hiv_df$admin2, levels = nameVec_adm2)
hiv_df$admin1Fac = as.factor(hiv_df$admin1.char)
hiv_df$admin1.char= as.factor(hiv_df$admin1.char)
hiv_df$admin2.char= as.factor(hiv_df$admin2.char)

hiv_df$admin3.char= as.factor(hiv_df$admin3.char)
hiv_df$admin3Fac = factor(hiv_df$admin3, levels = nameVec_adm3)


setwd(paste0(data_dir,country,'/HIV',sep=''))
#save.image(file = 'HIV_mwi.RData')
#save(hiv_df, file = "hiv_df_adm3.rda")
#load(file = "hiv_df.rda")
#hiv_df <- hiv_df[,-which(names(hiv_df) %in%c('admin2Fac','admin1Fac','admin1.char','admin2.char','admin2Idx',
#                                             'admin2','admin1'))]



################################################################
######### prepare anlaysis data set
################################################################

# remove missing values of responseVar so function will run -- DOUBLE CHECK THIS
hiv_df_no_na <- hiv_df %>% filter(!is.na(hiv_ind))
# only include young women
hiv_df_subgroup <- hiv_df_no_na %>% filter(sex == "female" & age < 30) 


hiv_df_subgroup <- hiv_df_subgroup %>% mutate(strata = as.character(strata))



