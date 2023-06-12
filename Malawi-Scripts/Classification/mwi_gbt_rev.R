################################################################
#########   load libraries
################################################################
#### Libraries ####

rm(list = ls())

#options(java.parameters = "-Xmx10g")
#gc()
library(SUMMER)
#help(package = "SUMMER", help_type = "html")
#utils::browseVignettes(package = "SUMMER")
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(INLA)
library(survey)
library(ggplot2)
library(maptools)
library(gridExtra)
library(mgcv)
library(caret)
library(geosphere)
library(rgeos)
library(sqldf)
library(raster)
library(parallel)
library(tictoc)
library(bartMachine)
#library(XLConnect)

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
frame.year<-2008
country.abbrev<-'mwi' #for population rasters
gadm.abbrev = "MWI"

data_dir<-paste0(home_dir,'/Data/')
res_dir<-paste0(home_dir,'/Results/')
pop_dir<-paste0(main_dir,'/Population/')
cluster_dir<-paste0(main_dir,'countryDataFolders/',country)

# create directories to store results
setwd(paste0(res_dir,country,sep=''))

if(!dir.exists(paths = paste0('GBT'))){
  dir.create(path = paste0('GBT'))
}



################################################################
######### load polygons
################################################################

setwd(paste(data_dir,country,sep=''))


poly.path <- paste0("shapeFiles_gadm")


poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_")
poly.layer.adm1 <- paste('new', gadm.abbrev,
                         '1', sep = "_")
poly.layer.adm2 <- paste('new', gadm.abbrev,
                         '2', sep = "_")

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

load(paste0('shapeFiles_gadm/', country, '_Amat.rda'))
load(paste0('shapeFiles_gadm/', country, '_Amat_Names.rda'))


################################################################
#########   load worldpop
################################################################

setwd(paste(data_dir,country,sep=''))

worldpop <- raster(paste0('worldpop/',country.abbrev,'_ppp_',frame.year,'_','1km_Aggregated_UNadj.tif',sep=''))

################################################################

################################################################
######### load national grid 
################################################################

setwd(paste0(data_dir,country,sep=''))
load("prepared_dat/natl_grid.rda")


complete_id <- complete.cases(urb_dat)

################################################################
#########  load cluster training data
################################################################

setwd(paste0(data_dir,country,'/prepared_dat/',sep=''))
load(file='crc_dat_nolikoma.rda')
load(file='uncrc_dat_nolikoma.rda')


################################################################
#########  GBT using uncorrected locations
################################################################


## set up mean structure
cov_list<-c('pop_den','light','ndvi','admin2.char','elev','prec','tavg','access')
mean_formula<-as.formula(paste('urban~pop_den+elev+ndvi+tavg+prec+access+
                               light+as.factor(admin2.char)'))

train_dat <- uncrc_dat_nolikoma

gbtGrid <- expand.grid(interaction.depth=c(3,5,7,9), 
                       n.trees = (20:40)*50,
                       shrinkage=c(0.01,0.005,0.001),
                       n.minobsinnode=10)


set.seed(2021)
fit_natl<-caret::train(mean_formula,
                       data = train_dat,
                       method = "gbm",
                       na.action = na.pass,
                       trControl = trainControl(method="cv", 
                                                number=10),
                       verbose = 0,
                       tuneGrid=gbtGrid)


summary(fit_natl)
fit_natl$bestTune
## prediction
complete_urb_dat<-urb_dat[complete_id,]

complete_urb_dat$admin2.char.rep <-complete_urb_dat$admin2.char
complete_urb_dat[complete_urb_dat$admin2.char=='admin2_10',]$admin2.char <-NA

pred_prob<-predict(fit_natl, complete_urb_dat,type = "prob",na.action = na.pass)




complete_urb_dat$pred_p <- pred_prob$urban

complete_urb_dat$admin2.char <- complete_urb_dat$admin2.char.rep


natl_pred_vec <- complete_urb_dat$pred_p

setwd(paste0(res_dir,country,'/GBT',sep=''))
save(natl_pred_vec,file=paste0('gbt_pred_p_uncrc.rda'))
save(complete_urb_dat,file=paste0('gbt_pred_p_frame_uncrc.rda'))

####
#complete_urb_dat$pred_p <- natl_pred_vec

complete_urb_dat$urban_ind <- as.numeric(complete_urb_dat$pred_p>0.5)

sum(complete_urb_dat$urban_ind*complete_urb_dat$pop_den)/sum(complete_urb_dat$pop_den)

###

tmp.likoma <- complete_urb_dat[complete_urb_dat$admin2.char=='admin2_10',]







################################################################
#########  GBT using corrected locations
################################################################


## set up mean structure
cov_list<-c('pop_den','light','ndvi','admin2.char','elev','prec','tavg','access')
mean_formula<-as.formula(paste('urban~pop_den+elev+ndvi+tavg+prec+access+
                               light+as.factor(admin2.char)'))

train_dat <- crc_dat_nolikoma

gbtGrid <- expand.grid(interaction.depth=c(3,5,7,9), 
                       n.trees = (20:40)*50,
                       shrinkage=c(0.01,0.005,0.001),
                       n.minobsinnode=10)


set.seed(2021)
fit_natl<-caret::train(mean_formula,
                       data = train_dat,
                       method = "gbm",
                       na.action = na.pass,
                       trControl = trainControl(method="cv", 
                                                number=10),
                       verbose = 0,
                       tuneGrid=gbtGrid)


summary(fit_natl)
fit_natl$bestTune
## prediction
complete_urb_dat<-urb_dat[complete_id,]

complete_urb_dat$admin2.char.rep <-complete_urb_dat$admin2.char
complete_urb_dat[complete_urb_dat$admin2.char=='admin2_10',]$admin2.char <-NA

pred_prob<-predict(fit_natl, complete_urb_dat,type = "prob",na.action = na.pass)




complete_urb_dat$pred_p <- pred_prob$urban

complete_urb_dat$admin2.char <- complete_urb_dat$admin2.char.rep


natl_pred_vec <- complete_urb_dat$pred_p

setwd(paste0(res_dir,country,'/GBT',sep=''))
save(natl_pred_vec,file=paste0('gbt_pred_p_crc.rda'))
save(complete_urb_dat,file=paste0('gbt_pred_p_frame_crc.rda'))

