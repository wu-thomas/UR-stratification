################################################################
#########   load libraries
################################################################
#### Libraries ####

rm(list = ls())

options(java.parameters = "-Xmx10g")
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
library(XLConnect)

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

if(!dir.exists(paths = paste0('BART_corrected'))){
  dir.create(path = paste0('BART_corrected'))
}


if(!dir.exists(paths = paste0('BART_uncorrected'))){
  dir.create(path = paste0('BART_uncorrected'))
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

################################################################
#########   load worldpop
################################################################

worldpop <- raster('worldpop/mwi_ppp_2008_1km_Aggregated_UNadj.tif')


################################################################
######### load national grid 
################################################################

setwd(paste0(data_dir,country,sep=''))
load("prepared_dat/natl_grid.rda")

################################################################
#########  load cluster training data
################################################################

setwd(paste0(data_dir,country,'/prepared_dat/',sep=''))
load(file='crc_dat_nolikoma.rda')
load(file='uncrc_dat_nolikoma.rda')




################################################################
#########  classification using corrected locations 
################################################################

cov_list<-c('pop_den','light','ndvi','admin2.char','elev','prec','tavg','access')

train_dat <- crc_dat_nolikoma

set_bart_machine_num_cores(12)

# prepare national prediction grid

num_groups <- 20
natl_comp<-urb_dat[complete.cases(urb_dat),]
natl_comp_group<-natl_comp %>% 
  group_by((row_number()-1) %/% (n()/num_groups)) %>%
  nest %>% pull(data)


ptm <- proc.time()




xlcFreeMemory()

bart_temp<-build_bart_machine(X=as.data.frame(train_dat[,cov_list]),
                              y=as.factor(train_dat$urban),
                              num_trees = 200,
                              num_burn_in = 300,
                              num_iterations_after_burn_in = 1000,
                              use_missing_data=TRUE)

xlcFreeMemory()



#temp_pred <- predict(bart_temp,train_dat[,cov_list])
#tmp_train <- train_dat
#tmp_train$pred <- temp_pred

# get full posterior 

natl_pred_post <- data.frame()

for ( grid_group in 1:num_groups){
  xlcFreeMemory()
  
  # print(grid_group)
  natl_sec <- natl_comp_group[[grid_group]]
  natl_sec_pred <- bart_machine_get_posterior(bart_temp, natl_sec[,cov_list])
  natl_sec_pred_draws <- natl_sec_pred$y_hat_posterior_samples
    
  natl_pred_post <- rbind(natl_pred_post,natl_sec_pred_draws)
  
}


setwd(paste0(res_dir,country,'/BART_corrected',sep=''))
save(natl_pred_post,file=paste0('nolikoma_natl_pred_crc_draws_1000.rda'))



if(FALSE){
  
  ### getting posterior mean
natl_pred_vec <- vector()

for ( grid_group in 1:num_groups){
  xlcFreeMemory()
  
  # print(grid_group)
  natl_sec <- natl_comp_group[[grid_group]]
  natl_sec_pred <- predict(bart_temp, natl_sec[,cov_list])
  
  natl_pred_vec <- c(natl_pred_vec,natl_sec_pred)
  
}



proc.time() - ptm


xlcFreeMemory()
}






################################################################
#########  classification using uncorrected locations 
################################################################

cov_list<-c('pop_den','light','ndvi','admin2.char','elev','prec','tavg','access')

train_dat <- uncrc_dat_nolikoma

set_bart_machine_num_cores(12)

# prepare national prediction grid

num_groups <- 20
natl_comp<-urb_dat[complete.cases(urb_dat),]
natl_comp_group<-natl_comp %>% 
  group_by((row_number()-1) %/% (n()/num_groups)) %>%
  nest %>% pull(data)


ptm <- proc.time()




xlcFreeMemory()

bart_temp<-build_bart_machine(X=as.data.frame(train_dat[,cov_list]),
                              y=as.factor(train_dat$urban),
                              num_trees = 200,
                              num_burn_in = 300,
                              num_iterations_after_burn_in = 1000,
                              use_missing_data=TRUE)

xlcFreeMemory()


# get full posterior 

natl_pred_post <- data.frame()

for ( grid_group in 1:num_groups){
  xlcFreeMemory()
  print(grid_group)
  # print(grid_group)
  natl_sec <- natl_comp_group[[grid_group]]
  natl_sec_pred <- bart_machine_get_posterior(bart_temp, natl_sec[,cov_list])
  natl_sec_pred_draws <- natl_sec_pred$y_hat_posterior_samples
  
  natl_pred_post <- rbind(natl_pred_post,natl_sec_pred_draws)
  
}


setwd(paste0(res_dir,country,'/BART_uncorrected',sep=''))
save(natl_pred_post,file=paste0('nolikoma_natl_pred_uncrc_draws_1000.rda'))



if(FALSE){

  # point estimate
natl_pred_vec <- vector()

for ( grid_group in 1:num_groups){
  xlcFreeMemory()
  
  # print(grid_group)
  natl_sec <- natl_comp_group[[grid_group]]
  natl_sec_pred <- predict(bart_temp, natl_sec[,cov_list])
  
  natl_pred_vec <- c(natl_pred_vec,natl_sec_pred)
  
}



setwd(paste0(res_dir,country,'/BART_uncorrected',sep=''))
save(natl_pred_vec,file=paste0('nolikoma_natl_pred_uncrc_500.rda'))


proc.time() - ptm


xlcFreeMemory()
}





