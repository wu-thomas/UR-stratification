################################################################
#########   load libraries
################################################################
#### Libraries ####

rm(list = ls())

options(java.parameters = "-Xmx10g")
# gc()
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

if(!dir.exists(paths = paste0('BART_jitter'))){
  dir.create(path = paste0('BART_jitter'))
}

if(!dir.exists(paths = paste0('BART_jitter/samp_location'))){
  dir.create(path = paste0('BART_jitter/samp_location'))
}

if(!dir.exists(paths = paste0('BART_jitter/natl_pred'))){
  dir.create(path = paste0('BART_jitter/natl_pred'))
}

################################################################
######### load polygons
################################################################


setwd(paste(data_dir,country,sep=''))

gadm.abbrev = "MWI"

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
setwd(paste(data_dir,country,sep=''))

worldpop <- raster(paste0('worldpop/',country.abbrev,
                          '_ppp_',frame.year,
                          '_1km_Aggregated_UNadj.tif',sep=''))


################################################################
#########  Load cluster data
################################################################

setwd(cluster_dir)
load(paste0(country,'_cluster_dat.rda'),
     envir = .GlobalEnv)

# clusters

cluster_list<-mod.dat[!duplicated(mod.dat[c('cluster','survey',
                                            'LONGNUM','LATNUM')]),]
cluster_list<-cluster_list[cluster_list$survey==2015|
                             cluster_list$survey==2010,]

################################################################
######### load national grid 
################################################################

setwd(paste0(data_dir,country,sep=''))
load("prepared_dat/natl_grid.rda")

################################################################
#########  load potential locations
################################################################

setwd(paste0(data_dir,country,'/prepared_dat/',sep=''))
load(file='cluster_jitter_frame_nolikoma.rda')
#urban_crc <- readRDS(file='crc_dat.rds')
#urban_crc$light <- urban_crc$noaa


clus_loc_frame <- do.call(rbind, clust_nolikoma)

################################################################
#########  function to draw random samples from the master list
################################################################

# assign location id
clus_loc_frame <- clus_loc_frame %>% group_by(id_svy) %>% mutate(loc_id = sequence(n()))

# test_samp <- clus_loc_frame %>% group_by(id_svy) %>% sample_n(1,weight=norm_w)

################################################################
#########  iteration
################################################################

cov_list<-c('pop_den','light','ndvi','admin2.char','elev','prec','tavg','access')


# prepare samples for first iteration
clus_loc_frame$new_w <- clus_loc_frame$norm_w 


set_bart_machine_num_cores(12)

# prepare national prediction grid

natl_comp<-urb_dat[complete.cases(urb_dat),]



for (iteration in 1:2000){
  
  print(iteration)
  ptm <- proc.time()
  
  
  ite_samp <- clus_loc_frame %>% group_by(id_svy) %>% sample_n(1,weight=new_w)
  
  setwd(paste0(res_dir,country,'/BART_jitter/samp_location',sep=''))
  
  save(ite_samp,file=paste0('sampled_loc_ite_',iteration,'.rda',sep=''))
  
  xlcFreeMemory()
  
  bart_temp<-build_bart_machine(X=as.data.frame(ite_samp[,cov_list]),
                                y=as.factor(ite_samp$urban),
                                num_trees = 200,
                                num_burn_in = 300,
                                num_iterations_after_burn_in = 1,
                                use_missing_data=TRUE)
  
  xlcFreeMemory()
  
  
  clus_loc_frame$urb_p <- predict(bart_temp,as.data.frame(clus_loc_frame[,cov_list]))
  clus_loc_frame$new_w <- ifelse(clus_loc_frame$urban == 'urban',
                                 clus_loc_frame$norm_w*clus_loc_frame$urb_p, 
                                 clus_loc_frame$norm_w*(1-clus_loc_frame$urb_p) )
  
  clus_loc_frame <- clus_loc_frame %>% group_by(id_svy)  %>% 
    mutate_at(vars(new_w), funs(./sum(.))) 
  
  
  xlcFreeMemory()
  
  
  
  if(FALSE) {
  natl_pred_vec <- vector()
  
  for ( grid_group in 1:num_groups){
    xlcFreeMemory()
    
    # print(grid_group)
    natl_sec <- natl_comp_group[[grid_group]]
    natl_sec_pred <- predict(bart_temp, natl_sec[,cov_list])
    
    natl_pred_vec <- c(natl_pred_vec,natl_sec_pred)
    
  }
  }
  
  
  natl_pred_vec <-  predict(bart_temp, natl_comp[,cov_list])
    
  setwd(paste0(res_dir,country,'/BART_jitter/natl_pred',sep=''))
  save(natl_pred_vec,file=paste0('natl_pred_ite_',iteration,'.rda',sep=''))
  
  
  proc.time() - ptm
  
  
  xlcFreeMemory()
 
  
}
