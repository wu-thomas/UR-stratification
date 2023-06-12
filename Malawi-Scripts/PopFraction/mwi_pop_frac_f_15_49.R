################################################################
#########   load libraries
################################################################
#### Libraries ####

rm(list = ls())

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
library(stringdist)
library(data.table)
library(nleqslv)

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

#main_dir<-'E:/Dropbox/AfricaAdmin2Estimates/Data/'
svy.year <- 2015
frame.year<-2008
country.abbrev<-'mwi' #for population rasters
gadm.abbrev = "MWI"


data_dir<-paste0(home_dir,'/Data/')
res_dir<-paste0(home_dir,'/Results/')
#pop_dir<-paste0(main_dir,'/Population/')
#cluster_dir<-paste0(main_dir,'countryDataFolders/',country)
script_dir<-paste0(home_dir,'/Scripts/')

### load functions
setwd(paste0(script_dir,country,'/PopFraction/',sep=''))
source('mwi_pop_frac_functions.R')

# national urban fraction for threshold

natl_urban_frac <- 2003309/13077160
# cutoff <- natl_urban_frac


# create directories to store results
setwd(paste0(res_dir,country,sep=''))

if(!dir.exists(paths = paste0('UR'))){
  dir.create(path = paste0('UR'))
}


if(!dir.exists(paths = paste0('UR/Fractions'))){
  dir.create(path = paste0('UR/Fractions'))
}

if(!dir.exists(paths = paste0('UR/Fractions/Samples'))){
  dir.create(path = paste0('UR/Fractions/Samples'))
}

if(!dir.exists(paths = paste0('UR/Fractions/Uncertainty'))){
  dir.create(path = paste0('UR/Fractions/Uncertainty'))
}

if(!dir.exists(paths = paste0('UR/Classification'))){
  dir.create(path = paste0('UR/Classification'))
}

if(!dir.exists(paths = paste0('UR/Classification/Thresh'))){
  dir.create(path = paste0('UR/Classification/Thresh'))
}

################################################################
######### load polygons
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

poly.adm3 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm3)) 

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


u5.pop <- raster(paste0('worldpop/',country.abbrev,'_u5_',svy.year,'_1km.tif',sep=''))


f.15.49.pop <- raster(paste0('worldpop/',country.abbrev,'_f15_49_',svy.year,'_1km.tif',sep=''))





################################################################
######### load national grid 
################################################################

setwd(paste0(data_dir,country,sep=''))
load("prepared_dat/natl_grid.rda")

urb_dat$complete <- complete.cases(urb_dat) 
urb_dat$pixel_id <- c(1:dim(urb_dat)[1])

complete_id <- complete.cases(urb_dat) 
complete_urb_dat <- urb_dat[complete_id,]

complete_urb_dat$f_15_49_pop <- raster::extract(f.15.49.pop,complete_urb_dat[c('x','y')])

setwd(paste0(data_dir,country,sep=''))

load('admin2_urb_ref_tab.rda')


#tmp_urb <- complete_urb_dat
#tmp_urb$tmp_pred <- natl_pred_post$V1

complete.natl.dat <- urb_dat[urb_dat$complete,]
pixel.grid <- complete.natl.dat[,c('x','y')]

##################################################################
######### uncorrected BART
##################################################################

## load draws from BART unjitter using uncorrected locations
method= 'BART_uncrc'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
load(file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm2 ,
                            admin.level = 'admin2',
                            admin.name=admin2.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm3 ,
                            admin.level = 'admin3',
                            admin.name=admin3.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))




##################################################################
######### corrected BART
##################################################################

## load draws from BART unjitter using uncorrected locations
method= 'BART_crc'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
load(file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm2 ,
                            admin.level = 'admin2',
                            admin.name=admin2.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm3 ,
                            admin.level = 'admin3',
                            admin.name=admin3.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))



##################################################################
######### BART MCMC
##################################################################

## load draws from BART unjitter using uncorrected locations
method= 'BART_MCMC'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
load(file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm2 ,
                            admin.level = 'admin2',
                            admin.name=admin2.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm3 ,
                            admin.level = 'admin3',
                            admin.name=admin3.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))




##################################################################
######### logistic regression with uncorrected locations
##################################################################


## load draws from BART unjitter using uncorrected locations
method= 'logistic_uncrc'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
load(file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm2 ,
                            admin.level = 'admin2',
                            admin.name=admin2.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm3 ,
                            admin.level = 'admin3',
                            admin.name=admin3.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))





##################################################################
######### logistic regression with corrected locations
##################################################################



## load draws from BART unjitter using uncorrected locations
method= 'logistic_crc'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
load(file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm2 ,
                            admin.level = 'admin2',
                            admin.name=admin2.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm3 ,
                            admin.level = 'admin3',
                            admin.name=admin3.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))




##################################################################
######### GBT with uncorrected locations
##################################################################


## load draws from BART unjitter using uncorrected locations
method= 'GBT_uncrc'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
load(file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm2 ,
                            admin.level = 'admin2',
                            admin.name=admin2.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm3 ,
                            admin.level = 'admin3',
                            admin.name=admin3.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))






##################################################################
######### GBT with corrected locations
##################################################################

## load draws from BART unjitter using uncorrected locations
method= 'GBT_crc'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
load(file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm2 ,
                            admin.level = 'admin2',
                            admin.name=admin2.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                            pop.ras = f.15.49.pop,
                            urb.prob.frame = calibrated.prob.samp,
                            admin.poly = poly.adm3 ,
                            admin.level = 'admin3',
                            admin.name=admin3.names)

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))

save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))




