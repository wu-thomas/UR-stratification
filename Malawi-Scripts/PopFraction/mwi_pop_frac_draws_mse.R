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

if(FALSE){
  setwd(paste(data_dir,country,sep=''))
  
### load 1-5 population
  pop.f0 <-  raster(paste0('worldpop/',country.abbrev,'_f_0_',svy.year,'.tif',sep=''))
  pop.m0 <-  raster(paste0('worldpop/',country.abbrev,'_m_0_',svy.year,'.tif',sep=''))
  
pop.f1 <-  raster(paste0('worldpop/',country.abbrev,'_f_1_',svy.year,'.tif',sep=''))
pop.m1 <-  raster(paste0('worldpop/',country.abbrev,'_m_1_',svy.year,'.tif',sep=''))

pop.u5 <- pop.f0+pop.m0+pop.f1+pop.m1

#plot(pop.u5)
sum(values(pop.m1),na.rm=T)
writeRaster(pop.u5, file=paste0('worldpop/',country.abbrev,'_u5_',svy.year,'_100m.tif',sep=''), 
            format="GTiff",overwrite=TRUE)

pop.u5.1km  <- aggregate(pop.u5, fact=10,fun=sum)
#sum(values(pop.fm1.1km),na.rm=TRUE)

pop.u5.1km.resample <- resample(pop.u5.1km, worldpop, method='bilinear')
#plot(pop.fm1.1km.resample)
#sum(values(pop.fm1.1km.resample),na.rm=TRUE)

writeRaster(pop.u5.1km.resample, file=paste0('worldpop/',country.abbrev,'_u5_',svy.year,'_1km.tif',sep=''), 
            format="GTiff",overwrite=TRUE)

}

u5.pop <- raster(paste0('worldpop/',country.abbrev,'_u5_',svy.year,'_1km.tif',sep=''))


f.15.49.pop <- raster(paste0('worldpop/',country.abbrev,'_f15_49_',svy.year,'_1km.tif',sep=''))
if(FALSE){
  
  setwd(paste(data_dir,country,sep=''))
  pop.f15 <-  raster(paste0('worldpop/',country.abbrev,'_f_15_',svy.year,'.tif',sep=''))
  pop.f20 <-  raster(paste0('worldpop/',country.abbrev,'_f_20_',svy.year,'.tif',sep=''))
  pop.f25 <-  raster(paste0('worldpop/',country.abbrev,'_f_25_',svy.year,'.tif',sep=''))
  

  pop.f15_29 <- pop.f15+pop.f20+pop.f25

  #plot(pop.f15_29)
  #sum(values(pop.f15_29),na.rm=T)
  writeRaster(pop.f15_29, file=paste0('worldpop/',country.abbrev,'_f15_29_',svy.year,'_100m.tif',sep=''), 
            format="GTiff",overwrite=TRUE)

  pop.f15_29.1km  <- aggregate(pop.f15_29, fact=10,fun=sum)
  #sum(values(pop.f15_29.1km),na.rm=TRUE)

  pop.f15_29.1km.resample <- resample(pop.f15_29.1km, worldpop, method='bilinear')
  #plot(pop.f15_29.1km.resample)
  #sum(values(pop.f15_29.1km.resample),na.rm=TRUE)

  writeRaster(pop.f15_29.1km.resample, file=paste0('worldpop/',country.abbrev,'_f15_29_',svy.year,'_1km.tif',sep=''), 
            format="GTiff",overwrite=TRUE)

}



if(FALSE){
  
  setwd(paste(data_dir,country,sep=''))
  pop.f15 <-  raster(paste0('worldpop/',country.abbrev,'_f_15_',svy.year,'.tif',sep=''))
  pop.f20 <-  raster(paste0('worldpop/',country.abbrev,'_f_20_',svy.year,'.tif',sep=''))
  pop.f25 <-  raster(paste0('worldpop/',country.abbrev,'_f_25_',svy.year,'.tif',sep=''))
  pop.f30 <-  raster(paste0('worldpop/',country.abbrev,'_f_30_',svy.year,'.tif',sep=''))
  pop.f35 <-  raster(paste0('worldpop/',country.abbrev,'_f_35_',svy.year,'.tif',sep=''))
  pop.f40 <-  raster(paste0('worldpop/',country.abbrev,'_f_40_',svy.year,'.tif',sep=''))
  pop.f45 <-  raster(paste0('worldpop/',country.abbrev,'_f_45_',svy.year,'.tif',sep=''))
  
  
  pop.f15_49 <- pop.f15+pop.f20+pop.f25+pop.f30+pop.f35+pop.f40+pop.f45
  
  #plot(pop.f15_49)
  #sum(values(pop.f15_49),na.rm=T)
  writeRaster(pop.f15_49, file=paste0('worldpop/',country.abbrev,'_f15_49_',svy.year,'_100m.tif',sep=''), 
              format="GTiff",overwrite=TRUE)
  
  pop.f15_49.1km  <- aggregate(pop.f15_49, fact=10,fun=sum)
  #sum(values(pop.f15_49.1km),na.rm=TRUE)
  
  pop.f15_49.1km.resample <- resample(pop.f15_49.1km, worldpop, method='bilinear')
  #plot(pop.f15_49.1km.resample)
  #sum(values(pop.f15_49.1km.resample),na.rm=TRUE)
  
  writeRaster(pop.f15_49.1km.resample, file=paste0('worldpop/',country.abbrev,'_f15_49_',svy.year,'_1km.tif',sep=''), 
              format="GTiff",overwrite=TRUE)
  
}


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

##################################################################
######### uncorrected BART
##################################################################

## load draws from BART unjitter using uncorrected locations

setwd(paste0(res_dir,country,'/BART_uncorrected',sep=''))
load(file=paste0('nolikoma_natl_pred_uncrc_draws_1000.rda'))

natl_pred_post[natl_pred_post> expit(15)] <-  expit(15)

get.uncali.surface(method= 'BART_uncrc',
                               natl.dat= urb_dat,
                               prob.frame = natl_pred_post,
                               ras.template = worldpop,
                               pop.ras = worldpop)
  
natl.cali.samp.pipeline(method= 'BART_uncrc',
                        natl.dat= urb_dat,
                        prob.frame = natl_pred_post,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)

urb_dat$bench_adm <- urb_dat$admin2.char

adm.cali.samp.pipeline(method= 'BART_uncrc',
                        natl.dat= urb_dat,
                        prob.frame = natl_pred_post,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        ref.tab=ref.tab,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)


##################################################################
######### corrected BART
##################################################################

## load draws from BART unjitter using uncorrected locations

setwd(paste0(res_dir,country,'/BART_corrected',sep=''))
load(file=paste0('nolikoma_natl_pred_crc_draws_1000.rda'))

natl_pred_post[natl_pred_post> expit(15)] <-  expit(15)

get.uncali.surface(method= 'BART_crc',
                   natl.dat= urb_dat,
                   prob.frame = natl_pred_post,
                   ras.template = worldpop,
                   pop.ras = worldpop)

natl.cali.samp.pipeline(method= 'BART_crc',
                        natl.dat= urb_dat,
                        prob.frame = natl_pred_post,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)

urb_dat$bench_adm <- urb_dat$admin2.char

adm.cali.samp.pipeline(method= 'BART_crc',
                       natl.dat= urb_dat,
                       prob.frame = natl_pred_post,
                       bench.frac = natl_urban_frac,
                       ras.template = worldpop,
                       pop.ras = worldpop,
                       ref.tab=ref.tab,
                       admin.poly=poly.adm2,
                       admin.level='admin2',
                       admin.name=admin2.names)

##################################################################
######### logistic regression with uncorrected locations
##################################################################

setwd(paste0(res_dir,country,'/logistic_regression',sep=''))
logis.post.draws <- readRDS('logistic_post_p_uncrc.rds')

logis.post.draws[logis.post.draws> expit(15)] <-  expit(15)
#rownames(logis.post.draws) <- c(1:nrow(logis.post.draws))

get.uncali.surface(method= 'logistic_uncrc',
                   natl.dat= urb_dat,
                   prob.frame = logis.post.draws,
                   ras.template = worldpop,
                   pop.ras = worldpop)

natl.cali.samp.pipeline(method= 'logistic_uncrc',
                        natl.dat= urb_dat,
                        prob.frame = logis.post.draws,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)

urb_dat$bench_adm <- urb_dat$admin2.char

adm.cali.samp.pipeline(method= 'logistic_uncrc',
                       natl.dat= urb_dat,
                       prob.frame = logis.post.draws,
                       bench.frac = natl_urban_frac,
                       ras.template = worldpop,
                       pop.ras = worldpop,
                       ref.tab=ref.tab,
                       admin.poly=poly.adm2,
                       admin.level='admin2',
                       admin.name=admin2.names)



### quick check
setwd(paste0(res_dir,country,'/UR/Fractions/','logistic_uncrc',sep=''))

load('logistic_uncrc_mse_natl_bench_total_admin2_frac.rda')

no.likoma.frac <- ref.tab[!ref.tab$GADM==
                            'Likoma',]
no.likoma.frac.sorted <- no.likoma.frac[order(no.likoma.frac$Internal),]

adm.frac.no.likoma <- frac.adm[[2]]
adm.frac.no.likoma <- adm.frac.no.likoma[!adm.frac.no.likoma$admin=='admin2_10',]

mean(abs(no.likoma.frac.sorted$urb_frac-adm.frac.no.likoma$Frac))


#setwd(paste0(res_dir,country,'/UR/Classification/','Calibration',sep=''))
#tmp.thres <- readRDS('logistic_uncrc_mse_natl_shift.rds')



##################################################################
######### logistic regression with corrected locations
##################################################################

setwd(paste0(res_dir,country,'/logistic_regression',sep=''))
logis.post.draws <- readRDS('logistic_post_p_crc.rds')

logis.post.draws[logis.post.draws> expit(15)] <-  expit(15)

get.uncali.surface(method= 'logistic_crc',
                   natl.dat= urb_dat,
                   prob.frame = logis.post.draws,
                   ras.template = worldpop,
                   pop.ras = worldpop)

natl.cali.samp.pipeline(method= 'logistic_crc',
                        natl.dat= urb_dat,
                        prob.frame = logis.post.draws,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)

urb_dat$bench_adm <- urb_dat$admin2.char

adm.cali.samp.pipeline(method= 'logistic_crc',
                       natl.dat= urb_dat,
                       prob.frame = logis.post.draws,
                       bench.frac = natl_urban_frac,
                       ras.template = worldpop,
                       pop.ras = worldpop,
                       ref.tab=ref.tab,
                       admin.poly=poly.adm2,
                       admin.level='admin2',
                       admin.name=admin2.names)



##################################################################
######### BART MCMC
##################################################################

if(FALSE){
setwd(paste0(res_dir,country,'/BART_jitter/natl_pred',sep=''))

nsamp=1000
nburn=500

load.MCMC.ite <- function(idx,nburn){
  
  load(paste0('natl_pred_ite_',idx+nburn,'.rda',sep=''))
  return(natl_pred_vec)
}
mcmc.draws.list<-lapply(c(1:nsamp), FUN=load.MCMC.ite,nburn=nburn)
mcmc.draws.frame <- do.call("cbind", mcmc.draws.list)

rm(mcmc.draws.list)

setwd(paste0(res_dir,country,'/BART_jitter/',sep=''))

saveRDS(mcmc.draws.frame,file='nolikoma_natl_pred_mcmc_draws_1000.rda')
}

setwd(paste0(res_dir,country,'/BART_jitter/',sep=''))

mcmc.draws.frame <- readRDS('nolikoma_natl_pred_mcmc_draws_1000.rda')
mcmc.draws.frame[mcmc.draws.frame> expit(15)] <-  expit(15)
mcmc.draws.frame <- as.data.frame(mcmc.draws.frame)

get.uncali.surface(method= 'BART_MCMC',
                   natl.dat= urb_dat,
                   prob.frame = mcmc.draws.frame,
                   ras.template = worldpop,
                   pop.ras = worldpop)

natl.cali.samp.pipeline(method= 'BART_MCMC',
                        natl.dat= urb_dat,
                        prob.frame = mcmc.draws.frame,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)

urb_dat$bench_adm <- urb_dat$admin2.char

adm.cali.samp.pipeline(method= 'BART_MCMC',
                       natl.dat= urb_dat,
                       prob.frame = mcmc.draws.frame,
                       bench.frac = natl_urban_frac,
                       ras.template = worldpop,
                       pop.ras = worldpop,
                       ref.tab=ref.tab,
                       admin.poly=poly.adm2,
                       admin.level='admin2',
                       admin.name=admin2.names)

#dim(mcmc.draws.frame)


##################################################################
######### GBT with uncorrected locations
##################################################################

setwd(paste0(res_dir,country,'/GBT/',sep=''))
load(file='gbt_pred_p_uncrc.rda')

gbt_pred_pvec <- natl_pred_vec
gbt_frame <- as.data.frame(gbt_pred_pvec)

get.uncali.surface(method= 'GBT_uncrc',
                   natl.dat= urb_dat,
                   prob.frame = gbt_frame,
                   ras.template = worldpop,
                   pop.ras = worldpop)

natl.cali.samp.pipeline(method= 'GBT_uncrc',
                        natl.dat= urb_dat,
                        prob.frame = gbt_frame,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)

urb_dat$bench_adm <- urb_dat$admin2.char

adm.cali.samp.pipeline(method= 'GBT_uncrc',
                       natl.dat= urb_dat,
                       prob.frame = gbt_frame,
                       bench.frac = natl_urban_frac,
                       ras.template = worldpop,
                       pop.ras = worldpop,
                       ref.tab=ref.tab,
                       admin.poly=poly.adm2,
                       admin.level='admin2',
                       admin.name=admin2.names)






##################################################################
######### GBT with corrected locations
##################################################################

setwd(paste0(res_dir,country,'/GBT/',sep=''))
load(file='gbt_pred_p_crc.rda')

gbt_pred_pvec <- natl_pred_vec
gbt_frame <- as.data.frame(gbt_pred_pvec)

get.uncali.surface(method= 'GBT_crc',
                   natl.dat= urb_dat,
                   prob.frame = gbt_frame,
                   ras.template = worldpop,
                   pop.ras = worldpop)

natl.cali.samp.pipeline(method= 'GBT_crc',
                        natl.dat= urb_dat,
                        prob.frame = gbt_frame,
                        bench.frac = natl_urban_frac,
                        ras.template = worldpop,
                        pop.ras = worldpop,
                        admin.poly=poly.adm2,
                        admin.level='admin2',
                        admin.name=admin2.names)

urb_dat$bench_adm <- urb_dat$admin2.char

adm.cali.samp.pipeline(method= 'GBT_crc',
                       natl.dat= urb_dat,
                       prob.frame = gbt_frame,
                       bench.frac = natl_urban_frac,
                       ras.template = worldpop,
                       pop.ras = worldpop,
                       ref.tab=ref.tab,
                       admin.poly=poly.adm2,
                       admin.level='admin2',
                       admin.name=admin2.names)
