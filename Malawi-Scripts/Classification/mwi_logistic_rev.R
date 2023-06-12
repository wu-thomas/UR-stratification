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

if(!dir.exists(paths = paste0('logistic_regression'))){
  dir.create(path = paste0('logistic_regression'))
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
######### load national grid 
################################################################

setwd(paste0(data_dir,country,sep=''))
load("prepared_dat/natl_grid.rda")

urb_dat$complete <- complete.cases(urb_dat) 
urb_dat$pixel_id <- c(1:dim(urb_dat)[1])

complete_id <- complete.cases(urb_dat) 
complete_urb_dat <- urb_dat[complete_id,]

################################################################
#########  load cluster training data
################################################################

setwd(paste0(data_dir,country,'/prepared_dat/',sep=''))
load(file='crc_dat_nolikoma.rda')
load(file='uncrc_dat_nolikoma.rda')


################################################################
#########  standardize ndvi, tavg
################################################################

set.seed(2022)

ndvi_mean <- mean(complete_urb_dat$ndvi,na.rm=T)
ndvi_sd <- sd(complete_urb_dat$ndvi,na.rm=T)

pop_mean <- mean(complete_urb_dat$pop_den,na.rm=T)
pop_sd <- sd(complete_urb_dat$pop_den,na.rm=T)

tavg_mean <- mean(complete_urb_dat$tavg,na.rm=T)
tavg_sd <- sd(complete_urb_dat$tavg,na.rm=T)

access_mean <- mean(complete_urb_dat$access,na.rm=T)
access_sd <- sd(complete_urb_dat$access,na.rm=T)

prec_mean <- mean(complete_urb_dat$prec,na.rm=T)
prec_sd <- sd(complete_urb_dat$prec,na.rm=T)


elev_mean <- mean(complete_urb_dat$elev,na.rm=T)
elev_sd <- sd(complete_urb_dat$elev,na.rm=T)

light_mean <- mean(complete_urb_dat$light,na.rm=T)
light_sd <- sd(complete_urb_dat$light,na.rm=T)


complete_urb_dat$ndvi_std <- (complete_urb_dat$ndvi-ndvi_mean)/ndvi_sd
complete_urb_dat$tavg_std <- (complete_urb_dat$tavg-tavg_mean)/tavg_sd
complete_urb_dat$pop_std <- (complete_urb_dat$pop_den-pop_mean)/pop_sd
complete_urb_dat$access_std <- (complete_urb_dat$access-access_mean)/access_sd
complete_urb_dat$prec_std <- (complete_urb_dat$prec-prec_mean)/prec_sd
complete_urb_dat$elev_std <- (complete_urb_dat$elev-elev_mean)/elev_sd
complete_urb_dat$light_std <- (complete_urb_dat$light-light_mean)/light_sd


complete_urb_dat$urban_ind <- NA

complete_urb_dat <-  complete_urb_dat[sample(1:nrow(complete_urb_dat)), ]


num_groups <- 10
natl_comp_group<-complete_urb_dat %>% 
  group_by((row_number()-1) %/% (n()/num_groups)) %>%
  nest %>% pull(data)

################################################################
#########  classification using uncorrected locations 
################################################################
nsamp=1000

cov_list<-c('pop_std','light_std','elev_std','prec_std','access_std','ndvi_std','tavg_std','admin2.char')

train_dat <- uncrc_dat_nolikoma
train_dat$urban <- relevel(as.factor(as.character(train_dat$urban)),ref='urban')
train_dat$urban_ind <- ifelse(train_dat$urban == "urban", 1,0)

train_dat$ndvi_std <- (train_dat$ndvi-ndvi_mean)/ndvi_sd
train_dat$tavg_std <- (train_dat$tavg-tavg_mean)/tavg_sd
train_dat$pop_std <- (train_dat$pop_den-pop_mean)/pop_sd
train_dat$access_std <- (train_dat$access-access_mean)/access_sd
train_dat$prec_std <- (train_dat$prec-prec_mean)/prec_sd
train_dat$elev_std <- (train_dat$elev-elev_mean)/elev_sd
train_dat$light_std <- (train_dat$light-light_mean)/light_sd


f0 <- as.formula(
  paste('urban_ind', 
        paste(cov_list, collapse = " + "), 
        sep = " ~ "))
print(f0)


natl.pred.draws <- data.frame()


for (i in 1:num_groups){
  
  print(i)
  train_pred <- rbind (train_dat[,c('urban_ind',cov_list)],
                       natl_comp_group[[i]][,c('urban_ind',cov_list),])
  
  
  
  logis.mod<-inla(formula=f0,data=train_pred[,c('urban_ind',cov_list),],
                  family = "binomial", Ntrials=1 , #for bernoulli, logistic regression
                  control.family=list(link='logit'),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,
                                       config = TRUE))
  
  #summary(logis.mod)
  post.sample = inla.posterior.sample(n = nsamp, result = logis.mod)
  
  pred.idx <- (dim(train_dat)[1]+1):(dim(natl_comp_group[[i]])[1]+dim(train_dat)[1])
  
  get.pred.draws <- function(idx){
    one.draw <- post.sample[[idx]]
    return(one.draw$latent[pred.idx])
  }
  
  pred.draws.frame.list<-lapply(c(1:nsamp), FUN=get.pred.draws)
  pred.draws.frame <- do.call("cbind", pred.draws.frame.list)
  
  natl.pred.draws <- rbind (natl.pred.draws,pred.draws.frame)
  
  #logis.mod$all.hyper
  #dim(logis.mod$summary.fitted.values)
  
  
}

#save(fitted_natl,file=paste0('logistic_pred_p_uncrc.rda'))
natl.pred.draws <- expit(natl.pred.draws)
match.order <- match(sort(complete_urb_dat$pixel_id),complete_urb_dat$pixel_id)

natl.pred.draws.ordered <- natl.pred.draws[match.order,]

setwd(paste0(res_dir,country,'/logistic_regression',sep=''))
saveRDS(natl.pred.draws.ordered,file=paste0('logistic_post_p_uncrc.rds'))





################################################################
#########  classification using corrected locations 
################################################################
nsamp=1000

cov_list<-c('pop_std','light_std','elev_std','prec_std','access_std','ndvi_std','tavg_std','admin2.char')

train_dat <- crc_dat_nolikoma
train_dat$urban <- relevel(as.factor(as.character(train_dat$urban)),ref='urban')
train_dat$urban_ind <- ifelse(train_dat$urban == "urban", 1,0)

train_dat$ndvi_std <- (train_dat$ndvi-ndvi_mean)/ndvi_sd
train_dat$tavg_std <- (train_dat$tavg-tavg_mean)/tavg_sd
train_dat$pop_std <- (train_dat$pop_den-pop_mean)/pop_sd
train_dat$access_std <- (train_dat$access-access_mean)/access_sd
train_dat$prec_std <- (train_dat$prec-prec_mean)/prec_sd
train_dat$elev_std <- (train_dat$elev-elev_mean)/elev_sd
train_dat$light_std <- (train_dat$light-light_mean)/light_sd


f0 <- as.formula(
  paste('urban_ind', 
        paste(cov_list, collapse = " + "), 
        sep = " ~ "))
print(f0)


natl.pred.draws <- data.frame()


for (i in 1:num_groups){
  
  print(i)
  train_pred <- rbind (train_dat[,c('urban_ind',cov_list)],
                       natl_comp_group[[i]][,c('urban_ind',cov_list),])
  
  
  
  logis.mod<-inla(formula=f0,data=train_pred[,c('urban_ind',cov_list),],
                  family = "binomial", Ntrials=1 , #for bernoulli, logistic regression
                  control.family=list(link='logit'),
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,
                                       config = TRUE))
  
  #summary(logis.mod)
  post.sample = inla.posterior.sample(n = nsamp, result = logis.mod)
  
  pred.idx <- (dim(train_dat)[1]+1):(dim(natl_comp_group[[i]])[1]+dim(train_dat)[1])
  
  get.pred.draws <- function(idx){
    one.draw <- post.sample[[idx]]
    return(one.draw$latent[pred.idx])
  }
  
  pred.draws.frame.list<-lapply(c(1:nsamp), FUN=get.pred.draws)
  pred.draws.frame <- do.call("cbind", pred.draws.frame.list)
  
  natl.pred.draws <- rbind (natl.pred.draws,pred.draws.frame)
  
  #logis.mod$all.hyper
  #dim(logis.mod$summary.fitted.values)
  
  
}

#save(fitted_natl,file=paste0('logistic_pred_p_uncrc.rda'))
natl.pred.draws <- expit(natl.pred.draws)
match.order <- match(sort(complete_urb_dat$pixel_id),complete_urb_dat$pixel_id)

natl.pred.draws.ordered <- natl.pred.draws[match.order,]

setwd(paste0(res_dir,country,'/logistic_regression',sep=''))
saveRDS(natl.pred.draws.ordered,file=paste0('logistic_post_p_crc.rds'))

