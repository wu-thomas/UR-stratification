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
svy.year <- 2015

country.abbrev<-'mwi' #for population rasters
gadm.abbrev = "MWI"


set.seed(2022)

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

poly.adm2 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm2))
poly.adm3 <- readOGR(dsn = poly.path,encoding = "UTF-8", use_iconv = TRUE,
                     layer = as.character(poly.layer.adm3)) 

#plot(poly.adm3
proj4string(poly.adm0) <- proj4string(poly.adm1)  <- proj4string(poly.adm2)


load(paste0('shapeFiles_gadm/', country, '_Amat.rda'))
load(paste0('shapeFiles_gadm/', country, '_Amat_Names.rda'))


################################################################
######### prepare spatial
################################################################


Admin1Map <- poly.adm1
Admin2Map <- poly.adm2
Admin3Map <- poly.adm3



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
### A: fixed urban fraction from uncrc BART
### B: fixed urban fraction from BART-MCMC
### C: fixed urban fraction from BART-MCMC
################################################################

################################################################
######### prepare anlaysis data set
################################################################
setwd(paste0(data_dir,country,'/HIV',sep=''))
#save.image(file = 'HIV_mwi.RData')
load( file = "hiv_df_adm3.rda")


# remove missing values of responseVar so function will run -- DOUBLE CHECK THIS
hiv_df_no_na <- hiv_df %>% filter(!is.na(hiv_ind))
# only include young women
#hiv_df_subgroup <- hiv_df_no_na %>% filter(sex == "female" & age < 30) 
hiv_df_subgroup <- hiv_df_no_na %>% filter(sex == "female" & age < 50) 
#hiv_df_subgroup <- hiv_df_no_na

hiv_df_subgroup <- hiv_df_subgroup %>% mutate(strata = as.character(strata))



#setwd(paste0(data_dir,country,'/HIV',sep=''))

#load( file = "hiv_df_subgroup.rda")

hiv_dat <- hiv_df_subgroup
hiv_dat$HIV <- hiv_dat$hiv_ind
hiv_dat$stratum <- hiv_dat$strata

# Add strata
hiv_dat$adm1_ur <- paste0(hiv_dat$admin1.char,':',hiv_dat$urban,sep='')
hiv_dat$adm2_ur <- paste0(hiv_dat$admin2.char,':',hiv_dat$urban,sep='')

# add other variables
hiv_dat$weight <- hiv_dat$hiv05
hiv_dat$Ntrials <- 1

### load models
setwd(paste0(res_dir,country,'/HIV',sep=''))
#load(file = 'HIV_female_15_49_updated.RData')

load(file = 'HIV_female_15_49_updated.RData')


######   load helper functions

setwd(paste0(code_dir,'/',country,'/HIV'))
source('helper_functions_HIV_updated.R')

rm(admin2.direct.aggre.A)
rm(admin2.direct.aggre.B)
rm(admin2.direct.aggre.C)

rm(area.iid.strat.admin2.res.A)
rm(area.iid.strat.admin2.res.B)
rm(area.iid.strat.admin2.res.C)

rm(area.iid.strat.int.admin2.res.A)
rm(area.iid.strat.int.admin2.res.B)
rm(area.iid.strat.int.admin2.res.C)

rm(BYM2.strat.int.admin2.res.A)
rm(BYM2.strat.int.admin2.res.B)
rm(BYM2.strat.int.admin2.res.C)

rm(BYM2.strat.noint.admin2.res.A)
rm(BYM2.strat.noint.admin2.res.B)
rm(BYM2.strat.noint.admin2.res.C)

rm(strat.admin2.res.A)
rm(strat.admin2.res.B)
rm(strat.admin2.res.C)

rm(strat.int.admin2.res.A)
rm(strat.int.admin2.res.B)
rm(strat.int.admin2.res.C)



################################################################
#########  load urban fraction for female 15-49 population 
################################################################

## load point estimates of urban/rural fractions
method= 'BART_MCMC'
loss.func= 'mse'

setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))


### load admin-2 urban fractions
load(file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin2','_frac.rda'))

adm2.urb.fixed <- frac.adm$post.median
adm2.urb.fixed[adm2.urb.fixed$admin=='admin2_10',]$Frac<- 692/5497
adm2.urb.fixed <- merge(adm2.urb.fixed,admin2.names,by.x='admin',by.y='Internal')
adm2.urb.fixed <- adm2.urb.fixed %>% rename('admin.name' = 'GADM')
adm2.urb.fixed <- adm2.urb.fixed %>% rename('adm.char' = 'admin')
adm2.urb.fixed <- adm2.urb.fixed[match(admin2.names$Internal,adm2.urb.fixed$adm.char),]


adm2.urb.samp <- as.data.frame(frac.adm$post.sample)
adm2.urb.samp$admin2_10 <- 692/5497
adm2.urb.samp  <- adm2.urb.samp[,admin2.names$Internal]

### load admin3 urban fractions
load(file=paste0(method,'_',loss.func,'_adm_bench_f_15_49_','admin3','_frac.rda'))

adm3.urb.fixed <- frac.adm$post.median
adm3.urb.fixed <- merge(adm3.urb.fixed,admin3.names,by.x='admin',by.y='Internal')
adm3.urb.fixed <- adm3.urb.fixed %>% rename('admin.name' = 'GADM')
adm3.urb.fixed <- adm3.urb.fixed %>% rename('adm.char' = 'admin')
adm3.urb.fixed <- adm3.urb.fixed[match(admin3.names$Internal,adm3.urb.fixed$adm.char),]

adm3.urb.samp <- as.data.frame(frac.adm$post.sample)
adm3.urb.samp  <- adm3.urb.samp[,admin3.names$Internal]

################################################################
#########  admin-2, direct estimates
################################################################

# Compute estimates
admin2.direct = getHIVDirectAdm2(hiv_dat,save.draws = T)
admin2.direct.ur = getHIVDirectAdm2(hiv_dat, byUrban = TRUE,save.draws = T)

admin2.direct.est = admin2.direct$overall
admin2.direct.est.ur = admin2.direct.ur$overall


# admin2.direct.aggre <- aggreDirectAdm2(urb.res=admin2.direct.est.ur,urb.frac=adm2.urb.fixed)

# fixed fractions (medians) CIs
admin2.direct.aggre.A <- aggreDirectAdm2_CI(ur.draws=admin2.direct.ur$draws,urb.frac=adm2.urb.fixed,
                                            admin.ref = admin2.names)


# sample fractions CIs
admin2.direct.aggre.B <- aggreDirectAdm2_CI(ur.draws=admin2.direct.ur$draws,urb.frac=adm2.urb.samp,
                                             admin.ref = admin2.names,fixed.urb=F)


if(FALSE){
#svyby(~measles, ~admin1, DHSdesign, svymean, vartype=c("se","ci"))
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=admin2.direct.aggre.B$p_Med)
abline(a=0,b=1)

par(pty="s")

plot(x=admin2.direct.est$p_Med,y=admin2.direct.aggre.samp.ur$p_Med)
abline(a=0,b=1)


mean(admin2.direct.aggre.A$p_Med-admin2.direct.est$p_Med)
mean(admin2.direct.aggre.B$p_Med-admin2.direct.est$p_Med)

}



################################################################
#########  admin-2 fixed effects nonstratified cluster level model
################################################################
# Compute admin-2 nonstratified cluster level model
iidPrior  = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.05)))
inla.nonstrat.admin2 = FitAdmin2Nonstrat(myData = hiv_dat,
                                         clustPrior = iidPrior)

summary(inla.nonstrat.admin2)
# Compute estimates

nonstrat.admin2.res = aggAdmin2Nonstrat(res.inla = inla.nonstrat.admin2,
                                        admin.ref = admin2.names,
                                        nSamp = 1000)


par(pty="s")

plot(x=admin2.direct.est$p_Med,y=nonstrat.admin2.res$p_Med)
abline(a=0,b=1)

mean(nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med)


################################################################
#########  admin-2 fixed effects stratified cluster level model, no int
################################################################


iidPrior  = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.05)))

inla.strat.admin2 = FitAdmin2Strat(myData = hiv_dat,
                                   clustPrior = iidPrior)

summary(inla.strat.admin2)


# get aggregated estimates 
strat.admin2.res.A = aggAdmin2Strat(res.inla = inla.strat.admin2,
                                  nSamp = 1000,
                                  urb.frac = adm2.urb.fixed,
                                  admin.ref = admin2.names)

strat.admin2.res.B = aggAdmin2Strat(res.inla = inla.strat.admin2,
                                    nSamp = 1000,
                                    urb.frac = adm2.urb.samp,
                                    fixed.urb = F,
                                    admin.ref = admin2.names)



################################################################
#########  admin-2 fixed effects cluster level model strat, int
################################################################

inla.strat.int.admin2 = FitAdmin2StratInt(myData = hiv_dat,
                                          clustPrior = iidPrior)


summary(inla.strat.int.admin2)


strat.int.admin2.res.A = aggAdmin2StratInt(res.inla = inla.strat.int.admin2,
                                         nSamp = 1000,
                                         urb.frac = adm2.urb.fixed,
                                         admin.ref=admin2.names)

strat.int.admin2.res.B = aggAdmin2StratInt(res.inla = inla.strat.int.admin2,
                                         nSamp = 1000,
                                         fixed.urb = F,
                                         urb.frac = adm2.urb.samp,
                                         admin.ref=admin2.names)
if(FALSE){
### compare overall estimates 
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=strat.int.admin2.res.A$Combined.est$p_Med)
abline(a=0,b=1)

mean(strat.int.admin2.res.A$Combined.est$p_Med-admin2.direct.est$p_Med)

par(pty="s")

plot(x=admin2.direct.est$p_Med,y=strat.int.admin2.res.C$Combined.est$p_Med)
abline(a=0,b=1)

mean(strat.int.admin2.res.C$Combined.est$p_Med-admin2.direct.est$p_Med)


par(pty="s")

plot(strat.int.admin2.res$Combined.est$p_Med,strat.int.admin2.res.CI$Combined.est$p_Med)
abline(a=0,b=1)


mean(admin2.direct.est$p_Upp-admin2.direct.est$p_Low)
mean(strat.int.admin2.res.C$Combined.est$p_Upp-strat.int.admin2.res.C$Combined.est$p_Low)
mean(strat.int.admin2.res.A$Combined.est$p_Upp-strat.int.admin2.res.A$Combined.est$p_Low)

}

################################################################
#########  admin-2 random effects cluster level nonstratified model
################################################################



### non-stratified model

# fit model

iidPrior  = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.05)))


area.iid.nonstrat.admin2 = FitAdmin2IIDNonstrat(myData = hiv_dat,
                                                clustPrior = iidPrior)


summary(area.iid.nonstrat.admin2)

# posterior sample
area.iid.nonstrat.admin2.res = aggAdmin2IIDNonstrat(res.inla = area.iid.nonstrat.admin2,
                                                    nSamp = 1000,admin.ref=admin2.names)

if(FALSE){
# comparison
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=area.iid.nonstrat.admin2.res$p_Med)
abline(a=0,b=1)

mean(area.iid.nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med)


setwd(paste0(res_dir,country,'/HIV/rev_res',sep=''))
save.image(file = paste0('HIV_female_15_49_',method,'_',loss.func,'.RData'))
}


################################################################
#########  admin-2 random effects cluster level stratified, no int
################################################################


### stratified no interaction
area.iid.strat.admin2 = FitAdmin2IIDStrat(myData = hiv_dat,
                                          clustPrior = iidPrior)

area.iid.strat.admin2.res.A = aggAdmin2IIDStrat(res.inla = area.iid.strat.admin2,
                                              nSamp = 1000,
                                              urb.frac = adm2.urb.fixed,
                                              admin.ref = admin2.names)

area.iid.strat.admin2.res.B = aggAdmin2IIDStrat(res.inla = area.iid.strat.admin2,
                                              nSamp = 1000,
                                              fixed.urb = F,
                                              urb.frac = adm2.urb.samp,
                                              admin.ref=admin2.names)

if(FALSE){
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=area.iid.strat.admin2.res.C$Combined.est$p_Med,
     xlim=c(0,0.15),ylim=c(0,0.15))
abline(a=0,b=1)

mean(area.iid.strat.admin2.res$Combined.est$p_Med-admin2.direct.est$p_Med)


par(pty="s")

plot(x=admin2.direct.est$p_Med,y=area.iid.strat.admin2.res.CI$Combined.est$p_Med)
abline(a=0,b=1)

mean(area.iid.strat.admin2.res.CI$Combined.est$p_Med-admin2.direct.est$p_Med)



mean(admin2.direct.est$p_Upp-admin2.direct.est$p_Low)
mean(strat.int.admin2.res$Combined.est$p_Upp-strat.int.admin2.res$Combined.est$p_Low)
mean(strat.int.admin2.res.CI$Combined.est$p_Upp-strat.int.admin2.res.CI$Combined.est$p_Low)
mean(area.iid.strat.admin2.res.CI$Combined.est$p_Upp-area.iid.strat.admin2.res.CI$Combined.est$p_Low)
}

################################################################
#########  admin-2 random effects cluster level stratified,  int
################################################################


area.iid.strat.int.admin2 = FitAdmin2IIDStratInt(myData = hiv_dat,
                                                 clustPrior = iidPrior)


summary(area.iid.strat.int.admin2)

area.iid.strat.int.admin2.res.A = aggAdmin2IIDStratInt(res.inla = area.iid.strat.int.admin2,
                                                     nSamp = 1000,
                                                     urb.frac = adm2.urb.fixed,
                                                     admin.ref = admin2.names)


area.iid.strat.int.admin2.res.B = aggAdmin2IIDStratInt(res.inla = area.iid.strat.int.admin2,
                                                     nSamp = 1000,
                                                     fixed.urb = F,
                                                     urb.frac = adm2.urb.samp,
                                                     admin.ref=admin2.names)



if(FALSE){
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=area.iid.strat.int.admin2.res.C$Combined.est$p_Med)
abline(a=0,b=1)

mean(area.iid.strat.int.admin2.res.A$Combined.est$p_Med-admin2.direct.est$p_Med)

mean(admin2.direct.est$p_Upp-admin2.direct.est$p_Low)
mean(strat.int.admin2.res$Combined.est$p_Upp-strat.int.admin2.res$Combined.est$p_Low)
mean(strat.int.admin2.res.CI$Combined.est$p_Upp-strat.int.admin2.res.CI$Combined.est$p_Low)
mean(area.iid.strat.int.admin2.res.CI$Combined.est$p_Upp-area.iid.strat.int.admin2.res.CI$Combined.est$p_Low)
}

################################################################
#########  admin-2 BYM2 cluster level nonstratified model
################################################################

# Set priors
bym2prior = list(prec = list(param = c(1, 0.05)),
                 phi  = list(param = c(0.5, 0.5)))
iidPrior  = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.05)))


BYM2.nonstrat.admin2 = getBYMAdmin2(myData = hiv_dat,
                                    Admin2Graph = Admin2Graph,
                                    bym2prior = bym2prior,
                                    clustPrior = iidPrior,
                                    strata=FALSE)

BYM2.nonstrat.admin2.res <- aggBYMAdmin2(res.inla=BYM2.nonstrat.admin2, 
                                         nSamp = 1000,
                                         urb.frac=NULL,strata=FALSE,
                                         admin.ref=admin2.names)

if(FALSE){
# compare with observed 
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=BYM2.nonstrat.admin2.res$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))

abline(a=0,b=1)

mean(BYM2.nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med)

mean(abs(BYM2.nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med))

}

################################################################
#########  admin-2 BYM2 cluster level stratified model
################################################################

# Set priors
bym2prior = list(prec = list(param = c(1, 0.05)),
                 phi  = list(param = c(0.5, 0.5)))
iidPrior  = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.05)))


BYM2.strat.noint.admin2 = getBYMAdmin2(myData = hiv_dat,
                                       Admin2Graph = Admin2Graph,
                                       bym2prior = bym2prior,
                                       clustPrior = iidPrior)

summary(BYM2.strat.noint.admin2)


BYM2.strat.noint.admin2.res.A <- aggBYMAdmin2(res.inla=BYM2.strat.noint.admin2, 
                                           nSamp = 1000,
                                            urb.frac=adm2.urb.fixed,int=FALSE,
                                            fixed.urb = T,
                                            admin.ref=admin2.names)


BYM2.strat.noint.admin2.res.B <- aggBYMAdmin2(res.inla=BYM2.strat.noint.admin2, 
                                            nSamp = 1000,
                                            urb.frac=adm2.urb.samp,int=FALSE,
                                            fixed.urb = F,
                                            admin.ref=admin2.names)

if(FALSE){
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=BYM2.strat.noint.admin2.res.A$Combined.est$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))
abline(a=0,b=1)


mean(BYM2.strat.noint.admin2.res.CI$Combined.est$p_Upp-BYM2.strat.noint.admin2.res.CI$Combined.est$p_Low)
}

################################################################
#########  admin-2 BYM2 cluster level stratified model interactions
################################################################
# Compute full estimate
bym2prior = list(prec = list(param = c(1, 0.05)),
                 phi  = list(param = c(0.5, 0.5)))
iidPrior  = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.05)))


BYM2.strat.int.admin2 = getBYMAdmin2(myData = hiv_dat,
                                     Admin2Graph = Admin2Graph,
                                     bym2prior = bym2prior,
                                     clustPrior = iidPrior,int=TRUE)

summary(BYM2.strat.int.admin2)

BYM2.strat.int.admin2.res.A <- aggBYMAdmin2(res.inla=BYM2.strat.int.admin2, 
                                          nSamp = 1000,
                                          urb.frac=adm2.urb.fixed,int=TRUE,
                                          fixed.urb = T,
                                          admin.ref=admin2.names)

BYM2.strat.int.admin2.res.B <- aggBYMAdmin2(res.inla=BYM2.strat.int.admin2, 
                                            nSamp = 1000,
                                            urb.frac=adm2.urb.samp,int=TRUE,
                                            fixed.urb = F,
                                            admin.ref=admin2.names)
# Calculate estimates

if(FALSE){
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=BYM2.strat.int.admin2.res.A$Combined.est$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))
abline(a=0,b=1)


mean(BYM2.strat.int.admin2.res.A$Combined.est$p_Med-admin2.direct.est$p_Med)
}

if(FALSE){
#setwd(paste0(res_dir,country,'/HIV/rev_res',sep=''))
#save.image(file = paste0('admin2_HIV_female_15_49_',method,'_',loss.func,'.RData'))


rm(area.iid.nonstrat.admin2)
rm(area.iid.strat.admin2)
rm(area.iid.strat.int.admin2)
rm(inla.nonstrat.admin2)
rm(inla.strat.admin2)
rm(inla.strat.int.admin2)
rm(BYM2.nonstrat.admin2)
rm(BYM2.strat.noint.admin2)
rm(BYM2.strat.int.admin2)

#setwd(paste0(res_dir,country,'/HIV/rev_res',sep=''))
#save.image(file = paste0('res_admin2_HIV_female_15_49_',method,'_',loss.func,'.RData'))
}
################################################################
#########  admin-3 BYM2 cluster level nonstratified model 
################################################################

# Compute full estimate
bym2prior = list(prec = list(param = c(1, 0.05)),
                 phi  = list(param = c(0.5, 0.5)))
iidPrior  = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.05)))


BYM2.nonstrat.admin3 = getBYMAdmin3(myData = hiv_dat,
                                    Admin3Graph = Admin3Graph,
                                    bym2prior = bym2prior,
                                    clustPrior = iidPrior,
                                    strata=FALSE)

BYM2.nonstrat.admin3.res <- aggBYMAdmin2(res.inla=BYM2.nonstrat.admin3, 
                                         nSamp = 1000,
                                         urb.frac=NULL,strata=FALSE,
                                         admin.ref=admin3.names)


summary(BYM2.nonstrat.admin3)

################################################################
#########  admin-3 BYM2 cluster level stratified model 
################################################################


BYM2.strat.admin3 = getBYMAdmin3(myData = hiv_dat,
                                    Admin3Graph = Admin3Graph,
                                    bym2prior = bym2prior,
                                    clustPrior = iidPrior,
                                    strata=TRUE)


BYM2.strat.admin3.res.A <- aggBYMAdmin2(res.inla=BYM2.strat.admin3, 
                                      nSamp = 1000,
                                      urb.frac=adm3.urb.fixed,int=FALSE,
                                      fixed.urb = T,
                                      admin.ref=admin3.names)

BYM2.strat.admin3.res.B <- aggBYMAdmin2(res.inla=BYM2.strat.admin3, 
                                        nSamp = 1000,
                                        urb.frac=adm3.urb.samp,int=FALSE,
                                        fixed.urb = F,
                                        admin.ref=admin3.names)

par(pty="s")

plot(x=BYM2.strat.admin3.res.B$Combined.est$p_Med,y=BYM2.nonstrat.admin3.res$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))
abline(a=0,b=1)

mean(BYM2.strat.admin3.res.B$Combined.est$p_Med-BYM2.nonstrat.admin3.res$p_Med)



par(pty="s")

plot(x=BYM2.strat.noint.admin2.res.A$Combined.est$p_Med,y=BYM2.nonstrat.admin2.res$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))
abline(a=0,b=1)



################################################################
#########  save results
################################################################

#setwd(paste0(res_dir,country,'/HIV/rev_res',sep=''))
#save.image(file = paste0('HIV_female_15_49_',method,'_',loss.func,'.RData'))



#rm(BYM2.nonstrat.admin3)
#rm(BYM2.strat.admin3)

#setwd(paste0(res_dir,country,'/HIV/rev_res',sep=''))
#save.image(file = paste0('res_HIV_female_15_49_',method,'_',loss.func,'.RData'))

setwd(paste0(res_dir,country,'/HIV/',sep=''))

#load('')