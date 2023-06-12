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

setwd(paste0(code_dir,'/',country,'/HIV'))
source('helper_functions_HIV_updated.R')


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

Admin2Map <- poly.adm2

Admin1Map <- poly.adm1


# Get graph
nameVec = Admin2Map$NAME_1
for(i in 1:length(nameVec)){
  nameVec[i] = paste(nameVec[i], ":", Admin2Map$NAME_2[i], sep = "")
}
Admin2Graph = getAmat(Admin2Map, names = nameVec)
Admin1Graph = getAmat(Admin1Map, Admin1Map$NAME_1)

# Make INLA graph object
Admin2GraphINLA = inla.read.graph(Admin2Graph)

################################################################
#########  load urban fraction for female 15-49 population 
################################################################

## load point estimates of urban/rural fractions

setwd(paste0(res_dir,country,'/UR/Fractions/',sep=''))
#load('uncrc_adm2_thresh_f15_49_natl_frac.rda')
#load('uncrc_adm2_thresh_f15_49_adm1_frac.rda')
load('uncrc_adm2_thresh_f15_49_adm2_frac.rda')

f15_49.adm2.frac[f15_49.adm2.frac$admin.name=='Likoma',]$Frac<- 692/5497


## load samples of urban/rural fractions


setwd(paste0(res_dir,country,'/UR/Fractions/Samples/',sep=''))
method <-'mcmc'

adm2.urb.samp <- readRDS(file=paste0(method,'_f_adm2_frac','_adm2_thresh','.rds',sep=''))


adm2.urb.samp$admin2_10 <- 692/5497




################################################################
######### prepare anlaysis data set
################################################################
setwd(paste0(data_dir,country,'/HIV',sep=''))
#save.image(file = 'HIV_mwi.RData')
load( file = "hiv_df.rda")


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

################################################################
#########  national, direct estimates
################################################################


DHSdesign <- svydesign(id= ~clusterIdx,
                       strata=~strata,
                       weights=~hiv05, data= hiv_dat)
#DHSdesign <- svydesign(id= ~clusterIdx + householdIdx,
##                          strata=~stratum, nest=T, 
#                          weights=~weight, data=SData)

natl.direct.ur <- svyby(~ HIV, ~urban, DHSdesign, svymean, vartype=c("se","ci"))

## national survey direct estimate
svymean(~HIV, DHSdesign)
mean(hiv_dat$HIV)

## national urban rural aggregated direct estimates
natl.direct.ur$HIV[1]*(1-f15_49.natl.frac)+natl.direct.ur$HIV[2]*f15_49.natl.frac

#svyby(~measles, ~admin1*urban, DHSdesign, svymean, vartype=c("se","ci"))

#confint(svymean(~HIV, DHSdesign))

# national  36.57% (34.93%, 38.21%)
# mean(SData$stunt)

table(hiv_dat$urban)


################################################################
#########   test direct estimates uncertainty
################################################################
myData <- hiv_dat
myData$adm2_ur <- as.factor(myData$adm2_ur )
# Get direct stratum-specific direct estimates
my.svydesign <- svydesign(id= ~clusterIdx ,
                          strata=~stratum, nest=T, 
                          weights=~weight, data=myData)

mean_tmp <- survey::svyby(formula=~HIV, by=~adm2_ur, design=my.svydesign, survey::svymean, drop.empty.groups=FALSE)

name.i <- mean$region0            
p.i <- mean_tmp$HIV
var.i <- mean_tmp$se^2
  ht <- log(p.i/(1-p.i))
  ht.v <- sqrt(var.i/(p.i^2*(1-p.i)^2))
  ht.prec <- 1/ht.v
  
  
  ur.draws <- admin2.direct.ur$draws
  
post.sample = inla.posterior.sample(n = 1000, result = inla.strat.admin2)

# Initialize storage
n_admin2 = 28
est.U =  matrix(0, nrow = 1000, ncol = n_admin2)
est.R =  matrix(0, nrow = 1000, ncol = n_admin2)
est.admin2 = matrix(0, nrow = 1000, ncol = n_admin2)

# Get indicies in the sample
areaIdx  = inla.strat.admin2$misc$configs$contents$start[seq(3,2+n_admin2)]
urbIdx  = inla.strat.admin2$misc$configs$contents$start[3+n_admin2]

# Sample urban and rural
nugSimStd = rnorm(1e5, mean = 0, sd = 1)
for(i in 1:1000){
  for(j in 1:n_admin2){
    # Nugget is measurement error
    etaUrb.tmp = post.sample[[i]]$latent[areaIdx[j]] + post.sample[[i]]$latent[urbIdx]
    etaRur.tmp = post.sample[[i]]$latent[areaIdx[j]]
    
    # Nugget is overdispersion
    cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
    est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
    est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
}


tmp.check <- aggre_Adm_samp(u.draws=expit(est.U),r.draws=expit(est.R),urb.frac=f15_49.adm2.frac,fixed.urb=T,nSamp=1000,
               admin.ref=admin2.names)


par(pty="s")

plot(admin2.direct.est$p_Med,tmp.check$p_Med)
abline(a=0,b=1)


################################################################
#########  admin-2, direct estimates
################################################################

# Compute estimates
admin2.direct = getHIVDirectAdm2(hiv_dat,save.draws = T)
admin2.direct.ur = getHIVDirectAdm2(hiv_dat, byUrban = TRUE,save.draws = T)

admin2.direct.est = admin2.direct$overall
admin2.direct.est.ur = admin2.direct.ur$overall


admin2.direct.aggre <- aggreDirectAdm2(urb.res=admin2.direct.est.ur,urb.frac=f15_49.adm2.frac)

# fixed fractions CIs
admin2.direct.aggre.CI <- aggreDirectAdm2_CI(ur.draws=admin2.direct.ur$draws,urb.frac=f15_49.adm2.frac)

# sample fractions CIs
admin2.direct.aggre.samp.ur <- aggreDirectAdm2_CI(ur.draws=admin2.direct.ur$draws,urb.frac=adm2.urb.samp,
                                             admin.ref = admin2.names,fixed.urb=F)


#svyby(~measles, ~admin1, DHSdesign, svymean, vartype=c("se","ci"))
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=admin2.direct.aggre.CI$p_Med)
abline(a=0,b=1)

par(pty="s")

plot(x=admin2.direct.est$p_Med,y=admin2.direct.aggre.samp.ur$p_Med)
abline(a=0,b=1)


mean(admin2.direct.aggre.samp.ur$p_Med-admin2.direct.est$p_Med)






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
                                        myData  = hiv_dat,
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
strat.admin2.res = aggAdmin2Strat(res.inla = inla.strat.admin2,
                                  myData = hiv_dat,
                                  nSamp = 1000,
                                  urb.frac = f15_49.adm2.frac,
                                  admin.ref = admin2.names)


par(pty="s")

plot(x=admin2.direct.est$p_Med,y=strat.admin2.res$Combined.est$p_Med)
abline(a=0,b=1)

mean(strat.admin2.res$Combined.est$p_Med-admin2.direct.est$p_Med)

#mean(abs(strat.admin2.res$Combined.est$p_Med-admin2.direct.est$p_Med))
#mean(abs(nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med))


urban.direct.adm2 <- admin2.direct.est.ur[admin2.direct.est.ur$urb=='U',]
par(pty="s")
plot(x=urban.direct.adm2$p_Med,y=strat.admin2.res$Urban.res$p_Med)
abline(a=0,b=1)


par(pty="s")
plot(x=admin2.direct.est.ur[admin2.direct.est.ur$urb=='R',]$p_Med,y=strat.admin2.res$Rural.res$p_Med)
abline(a=0,b=1)

################################################################
#########  admin-2 fixed effects cluster level model strat, int
################################################################

inla.strat.int.admin2 = FitAdmin2StratInt(myData = hiv_dat,
                                          clustPrior = iidPrior)


summary(inla.strat.int.admin2)


strat.int.admin2.res = aggAdmin2StratInt(res.inla = inla.strat.int.admin2,
                                         myData = hiv_dat,
                                         nSamp = 1000,
                                         urb.frac = f15_49.adm2.frac,
                                         admin.ref=admin2.names)

strat.int.admin2.res.CI = aggAdmin2StratInt(res.inla = inla.strat.int.admin2,
                                         myData = hiv_dat,
                                         nSamp = 1000,
                                         fixed.urb = F,
                                         urb.frac = adm2.urb.samp,
                                         admin.ref=admin2.names)

### compare overall estimates 
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=strat.int.admin2.res$Combined.est$p_Med)
abline(a=0,b=1)

mean(strat.int.admin2.res$Combined.est$p_Med-admin2.direct.est$p_Med)

par(pty="s")

plot(x=admin2.direct.est$p_Med,y=strat.int.admin2.res.CI$Combined.est$p_Med)
abline(a=0,b=1)

mean(strat.int.admin2.res.CI$Combined.est$p_Med-admin2.direct.est$p_Med)


par(pty="s")

plot(strat.int.admin2.res$Combined.est$p_Med,strat.int.admin2.res.CI$Combined.est$p_Med)
abline(a=0,b=1)


mean(admin2.direct.est$p_Upp-admin2.direct.est$p_Low)
mean(strat.int.admin2.res$Combined.est$p_Upp-strat.int.admin2.res$Combined.est$p_Low)
mean(strat.int.admin2.res.CI$Combined.est$p_Upp-strat.int.admin2.res.CI$Combined.est$p_Low)


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
                                                    myData  = hiv_dat,
                                                    nSamp = 1000)

# comparison
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=area.iid.nonstrat.admin2.res$p_Med,
     xlim=c(0,0.15),ylim=c(0,0.15))
abline(a=0,b=1)

mean(area.iid.nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med)


setwd(paste0(res_dir,country,'/HIV',sep=''))
#save.image(file = 'HIV_female_15_49.RData')



################################################################
#########  admin-2 random effects cluster level stratified, no int
################################################################


### stratified no interaction
area.iid.strat.admin2 = FitAdmin2IIDStrat(myData = hiv_dat,
                                          clustPrior = iidPrior)

area.iid.strat.admin2.res = aggAdmin2IIDStrat(res.inla = area.iid.strat.admin2,
                                              myData = hiv_dat,
                                              nSamp = 1000,
                                              urb.frac = f15_49.adm2.frac,
                                              admin.ref = admin2.names)

area.iid.strat.admin2.res.CI = aggAdmin2IIDStrat(res.inla = area.iid.strat.admin2,
                                              myData = hiv_dat,
                                              nSamp = 1000,
                                              fixed.urb = F,
                                              urb.frac = adm2.urb.samp,
                                              admin.ref=admin2.names)


par(pty="s")

plot(x=admin2.direct.est$p_Med,y=area.iid.strat.admin2.res$Combined.est$p_Med,
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


################################################################
#########  admin-2 random effects cluster level stratified,  int
################################################################


area.iid.strat.int.admin2 = FitAdmin2IIDStratInt(myData = hiv_dat,
                                                 clustPrior = iidPrior)


summary(area.iid.strat.int.admin2)

area.iid.strat.int.admin2.res = aggAdmin2IIDStratInt(res.inla = area.iid.strat.int.admin2,
                                                     myData = hiv_dat,
                                                     nSamp = 1000,
                                                     urb.frac = f15_49.adm2.frac,
                                                     admin.ref = admin2.names)

area.iid.strat.int.admin2.res.CI = aggAdmin2IIDStratInt(res.inla = area.iid.strat.int.admin2,
                                                     myData = hiv_dat,
                                                     nSamp = 1000,
                                                     fixed.urb = F,
                                                     urb.frac = adm2.urb.samp,
                                                     admin.ref=admin2.names)



par(pty="s")

plot(x=admin2.direct.est$p_Med,y=area.iid.strat.int.admin2.res$Combined.est$p_Med)
abline(a=0,b=1)

mean(area.iid.strat.int.admin2.res$Combined.est$p_Med-admin2.direct.est$p_Med)

mean(admin2.direct.est$p_Upp-admin2.direct.est$p_Low)
mean(strat.int.admin2.res$Combined.est$p_Upp-strat.int.admin2.res$Combined.est$p_Low)
mean(strat.int.admin2.res.CI$Combined.est$p_Upp-strat.int.admin2.res.CI$Combined.est$p_Low)
mean(area.iid.strat.int.admin2.res.CI$Combined.est$p_Upp-area.iid.strat.int.admin2.res.CI$Combined.est$p_Low)


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
                                         myData=hiv_dat, nSamp = 1000,
                                         urb.frac=NULL,strata=FALSE,
                                         admin.ref=admin2.names)

# compare with observed 
par(pty="s")

plot(x=admin2.direct.est$p_Med,y=BYM2.nonstrat.admin2.res$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))

abline(a=0,b=1)

mean(BYM2.nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med)

mean(abs(BYM2.nonstrat.admin2.res$p_Med-admin2.direct.est$p_Med))



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


BYM2.strat.noint.admin2.res <- aggBYMAdmin2(res.inla=BYM2.strat.noint.admin2, 
                                            myData=hiv_dat, nSamp = 1000,
                                            urb.frac=f15_49.adm2.frac,int=FALSE,
                                            fixed.urb = T,
                                            admin.ref=admin2.names)


BYM2.strat.noint.admin2.res.CI <- aggBYMAdmin2(res.inla=BYM2.strat.noint.admin2, 
                                            myData=hiv_dat, nSamp = 1000,
                                            urb.frac=adm2.urb.samp,int=FALSE,
                                            fixed.urb = F,
                                            admin.ref=admin2.names)

par(pty="s")

plot(x=admin2.direct.est$p_Med,y=BYM2.strat.noint.admin2.res.CI$Combined.est$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))
abline(a=0,b=1)


mean(BYM2.strat.noint.admin2.res.CI$Combined.est$p_Upp-BYM2.strat.noint.admin2.res.CI$Combined.est$p_Low)
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

BYM2.strat.int.admin2.res <- aggBYMAdmin2(res.inla=BYM2.strat.int.admin2, 
                                          myData=hiv_dat, nSamp = 1000,
                                          urb.frac=f15_49.adm2.frac,int=TRUE)
# Calculate estimates


par(pty="s")

plot(x=admin2.direct.est$p_Med,y=BYM2.strat.int.admin2.res$Combined.est$p_Med,
     xlim=c(0,0.275),ylim=c(0,0.275))
abline(a=0,b=1)


mean(BYM2.strat.int.admin2.res$Combined.est$p_Med-admin2.direct.est$p_Med)



################################################################
#########  save results
################################################################

setwd(paste0(res_dir,country,'/HIV',sep=''))
load(file = 'HIV_female_15_49.RData')
#save.image(file = 'HIV_female_15_49.RData')