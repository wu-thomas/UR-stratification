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
library(xtable)

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


## load point estimates of urban/rural fractions
method= 'BART_MCMC'
loss.func= 'MSE'


setwd(paste0(res_dir,country,'/HIV/rev_res',sep=''))
load(file = paste0('HIV_female_15_49_',method,'_',loss.func,'.RData'))


setwd(paste0(res_dir,country,'/HIV',sep=''))

if(!dir.exists(paths = paste0('Visualization'))){
  dir.create(path = paste0('Visualization'))
}


################################################################
#########   load functions for measurements
################################################################

cal_relative <- function(est,ref,digits=2,percent=T,abs=FALSE){
  relative_bias <- (est-ref)/ref
  if(percent){
    scale=100
  }else{scale=1}
  if(abs){
    return(paste0(round(mean(abs(relative_bias))*scale,digits=digits),'%'))
  }else{
    return(paste0(round(mean(relative_bias)*scale,digits=digits),'%'))
    
  }
}

cal_bias <- function(est,ref,digits=2,percent=T,abs=FALSE){
  bias <- est-ref
  if(percent){
    scale=100
  }else{scale=1}
  if(abs){
    return(round(mean(abs(bias))*scale,digits=digits))
  }else{
    return(round(mean(bias)*scale,digits=digits))
    
  }
}  

################################################################
#########   prepare table
################################################################


## prepare results data frame

DIC_vec <- vector()
bias_vec <- vector()
abs_bias_vec <- vector()
relative_bias_vec <- vector()
abs_relative_bias_vec <- vector()
CI_width_vec <- vector()

res.list <- list(nonstrat.admin2.res$Combined.est,strat.admin2.res.A$Combined.est,strat.int.admin2.res.A$Combined.est,
                 area.iid.nonstrat.admin2.res$Combined.est,area.iid.strat.admin2.res.A$Combined.est,
                 area.iid.strat.int.admin2.res.A$Combined.est,BYM2.nonstrat.admin2.res$Combined.est,
                 BYM2.strat.noint.admin2.res.A$Combined.est,BYM2.strat.int.admin2.res.A$Combined.est)

res.samp.list <- list(nonstrat.admin2.res$Combined.draws,strat.admin2.res.A$Combined.draws,strat.int.admin2.res.A$Combined.draws,
                 area.iid.nonstrat.admin2.res$Combined.draws,area.iid.strat.admin2.res.A$Combined.draws,
                 area.iid.strat.int.admin2.res.A$Combined.draws,BYM2.nonstrat.admin2.res$Combined.draws,
                 BYM2.strat.noint.admin2.res.A$Combined.draws,BYM2.strat.int.admin2.res.A$Combined.draws)


if(FALSE){
  res.list <- list(nonstrat.admin2.res,strat.admin2.res.B$Combined.est,strat.int.admin2.res.B$Combined.est,
                   area.iid.nonstrat.admin2.res,area.iid.strat.admin2.res.B$Combined.est,
                   area.iid.strat.int.admin2.res.B$Combined.est,BYM2.nonstrat.admin2.res,
                   BYM2.strat.noint.admin2.res.B$Combined.est,BYM2.strat.int.admin2.res.B$Combined.est)
  
}

for(idx in 1:length(res.list)){

  est_vec <- res.list[[idx]]$p_Med
  
  bias_vec[idx]<-cal_bias(est=est_vec,ref=admin2.direct.est$p_Med,abs=F)
  abs_bias_vec[idx] <- cal_bias(est=est_vec,ref=admin2.direct.est$p_Med,abs=T)
  relative_bias_vec[idx]<- cal_relative(est=est_vec,ref=admin2.direct.est$p_Med,abs=F)
  abs_relative_bias_vec[idx] <- cal_relative(est=est_vec,ref=admin2.direct.est$p_Med,abs=T)
  CI_width_vec[idx] <- mean( res.list[[idx]]$p_Upp- res.list[[idx]]$p_Low)
  
  #res.list[[idx]]$Combined.est$coeff_var <- apply(res.samp.list[[idx]],2,function(x) {sd(x)/mean(x)} )   # not generally a useful result

}


## fixed effects model
DIC_vec[1]<-inla.nonstrat.admin2$dic$dic
DIC_vec[2]<-inla.strat.admin2$dic$dic
DIC_vec[3]<-inla.strat.int.admin2$dic$dic

## area iid random effects model
DIC_vec[4]<-area.iid.nonstrat.admin2$dic$dic
DIC_vec[5]<-area.iid.strat.admin2$dic$dic
DIC_vec[6]<-area.iid.strat.int.admin2$dic$dic

## BYM2 models
DIC_vec[7]<-BYM2.nonstrat.admin2$dic$dic
DIC_vec[8]<-BYM2.strat.noint.admin2$dic$dic
DIC_vec[9]<-BYM2.strat.int.admin2$dic$dic


model_type <- rep(c('Fixed effects','Area IID random effects','BYM2'),each=3)
stratification <- rep(c('Area','U/R+Area','U/R*Area'),times=3)

summary.res <- data.frame (model_type = model_type, 
                           stratification= stratification,
                           #DIC = DIC_vec,
                           bias=bias_vec,
                           abs_bias= abs_bias_vec,
                           relative_bias = relative_bias_vec,
                           abs_relative_bias =abs_relative_bias_vec,
                           CI_width_avg = CI_width_vec*100)


# kable(summary.res)
xtable(summary.res)
################################################################
#########   prepare random urban sample table
################################################################


## prepare results data frame
res.list.A <- list(strat.admin2.res.A$Combined.est,
                   strat.int.admin2.res.A$Combined.est,
                   area.iid.strat.admin2.res.A$Combined.est,
                   area.iid.strat.int.admin2.res.A$Combined.est,
                   BYM2.strat.noint.admin2.res.A$Combined.est,
                   BYM2.strat.int.admin2.res.A$Combined.est)


res.list.B <- list(strat.admin2.res.B$Combined.est,
                 strat.int.admin2.res.B$Combined.est,
                 area.iid.strat.admin2.res.B$Combined.est,
                 area.iid.strat.int.admin2.res.B$Combined.est,
                 BYM2.strat.noint.admin2.res.B$Combined.est,
                 BYM2.strat.int.admin2.res.B$Combined.est)

ci_width_A <- vector()
ci_width_B <- vector()
#ci_width_diff_percent <- vector()

for(idx in 1:length(res.list.B)){
  
  #ci_width_B[idx] <- mean( res.list.B[[idx]]$p_Upp- res.list.B[[idx]]$p_Low)
  #ci_width_A[idx] <- mean( res.list.A[[idx]]$p_Upp- res.list.A[[idx]]$p_Low)
  
  # exclude Likoma
  ci_width_B[idx] <- mean( res.list.B[[idx]]$p_Upp[c(1:9,11:28)]- res.list.B[[idx]]$p_Low[c(1:9,11:28)])
  ci_width_A[idx] <- mean( res.list.A[[idx]]$p_Upp[c(1:9,11:28)]- res.list.A[[idx]]$p_Low[c(1:9,11:28)])
  
  #tmp_CI_A <- res.list.B[[idx]]$p_Upp[c(1:9,11:28)]- res.list.B[[idx]]$p_Low[c(1:9,11:28)]
  #tmp_CI_B <- res.list.A[[idx]]$p_Upp[c(1:9,11:28)]- res.list.A[[idx]]$p_Low[c(1:9,11:28)]
    
  #ci_width_diff_percent[idx] <- mean( (tmp_CI_B-tmp_CI_A)/tmp_CI_A)
  
}


model_type <- rep(c('Fixed effects','Area IID random effects','BYM2'),each=2)
stratification <- rep(c('U/R+Area','U/R*Area'),times=3)

summary.width <- data.frame (model_type = model_type, 
                           stratification= stratification,
                          ci_width_fixed_r= ci_width_A*100,
                            ci_width_samp_r= ci_width_B*100,
                          ci_width_diff = ci_width_B*100-ci_width_A*100)

summary.width$percent_diff <- round((summary.width$ci_width_diff/summary.width$ci_width_fixed_r)*100,digits=2)
  
xtable(summary.width)

################################################################
#########   prepare data for plot
################################################################

### fixed urban fraction plots 
fixed.res.A <- strat.admin2.res.A$Combined.est

fixed.res.A$method <- 'fixed effects, no int'
fixed.int.res.A <- strat.int.admin2.res.A$Combined.est
fixed.int.res.A$method <- 'fixed effects, int'
area.iid.res.A <- area.iid.strat.admin2.res.A$Combined.est
area.iid.res.A$method <-'area iid, no int'
area.iid.int.res.A <- area.iid.strat.int.admin2.res.A$Combined.est
area.iid.int.res.A$method <-'area iid, int'
area.sp.res.A <- BYM2.strat.noint.admin2.res.A$Combined.est
area.sp.res.A$method <-'BYM2, no int'
area.sp.int.res.A <- BYM2.strat.int.admin2.res.A$Combined.est
area.sp.int.res.A$method <-'BYM2, int'

###############################################################################
#########   admin2 direct vs. stratitifed and nonstratified fixed effects
###############################################################################

adm2.comp <- admin2.direct$overall
adm2.comp <- adm2.comp[,c('admin2','p_Med')]

adm2.comp.nonstrat <- merge(adm2.comp,nonstrat.admin2.res$Combined.est[,c('admin.name','p_Med')], 
                            by.x=c('admin2'),
                            by.y=c('admin.name'))

adm2.comp.nonstrat <- adm2.comp.nonstrat[,c('admin2','p_Med.x','p_Med.y')]
colnames(adm2.comp.nonstrat) <- c('region','direct_est','clust_est')
adm2.comp.nonstrat$method <- 'Without U/R stratification'

adm2.comp.strat <- merge(adm2.comp,fixed.int.res.A [,c('admin.name','p_Med')], 
                         by.x=c('admin2'),
                         by.y=c('admin.name'))

adm2.comp.strat <- adm2.comp.strat[,c('admin2','p_Med.x','p_Med.y')]
colnames(adm2.comp.strat) <- c('region','direct_est','clust_est')
adm2.comp.strat$method <- 'With U/R stratification'


adm2.comp.toplot <- rbind(adm2.comp.strat,adm2.comp.nonstrat)
adm2.comp.toplot$method <- factor(adm2.comp.toplot$method,levels=c('Without U/R stratification',
                                                                   'With U/R stratification'))



range.toplot <- c(min(c(adm2.comp.toplot$direct_est,adm2.comp.toplot$clust_est)),
                  max(c(adm2.comp.toplot$direct_est,adm2.comp.toplot$clust_est)))


adm2.comp.scatter.g <- ggplot(data = adm2.comp.toplot)+geom_point(aes(x = direct_est, y = clust_est),
                                                                 color=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.7))+
  facet_wrap(vars(method))+ coord_fixed() +
  xlim(range.toplot)+ylim(range.toplot)+
  geom_abline(slope=1, intercept=0)+
  xlab('Direct estimates')+ylab('Cluster-level model estimates')+
  theme_bw()+
  theme ( strip.text.x = element_text(size = 15),
          axis.title  = element_text(size = 15),
          axis.text  = element_text(size = 12) )+
  scale_x_continuous(breaks = c(1:5)*0.05)+
  scale_y_continuous(breaks = c(1:5)*0.05)


setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(adm2.comp.scatter.g, width=10, height = 8, file = paste0("compare-admin2-est-scatter.pdf"))

################################################################
#########   admin2 estimates compare
################################################################

if(FALSE){
### prepare plot object
#to_plot <- rbind(fixed.res.A,fixed.int.res.A,area.iid.res.A,area.iid.int.res.A,area.sp.res.A,area.sp.int.res.A)
to_plot <- rbind(fixed.int.res.A,area.iid.int.res.A,area.sp.int.res.A)


colnames(to_plot) <- c('region','admin_name','lower','median','upper','urb','method')
to_plot <- to_plot[,c('region','lower','upper','median','method')]
to_plot$years.num <- NA
class(to_plot) <- c("SUMMERproj", "data.frame")

my.dodge <- position_dodge(width = 1)

ggplot(data = to_plot)+geom_point(aes(x = region, y = median,color=method), position = my.dodge) + 
  xlab("Region") + ylab("")+
  geom_errorbar(aes(x = region,ymin = lower, ymax = upper,color=method),
                position = my.dodge, size = .5, width = .05, alpha = 0.5)+
  labs(color='method') 


plot(to_plot)+geom_errorbar(aes(ymin = lower, ymax = upper), size = .5, width = .05, alpha = 0.5)

plot(to_plot,plot.CI = TRUE, dodge.width = 0.25, per1000=F) +
  theme(legend.position = 'bottom') + scale_linetype(guide='none') + ggtitle(i) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = c(beg.year:end.year))+
  labs(color='Method') 

}


################################################################
#########   Maps: Admin-2, stratified estimates
################################################################

# direct estimates
adm2.direct.tp <- admin2.direct$overall
adm2.direct.tp$admin.name <- adm2.direct.tp$admin2

# fixed effects model
adm2.fixed.tp <- fixed.int.res.A

# BYM2 model
adm2.strat.int.tp <- BYM2.strat.int.admin2.res.A$Combined.est

#class(adm2.strat.noint.tp) <- c("SUMMERproj", "data.frame")
adm2.direct.tp$method<- 'Survey direct estimates'
adm2.fixed.tp$method<- 'Fixed effects model'
adm2.strat.int.tp$method <- 'BYM2 model'

adm2.map.tp <-rbind(adm2.direct.tp[,c('p_Med','method','admin.name')],
                    adm2.fixed.tp[,c('p_Med','method','admin.name')],
                    adm2.strat.int.tp[,c('p_Med','method','admin.name')])

adm2.map.tp$method <- factor(adm2.map.tp$method,levels=c('Survey direct estimates', 'Fixed effects model',
                                                            'BYM2 model'))

ymax <- ceiling(  max(adm2.map.tp$p_Med)*100 )/100
ymin <- floor(  min(adm2.map.tp$p_Med)*100)/100


g.compare <- mapPlot(data = adm2.map.tp, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
        values = c("p_Med"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
        is.long = TRUE, ncol = 3)


if(FALSE){
  # individual plots
  g1a <- mapPlot(adm2.direct.tp, geo = poly.adm2, by.data = "admin2", by.geo = "NAME_2", variable = "p_Med", 
                 is.long=FALSE, per1000=FALSE, removetab=TRUE, legend.label = "HIV", direction = -1, 
                 size= 0.1,ylim= c(ymin,ymax))
  
  
  g1a <- mapPlot(adm2.strat.noint.tp, geo = poly.adm2, by.data = "admin.name", by.geo = "NAME_2", variable = "p_Med", 
                 is.long=FALSE, per1000=FALSE, removetab=TRUE, legend.label = "HIV", direction = -1, 
                 size= 0.1)
}


g.compare <- g.compare+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='HIV prevalence',
                                label.position = "bottom"))
  

setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(g.compare, width=8, height = 10, file = paste0("compare-admin2-strat-est.pdf"))




################################################################
#########   Maps: Admin-2, unstratified estimates
################################################################

# direct estimates
adm2.direct.tp <- admin2.direct$overall
adm2.direct.tp$admin.name <- adm2.direct.tp$admin2

# fixed effects model
adm2.nonstrat.fixed.tp <- nonstrat.admin2.res
adm2.strat.fixed.tp <- strat.admin2.res.A$Combined.est

# BYM2 model
adm2.nonstrat.tp <- BYM2.nonstrat.admin2.res

#class(adm2.strat.noint.tp) <- c("SUMMERproj", "data.frame")
adm2.direct.tp$method<- 'Survey direct estimates'
adm2.nonstrat.fixed.tp$method<- 'Fixed effects without U/R stratification'
adm2.strat.fixed.tp$method <-  'Fixed effects with U/R stratification'

adm2.nonstrat.tp$method <- 'BYM2 model'

adm2.map.fixed.tp <-rbind(adm2.strat.fixed.tp[,c('p_Med','method','admin.name')],
                    adm2.nonstrat.fixed.tp[,c('p_Med','method','admin.name')])

adm2.map.fixed.tp$method <- factor(adm2.map.fixed.tp$method,levels=c('Fixed effects without U/R stratification',
                                                         'Fixed effects with U/R stratification'))

ymax <- ceiling(  max(adm2.map.tp$p_Med)*100 )/100
ymin <- floor(  min(adm2.map.tp$p_Med)*100)/100


g.compare <- mapPlot(data = adm2.map.fixed.tp, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
                     values = c("p_Med"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
                     is.long = TRUE, ncol = 3)


if(FALSE){
  # individual plots
  g1a <- mapPlot(adm2.direct.tp, geo = poly.adm2, by.data = "admin2", by.geo = "NAME_2", variable = "p_Med", 
                 is.long=FALSE, per1000=FALSE, removetab=TRUE, legend.label = "HIV", direction = -1, 
                 size= 0.1,ylim= c(ymin,ymax))
  
  
  g1a <- mapPlot(adm2.strat.noint.tp, geo = poly.adm2, by.data = "admin.name", by.geo = "NAME_2", variable = "p_Med", 
                 is.long=FALSE, per1000=FALSE, removetab=TRUE, legend.label = "HIV", direction = -1, 
                 size= 0.1)
}


g.compare <- g.compare+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='HIV prevalence',
                                label.position = "bottom"))


setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(g.compare, width=8, height = 10, file = paste0("compare-admin2-fixed-est.pdf"))



################################################################
#########   Maps: Admin-2, coefficient of variation
################################################################

# direct estimates
adm2.direct.tp <- admin2.direct.aggre.A$Combined.est
adm2.direct.tp$coeff_var <- apply(admin2.direct.aggre.A$Combined.draws,2,function(x) {sd(x)/mean(x)} )   


# fixed effects model
adm2.fixed.nonstrat.tp <- nonstrat.admin2.res$Combined.est
adm2.fixed.nonstrat.tp$coeff_var <- apply(nonstrat.admin2.res$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.fixed.strat.noint.tp <- strat.admin2.res.A$Combined.est
adm2.fixed.strat.noint.tp$coeff_var <- apply(strat.admin2.res.A$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.fixed.strat.int.tp <- strat.int.admin2.res.A$Combined.est
adm2.fixed.strat.int.tp$coeff_var <- apply(strat.int.admin2.res.A$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.map.cov.fixed.tp <-rbind(cbind(adm2.fixed.nonstrat.tp[,c('coeff_var','admin.name')],method='Area-only'),
                              cbind(adm2.fixed.strat.noint.tp[,c('coeff_var','admin.name')],method='Area + U/R'),
                              cbind(adm2.fixed.strat.int.tp[,c('coeff_var','admin.name')],method='Area * U/R'))

adm2.map.cov.fixed.tp$coeff_var_percent <- adm2.map.cov.fixed.tp$coeff_var*100

ymax <- max(adm2.map.cov.fixed.tp$coeff_var_percent)
ymin <- min(adm2.map.cov.fixed.tp$coeff_var_percent)
adm2.map.cov.fixed.tp$method <- factor(adm2.map.cov.fixed.tp$method,levels=c('Area-only',
                                                                   'Area + U/R','Area * U/R'))

g.compare.fixed.cov <- mapPlot(data = adm2.map.cov.fixed.tp, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
                         values = c("coeff_var_percent"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
                         is.long = TRUE, ncol = 3)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2.5,'cm'),legend.title = element_text(size=15),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='Coefficient of variation for HIV prevalence (%)',
                                label.position = "bottom"))

ggsave(g.compare.fixed.cov, width=8, height = 10, file = paste0("compare-admin2-fixed-cov-method.pdf"))

# area IID model
adm2.area.nonstrat.tp <- area.iid.nonstrat.admin2.res$Combined.est
adm2.area.nonstrat.tp$coeff_var <- apply(area.iid.nonstrat.admin2.res$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.area.strat.noint.tp <- area.iid.strat.admin2.res.A$Combined.est
adm2.area.strat.noint.tp$coeff_var <- apply(area.iid.strat.admin2.res.A$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.area.strat.int.tp <- area.iid.strat.int.admin2.res.A$Combined.est
adm2.area.strat.int.tp$coeff_var <- apply(area.iid.strat.int.admin2.res.A$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.map.cov.area.tp <-rbind(cbind(adm2.area.nonstrat.tp[,c('coeff_var','admin.name')],method='Area-only'),
                              cbind(adm2.area.strat.noint.tp[,c('coeff_var','admin.name')],method='Area + U/R'),
                              cbind(adm2.area.strat.int.tp[,c('coeff_var','admin.name')],method='Area * U/R'))

adm2.map.cov.area.tp$coeff_var_percent <- adm2.map.cov.area.tp$coeff_var*100

ymax <- max(adm2.map.cov.area.tp$coeff_var_percent)
ymin <- min(adm2.map.cov.area.tp$coeff_var_percent)

adm2.map.cov.area.tp$method <- factor(adm2.map.cov.area.tp$method,levels=c('Area-only',
                                                                             'Area + U/R','Area * U/R'))

g.compare.area.cov <- mapPlot(data = adm2.map.cov.area.tp, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
                               values = c("coeff_var_percent"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
                               is.long = TRUE, ncol = 3)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2.5,'cm'),legend.title = element_text(size=15),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='Coefficient of variation for HIV prevalence (%)',
                                label.position = "bottom"))


ggsave(g.compare.area.cov, width=8, height = 10, file = paste0("compare-admin2-area-cov-method.pdf"))

# BYM2 model
adm2.BYM2.nonstrat.tp <- BYM2.nonstrat.admin2.res$Combined.est
adm2.BYM2.nonstrat.tp$coeff_var <- apply(BYM2.nonstrat.admin2.res$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.BYM2.strat.noint.tp <- BYM2.strat.noint.admin2.res.A$Combined.est
adm2.BYM2.strat.noint.tp$coeff_var <- apply(BYM2.strat.noint.admin2.res.A$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.BYM2.strat.int.tp <- BYM2.strat.int.admin2.res.A$Combined.est
adm2.BYM2.strat.int.tp$coeff_var <- apply(BYM2.strat.int.admin2.res.A$Combined.draws,2,function(x) {sd(x)/mean(x)} )   

adm2.map.cov.BYM2.tp <-rbind(cbind(adm2.BYM2.nonstrat.tp[,c('coeff_var','admin.name')],method='Area-only'),
                              cbind(adm2.BYM2.strat.noint.tp[,c('coeff_var','admin.name')],method='Area + U/R'),
                              cbind(adm2.BYM2.strat.int.tp[,c('coeff_var','admin.name')],method='Area * U/R'))

adm2.map.cov.BYM2.tp$method <- factor(adm2.map.cov.BYM2.tp$method,levels=c('Area-only',
                                                                             'Area + U/R','Area * U/R'))

adm2.map.cov.BYM2.tp$coeff_var_percent <- adm2.map.cov.BYM2.tp$coeff_var*100

ymax <- max(adm2.map.cov.BYM2.tp$coeff_var_percent)
ymin <- min(adm2.map.cov.BYM2.tp$coeff_var_percent)

g.compare.BYM2.cov <- mapPlot(data = adm2.map.cov.BYM2.tp, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
                               values = c("coeff_var_percent"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
                               is.long = TRUE, ncol = 3)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2.5,'cm'),legend.title = element_text(size=15),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='Coefficient of variation for HIV prevalence (%)',
                                label.position = "bottom"))


ggsave(g.compare.BYM2.cov, width=8, height = 10, file = paste0("compare-admin2-BYM2-cov-method.pdf"))


# Comparing stratified models
adm2.map.cov.st.tp <-rbind(cbind(adm2.direct.tp[,c('coeff_var','admin.name')],method='Survey direct estimates'),
                           cbind(adm2.fixed.strat.int.tp[,c('coeff_var','admin.name')],method='Fixed effects model'),
                             cbind(adm2.area.strat.int.tp[,c('coeff_var','admin.name')],method='Area IID RE model'),
                             cbind(adm2.BYM2.strat.int.tp[,c('coeff_var','admin.name')],method= 'BYM2 model'))

adm2.map.cov.st.tp$method <- factor(adm2.map.cov.st.tp$method,levels=c('Survey direct estimates', 'Fixed effects model',
                                                                 'Area IID RE model',
                                                            'BYM2 model'))


adm2.map.cov.st.tp$coeff_var_percent <- adm2.map.cov.st.tp$coeff_var*100

ymax <- max(adm2.map.cov.st.tp$coeff_var_percent)
ymin <- min(adm2.map.cov.st.tp$coeff_var_percent)


g.compare.st.cov <- mapPlot(data = adm2.map.cov.st.tp, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
                        values = c("coeff_var_percent"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
                        is.long = TRUE, ncol = 4)

g.compare.st.cov <- g.compare.st.cov+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2.5,'cm'),legend.title = element_text(size=15),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='Coefficient of variation for HIV prevalence (%)',
                                label.position = "bottom"))


setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(g.compare.st.cov, width=10, height = 10, file = paste0("compare-admin2-cov-method.pdf"))



################################################################
#########   Maps: Admin-2, CI width
################################################################

# direct estimates
adm2.direct.tp <- admin2.direct$overall
adm2.direct.tp$admin.name <- adm2.direct.tp$admin2
adm2.direct.tp$ci_width <- adm2.direct.tp$p_Upp-adm2.direct.tp$p_Low

# fixed effects model
adm2.fixed.tp <- fixed.int.res.A
adm2.fixed.tp$ci_width <- adm2.fixed.tp$p_Upp-adm2.fixed.tp$p_Low

# fixed effects model
adm2.iid.tp <- area.iid.strat.int.admin2.res.A$Combined.est
adm2.iid.tp$ci_width <- adm2.iid.tp$p_Upp-adm2.iid.tp$p_Low

# BYM2 model
adm2.strat.int.tp <- BYM2.strat.int.admin2.res.A$Combined.est
adm2.strat.int.tp$ci_width <- adm2.strat.int.tp$p_Upp-adm2.strat.int.tp$p_Low

#class(adm2.strat.noint.tp) <- c("SUMMERproj", "data.frame")
adm2.direct.tp$method<- 'Survey direct estimates'
adm2.fixed.tp$method<- 'Fixed effects model'
adm2.iid.tp$method <- 'Area IID RE model'
adm2.strat.int.tp$method <- 'BYM2 model'

adm2.map.ci.tp <-rbind(adm2.direct.tp[,c('ci_width','method','admin.name')],
                    adm2.fixed.tp[,c('ci_width','method','admin.name')],
                    adm2.iid.tp[,c('ci_width','method','admin.name')],
                    adm2.strat.int.tp[,c('ci_width','method','admin.name')])

adm2.map.ci.tp$method <- factor(adm2.map.ci.tp$method,levels=c('Survey direct estimates', 'Fixed effects model','Area IID RE model',
                                                         'BYM2 model'))

ymax <- ceiling(  max(adm2.map.ci.tp$ci_width)*100 )/100
ymin <- floor(  min(adm2.map.ci.tp$ci_width)*100)/100


g.compare.ci <- mapPlot(data = adm2.map.ci.tp, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
                     values = c("ci_width"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
                     is.long = TRUE, ncol = 4)

g.compare.ci <- g.compare.ci+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='95% CI width for HIV prevalence',
                                label.position = "bottom"))


setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(g.compare.ci, width=10, height = 10, file = paste0("compare-admin2-CI-method.pdf"))

################################################################
#########   Maps: Admin-2, fixed effects
################################################################

# nonstratified fixed effects model
adm2.nonstrat.tp <- nonstrat.admin2.res$Combined.est
adm2.nonstrat.tp$admin.char <- adm2.nonstrat.tp$admin.ref.Internal

# stratified fixed effects model
adm2.strat.tp <- strat.int.admin2.res.A$Combined.est

#class(adm2.strat.noint.tp) <- c("SUMMERproj", "data.frame")
adm2.nonstrat.tp$method<- 'Fixed Effect without U/R stratification'
adm2.strat.tp$method<- 'Fixed Effect with U/R stratification'


adm2.map.tp <-rbind(adm2.nonstrat.tp[,c('p_Med','method','admin.char')],
                    adm2.strat.tp[,c('p_Med','method','admin.char')])


poly.adm2@data$Internal<-admin2.names$Internal[c(1:dim(admin2.names)[1])]

adm2.map.tp$method <- factor(adm2.map.tp$method,levels=c('Fixed Effect without U/R stratification',
                                                         'Fixed Effect with U/R stratification'))

ymax <- ceiling(  max(adm2.map.tp$p_Med)*100 )/100
ymin <- floor(  min(adm2.map.tp$p_Med)*100)/100


g.compare.FE.adm2 <- mapPlot(data = adm2.map.tp, geo = poly.adm2, variables = c("method"), 
                               size=0.1,ylim= c(ymin,ymax),
                               values = c("p_Med"), by.data = "admin.char", by.geo = "Internal",
                               direction = -1, legend.label = "HIV",
                               is.long = TRUE, ncol = 2)

g.compare.FE.adm2 <- g.compare.FE.adm2+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='HIV prevalence',
                                label.position = "bottom"))


setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(g.compare.FE.adm2, width=8, height = 10, file = paste0("compare-admin2-fixed-effects-map.pdf"))

################################################################
#########   Maps: Admin-2, BYM2
################################################################

# nonstratified BYM2
adm2.nonstrat.tp <- BYM2.nonstrat.admin2.res$Combined.est
adm2.nonstrat.tp$admin.char <- adm2.nonstrat.tp$admin.ref.Internal

# stratified BYM2
adm2.strat.tp <- BYM2.strat.int.admin2.res.A$Combined.est 

#class(adm2.strat.noint.tp) <- c("SUMMERproj", "data.frame")
adm2.nonstrat.tp$method<- 'BYM2 without U/R stratification'
adm2.strat.tp$method<- 'BYM2 with U/R stratification'


adm2.map.tp <-rbind(adm2.nonstrat.tp[,c('p_Med','method','admin.char')],
                    adm2.strat.tp[,c('p_Med','method','admin.char')])


poly.adm2@data$Internal<-admin2.names$Internal[c(1:dim(admin2.names)[1])]

adm2.map.tp$method <- factor(adm2.map.tp$method,levels=c('BYM2 without U/R stratification',
                                                         'BYM2 with U/R stratification'))

ymax <- ceiling(  max(adm2.map.tp$p_Med)*100 )/100
ymin <- floor(  min(adm2.map.tp$p_Med)*100)/100


g.compare.BYM2.adm2 <- mapPlot(data = adm2.map.tp, geo = poly.adm2, variables = c("method"), 
                          size=0.1,ylim= c(ymin,ymax),
                          values = c("p_Med"), by.data = "admin.char", by.geo = "Internal",
                          direction = -1, legend.label = "HIV",
                          is.long = TRUE, ncol = 2)

g.compare.BYM2.adm2 <- g.compare.BYM2.adm2+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='HIV prevalence',
                                label.position = "bottom"))


setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(g.compare.BYM2.adm2, width=8, height = 10, file = paste0("compare-admin2-BYM2-map.pdf"))

################################################################
#########   Maps: Admin-3
################################################################

# nonstratified BYM2
adm3.nonstrat.tp <- BYM2.nonstrat.admin3.res

# stratified BYM2
adm3.strat.tp <- BYM2.strat.admin3.res$Combined.est

#class(adm2.strat.noint.tp) <- c("SUMMERproj", "data.frame")
adm3.nonstrat.tp$method<- 'BYM2 without U/R stratification'
adm3.strat.tp$method<- 'BYM2 with U/R stratification'


adm3.map.tp <-rbind(adm3.nonstrat.tp[,c('p_Med','method','admin.char')],
                    adm3.strat.tp[,c('p_Med','method','admin.char')])


poly.adm3@data$Internal<-admin3.names$Internal[c(1:dim(admin3.names)[1])]

adm3.map.tp$method <- factor(adm3.map.tp$method,levels=c('BYM2 without U/R stratification',
                                                         'BYM2 with U/R stratification'))

ymax <- ceiling(  max(adm3.map.tp$p_Med)*100 )/100
ymin <- floor(  min(adm3.map.tp$p_Med)*100)/100


g.compare.adm3 <- mapPlot(data = adm3.map.tp, geo = poly.adm3, variables = c("method"), 
                     size=0.1,ylim= c(ymin,ymax),
                     values = c("p_Med"), by.data = "admin.char", by.geo = "Internal",
                     direction = -1, legend.label = "HIV",
                     is.long = TRUE, ncol = 2)

g.compare.adm3 <- g.compare.adm3+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='HIV prevalence',
                                label.position = "bottom"))


setwd(paste0(res_dir,country,'/HIV/Visualization',sep=''))

ggsave(g.compare.adm3, width=8, height = 10, file = paste0("compare-admin3-est.pdf"))
