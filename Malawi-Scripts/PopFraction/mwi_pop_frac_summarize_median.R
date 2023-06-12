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
svy.year <- 2015
frame.year<-2008
country.abbrev<-'mwi' #for population rasters
gadm.abbrev = "MWI"


data_dir<-paste0(home_dir,'/Data/')
res_dir<-paste0(home_dir,'/Results/')
pop_dir<-paste0(main_dir,'/Population/')
cluster_dir<-paste0(main_dir,'countryDataFolders/',country)



# national urban fraction for threshold

natl_urban_frac <- 2003309/13077160
cutoff <- natl_urban_frac


# create directories to store results
setwd(paste0(res_dir,country,sep=''))

if(!dir.exists(paths = paste0('UR'))){
  dir.create(path = paste0('UR'))
}


if(!dir.exists(paths = paste0('UR/Fractions'))){
  dir.create(path = paste0('UR/Fractions'))
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
######### load classifications
################################################################

score.method <- 'median'
admin.level <- 'admin2'
setwd(paste0(res_dir,country,'/UR/Fractions',sep=''))

class.method <- 'BART_uncrc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_natl_bench_total_',admin.level,'_frac.rda'))
bart.uncrc.frac <- frac.adm$post.median
#bart.uncrc.frac$Frac<- apply(frac.adm$post.sample, 2, mean, na.rm=T)

class.method <- 'BART_crc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_natl_bench_total_',admin.level,'_frac.rda'))
bart.crc.frac <- frac.adm$post.median
#bart.crc.frac$Frac<- apply(frac.adm$post.sample, 2, mean, na.rm=T)

class.method <- 'BART_mcmc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_natl_bench_total_',admin.level,'_frac.rda'))
bart.mcmc.frac <- frac.adm$post.median
#bart.mcmc.frac$Frac<- apply(frac.adm$post.sample, 2, mean, na.rm=T)



class.method <- 'pop'
load(file=paste0(class.method,'/',class.method,'_natl_thresh_total_',admin.level,'_frac.rda'))
pop.thresh.frac <- frac.adm

class.method <- 'logistic_uncrc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_natl_bench_total_',admin.level,'_frac.rda'))
logis.uncrc.frac <- frac.adm$post.median
#logis.uncrc.frac$Frac<- apply(frac.adm$post.sample, 2, mean, na.rm=T)

class.method <- 'logistic_crc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_natl_bench_total_',admin.level,'_frac.rda'))
logis.crc.frac <- frac.adm$post.median
#logis.crc.frac$Frac<- apply(frac.adm$post.sample, 2, mean, na.rm=T)

class.method <- 'GBT_uncrc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_natl_bench_total_',admin.level,'_frac.rda'))
gbt.uncrc.frac <- frac.adm$post.median
#gbt.uncrc.frac <- gbt.uncrc.frac
#  rename("urb_frac" = "Frac")

class.method <- 'GBT_crc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_natl_bench_total_',admin.level,'_frac.rda'))
gbt.crc.frac <- frac.adm$post.median



################################################################
######### load calibration parameter
################################################################

score.method <- 'median'
admin.level <- 'admin2'
setwd(paste0(res_dir,country,'/UR/Classification/Calibration',sep=''))

class.method <- 'BART_uncrc'
bart.uncrc.tau <- readRDS(file=paste0(class.method,'_',score.method,'_adm_shift.rds'))
bart.uncrc.tau$method <- 'BART-uncorrected'

class.method <- 'BART_crc'
bart.crc.tau <- readRDS(file=paste0(class.method,'_',score.method,'_adm_shift.rds'))
bart.crc.tau$method <-  'BART-corrected'

class.method <- 'BART_mcmc'
bart.mcmc.tau <- readRDS(file=paste0(class.method,'_',score.method,'_adm_shift.rds'))
bart.mcmc.tau$method <-  'BART-MCMC'

class.method <- 'logistic_uncrc'
logis.uncrc.tau <- readRDS(file=paste0(class.method,'_',score.method,'_adm_shift.rds'))
logis.uncrc.tau$method <- 'Logistic Regression-uncorrected'

class.method <- 'logistic_crc'
logis.crc.tau <- readRDS(file=paste0(class.method,'_',score.method,'_adm_shift.rds'))
logis.crc.tau$method <- 'Logistic Regression-corrected'

class.method <- 'GBT_crc'
gbt.crc.tau <- readRDS(file=paste0(class.method,'_',score.method,'_adm_shift.rds'))
gbt.crc.tau$method <-  'GBT-corrected'

class.method <- 'GBT_uncrc'
gbt.uncrc.tau <- readRDS(file=paste0(class.method,'_',score.method,'_adm_shift.rds'))
gbt.uncrc.tau$method <-  'GBT-uncorrected'


################################################################
######### prepare fractions results 
################################################################


setwd(paste0(data_dir,country,sep=''))

load('admin2_urb_ref_tab.rda')



# compare admin2 fractions
match_id <- match(ref.tab$Internal,admin2.names$Internal)
ref.tab <-  ref.tab[match_id,]
ref.tab$method <- 'Census'
ref.tab <- ref.tab %>% rename("known_frac" = "urb_frac")



# prepare truth from census
ref.tab.census <- ref.tab
ref.tab.census$Frac <- ref.tab.census$known_frac

# prepare BART uncorrected 
ref.bart.uncrc <- ref.tab
ref.bart.uncrc <- 
  merge(ref.bart.uncrc,bart.uncrc.frac,by.x='Internal',
        by.y='admin')
ref.bart.uncrc$method <- 'BART-uncorrected'

# prepare BART corrected 
ref.bart.crc <- ref.tab
ref.bart.crc <- 
  merge(ref.bart.crc,bart.crc.frac,by.x='Internal',
        by.y='admin')
ref.bart.crc$method <- 'BART-corrected'

# prepare BART MCMC 
ref.bart.mcmc <- ref.tab
ref.bart.mcmc <- 
  merge(ref.bart.mcmc,bart.mcmc.frac,by.x='Internal',
        by.y='admin')
ref.bart.mcmc$method <- 'BART-MCMC'

# prepare population thresholding
ref.pop.thresh <- ref.tab
ref.pop.thresh <- merge(ref.pop.thresh,pop.thresh.frac[,c('adm.char','Frac')],by.x='Internal',
      by.y='adm.char')
ref.pop.thresh$method <- 'Population Density'

# prepare logistic regression, uncorrected
ref.logis.uncrc <- ref.tab
#ref.logis.uncrc$urb_frac <- logis.uncrc.frac$Frac
ref.logis.uncrc <- merge(ref.logis.uncrc,logis.uncrc.frac,by.x='Internal',
                           by.y='admin')

ref.logis.uncrc$method <- 'Logistic Regression-uncorrected'

# prepare logistic regression, corrected
ref.logis.crc <- ref.tab
#ref.logis.uncrc$urb_frac <- logis.uncrc.frac$Frac
ref.logis.crc <- merge(ref.logis.crc,logis.crc.frac,by.x='Internal',
                         by.y='admin')

ref.logis.crc$method <- 'Logistic Regression-corrected'

# prepare GBT, uncorrected
ref.gbt.uncrc <- ref.tab
ref.gbt.uncrc <- 
  merge(ref.gbt.uncrc,gbt.uncrc.frac,by.x='Internal',
        by.y='admin')
ref.gbt.uncrc$method <- 'GBT-uncorrected'

# prepare GBT, corrected
ref.gbt.crc <- ref.tab
ref.gbt.crc <- 
  merge(ref.gbt.crc,gbt.crc.frac,by.x='Internal',
        by.y='admin')
ref.gbt.crc$method <- 'GBT-corrected'

### make comparsion dataframe

compare.tab <- rbind(ref.tab.census,ref.bart.uncrc,ref.bart.crc,
                     ref.bart.mcmc,
                     ref.pop.thresh,ref.logis.uncrc,
                     ref.gbt.uncrc)


ref.bart.uncrc$diff <- ref.bart.uncrc$Frac - ref.bart.uncrc$known_frac
ref.bart.crc$diff <- ref.bart.crc$Frac - ref.bart.crc$known_frac
ref.bart.mcmc$diff <- ref.bart.mcmc$Frac - ref.bart.mcmc$known_frac
ref.pop.thresh$diff <- ref.pop.thresh$Frac - ref.pop.thresh$known_frac
ref.logis.uncrc$diff <- ref.logis.uncrc$Frac - ref.logis.uncrc$known_frac
ref.logis.crc$diff <- ref.logis.crc$Frac - ref.logis.crc$known_frac
ref.gbt.uncrc$diff <- ref.gbt.uncrc$Frac - ref.gbt.uncrc$known_frac
ref.gbt.crc$diff <- ref.gbt.crc$Frac - ref.gbt.crc$known_frac

compare.diff.tab <- rbind(ref.bart.crc,ref.bart.uncrc,
                          ref.bart.mcmc,
                          ref.pop.thresh,ref.logis.uncrc,
                          ref.logis.crc,
                          ref.gbt.uncrc,ref.gbt.crc
                          )

compare.diff.tab$abs_diff <- abs(compare.diff.tab$diff)

compare.diff.tab <- compare.diff.tab[compare.diff.tab$matched_name!='Likoma',]


################################################################
######### prepare tau results 
################################################################

tau.method.frame <- rbind(gbt.uncrc.tau,gbt.crc.tau,logis.crc.tau,
                          logis.uncrc.tau,bart.uncrc.tau,
                          bart.crc.tau,bart.mcmc.tau)

tau.method.frame <- merge(tau.method.frame,ref.tab[,c('Internal',
                                                      'known_frac','matched_name')],
                          by.x='admin',by.y='Internal')

# exclude Likoma
tau.method.frame <- tau.method.frame[tau.method.frame$admin!='admin2_10',]
  
sort_ref <- ref.tab[order(ref.tab$known_frac),]


tau.method.frame$matched_name <- factor(tau.method.frame$matched_name, levels=sort_ref$matched_name)


################################################################
######### summary fractions results
################################################################

method_vec <- unique(sort(compare.diff.tab$method))
mean_err <- compare.diff.tab %>%
  group_by(method) %>%
  summarise_at(vars(diff), funs(mean(., na.rm=TRUE)))
mean_abs_err <- compare.diff.tab %>%
  group_by(method) %>%
  summarise_at(vars(abs_diff), funs(mean(., na.rm=TRUE)))

mean_err <- paste0(round(as.numeric(mean_err$diff)*100,digits=2),'%')
mean_abs_err <- paste0(round(as.numeric(mean_abs_err$abs_diff)*100,digits=2),'%')


summary_res <- data.frame(method=method_vec,
                          mean_err =mean_err,
                          mean_abs_err= mean_abs_err)
library(xtable)
xtable(summary_res)



################################################################
######### Calibration parameter plot
################################################################

setwd(paste0(res_dir,country,'/Visualization/Maps',sep=''))

pdf('adm2_tau_compare.pdf', height = 8, 
    width = 10)

{
ggplot(tau.method.frame, aes(x=logit(known_frac), y=tau, shape=method, color=method)) +
  scale_shape_manual(values=c(0,1,15:19)) +
  geom_point()+#theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
  xlab('Urban fraction (logit)')+ylab('tau')+
  scale_colour_viridis_d()+
  theme(legend.position="bottom")

}

dev.off()

pdf('adm2_tau_compare_admin_axis.pdf', height = 8, 
    width = 10)

{
  ggplot(tau.method.frame, aes(x=matched_name, y=tau, shape=method, color=method)) +
    scale_shape_manual(values=c(0,1,15:19)) +
    geom_point()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
    xlab('Admin-2 region, ordererd by increasing urban fractions')+ylab('tau')+
    scale_colour_viridis_d()+
    theme(legend.position="bottom")
  
}

dev.off()

# Map
if(FALSE){
tau.method.frame.t <- tau.method.frame[tau.method.frame$method=='GBT-uncorrected',]
ymax <- max(tau.method.frame.t$tau)
ymin <-  min(tau.method.frame.t$tau)

g.1 <- mapPlot(tau.method.frame.t, geo = poly.adm2, by.data = "matched_name", by.geo = "NAME_2", variable = "tau", 
                     is.long=FALSE, per1000=FALSE, removetab=TRUE, legend.label = "tau", direction = -1, 
                     size= 0.1,ylim= c(ymin,ymax))

ymax <- max(logit(tau.method.frame.t$known_frac))+0.1
ymin <-  min(logit(tau.method.frame.t$known_frac))-0.1
tau.method.frame.t$logit_frac <- logit(tau.method.frame.t$known_frac)
g.2 <- mapPlot(tau.method.frame.t, geo = poly.adm2, by.data = "matched_name", by.geo = "NAME_2", variable = "logit_frac", 
               is.long=FALSE, per1000=FALSE, removetab=TRUE, legend.label = "urban fraction", direction = -1, 
               size= 0.1,ylim= c(ymin,ymax))

}

################################################################
######### Plots
################################################################


# sort by urban fraction
sort_ref <- ref.tab[order(ref.tab$known_frac),]


compare.tab$matched_name <- factor(compare.tab$matched_name, levels=sort_ref$matched_name)
compare.diff.tab$matched_name <- factor(compare.diff.tab$matched_name, levels=sort_ref$matched_name)

# remove Likoma
compare.tab <- compare.tab[compare.tab$matched_name!='Likoma',]
compare.tab[compare.tab$method=='Population Density',]$method <- 'Population Thresholding'
compare.diff.tab[compare.diff.tab$method=='Population Density',]$method <- 'Population Thresholding'


setwd(paste0(res_dir,country,'/Visualization/Maps',sep=''))

if(FALSE){
pdf('natl_thresh_adm2_compare.pdf', height = 6, 
    width = 10)

{
  
# plot
ggplot(compare.tab, aes(x=matched_name, y=urb_frac, shape=method, color=method)) +
  scale_shape_manual(values=1:7) +
  geom_point()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
    xlab('Admin-2 region')+ylab('Urban fractions (total population)')
  

}

dev.off()
}


setwd(paste0(res_dir,country,'/Visualization/Maps',sep=''))

pdf('natl_thresh_adm2_compare_err.pdf', height = 8, 
    width = 10)

{
  
  # plot
  ggplot(compare.diff.tab, aes(x=matched_name, y=diff, shape=method, color=method)) +
    scale_shape_manual(values=c(0,1,2,15:19)) +
    geom_point()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
    xlab('Admin-2 region')+ylab('Difference in urban fractions compared with census')+
    geom_hline(yintercept=0, linetype="dashed")+
    scale_colour_viridis_d()+
    theme(legend.position="bottom")
  
  
  
  
}

dev.off()



setwd(paste0(res_dir,country,'/Visualization/Maps',sep=''))

pdf('natl_thresh_adm2_compare_abs_err.pdf', height = 8, 
    width = 10)

{
  
  # plot
  ggplot(compare.diff.tab, aes(x=matched_name, y=abs_diff, shape=method, color=method)) +
    scale_shape_manual(values=c(0,1,2,15:19)) +
    geom_point()+theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+
    xlab('Admin-2 region')+ylab('Absolute difference in urban fractions compared with census')+
    scale_colour_viridis_d()+
    theme(legend.position="bottom")

  
  
  
}

dev.off()



