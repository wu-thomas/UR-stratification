################################################################
#########   load libraries
################################################################
#### Libraries ####

rm(list = ls())

library(SUMMER)
#help(package = "SUMMER", help_type = "html")
#utils::browseVignettes(package = "SUMMER")
library(RColorBrewer)
library(dplyr)

library(INLA)
library(ggplot2)
library(maptools)
library(gridExtra)

library(raster)
library(parallel)
library(tictoc)
library(bartMachine)
library(stringdist)
library(reshape)
library(data.table)
library(rgdal)
library(viridis)

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
script_dir<-paste0(home_dir,'/Scripts/')



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


setwd(paste0(script_dir,country,'/PopFraction/',sep=''))

source(paste0(country.abbrev,'_pop_frac_functions.R'))

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

score.method <- 'mse'
admin.level <- 'admin2'
setwd(paste0(res_dir,country,'/UR/Fractions',sep=''))

class.method <- 'BART_uncrc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_adm_bench_total_',admin.level,'_frac.rda'))
bart.uncrc.samp <- as.data.frame(frac.adm$post.sample)

class.method <- 'BART_crc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_adm_bench_total_',admin.level,'_frac.rda'))
bart.crc.samp <-  as.data.frame(frac.adm$post.sample)

class.method <- 'BART_mcmc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_adm_bench_total_',admin.level,'_frac.rda'))
bart.mcmc.samp <- as.data.frame(frac.adm$post.sample)


class.method <- 'logistic_uncrc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_adm_bench_total_',admin.level,'_frac.rda'))
logis.uncrc.samp <- as.data.frame(frac.adm$post.sample)

class.method <- 'logistic_crc'
load(file=paste0(class.method,'/',class.method,'_',score.method,'_adm_bench_total_',admin.level,'_frac.rda'))
logis.crc.samp <-  as.data.frame(frac.adm$post.sample)


################################################################
######### prepare fractions results 
################################################################


setwd(paste0(data_dir,country,sep=''))

load('admin2_urb_ref_tab.rda')



# compare admin2 fractions
match_id <- match(ref.tab$Internal,admin2.names$Internal)
ref.tab <-  ref.tab[match_id,]
ref.tab$method <- 'Census'
#ref.tab <- ref.tab %>% rename("known_frac" = "urb_frac")


setwd(paste0(res_dir,country,'/Visualization/',sep=''))

method.vec <- c('BART_crc','BART_uncrc','BART_MCMC','logistic_crc','logistic_uncrc')
frac.samp.list <- list(bart.crc.samp,bart.uncrc.samp,bart.mcmc.samp,logis.crc.samp,logis.uncrc.samp)

score.method <- 'mse'

################################################################
######### plot uncertainties in admin regions
################################################################

for(i in 1:length(method.vec)){


tmp.samp <- frac.samp.list[[i]]
tmp.samp.long <- melt(setDT(tmp.samp), measure.vars = 1:28, variable.name = "admin2_char")
tmp.samp.long$urb_frac <- tmp.samp.long$value
match.id <- match(tmp.samp.long$admin2_char,ref.tab$Internal)
tmp.samp.long$matched_name <- ref.tab$matched_name[match.id]

adm.samp.q50 <- apply(tmp.samp, 2, median, na.rm=T)
match.id <- match(ref.tab$Internal,names(adm.samp.q50))

### frame for vertical bar
v_bar <- ref.tab
v_bar$admin2_char <- v_bar$Internal
v_bar$samp_q50 <- adm.samp.q50[match.id]


adm.density.plot <- tmp.samp.long[tmp.samp.long$matched_name!='Likoma',] %>% 
  ggplot(aes(x=urb_frac)) +
  geom_density( fill="dodgerblue", alpha=0.5)+
  xlab('urban fraction')+#geom_vline(xintercept=urb_frac, linetype="dashed", size=1,color="darkblue")+
  geom_vline(data = v_bar[v_bar$matched_name!='Likoma',],linetype="dashed", size=0.35,color="darkblue",
             mapping = aes(xintercept = urb_frac)) +
  geom_vline(data = v_bar[v_bar$matched_name!='Likoma',],linetype="twodash", size=0.35,color="red",
             mapping = aes(xintercept = samp_q50)) +
  ggtitle("Admin-2 total population urban fraction, Admin-2 calibrated") +
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~matched_name,scales = "free",
             ncol=3)


pdf(paste0(method.vec[i],'_',score.method,'_adm2_density_adm2_calibrated.pdf'), height = 16, 
    width = 9)
{
  
  print(adm.density.plot)
  
}
dev.off()

}






################################################################
######### load fractions for distribution shift plot
################################################################

######### load national grid 
setwd(paste0(data_dir,country,sep=''))
load("prepared_dat/natl_grid.rda")

urb_dat$complete <- complete.cases(urb_dat) 
urb_dat$pixel_id <- c(1:dim(urb_dat)[1])

complete_id <- complete.cases(urb_dat) 
complete_urb_dat <- urb_dat[complete_id,]

setwd(paste0(data_dir,country,sep=''))

load('admin2_urb_ref_tab.rda')

score.method <- 'mse'
admin.level <- 'admin2'
#setwd(paste0(res_dir,country,'/UR/Fractions',sep=''))

setwd(paste0(res_dir,country,'/BART_jitter/',sep=''))

mcmc.draws.frame <- readRDS('nolikoma_natl_pred_mcmc_draws_1000.rda')
mcmc.draws.frame[mcmc.draws.frame> expit(15)] <-  expit(15)
mcmc.draws.frame <- as.data.frame(mcmc.draws.frame)


complete.natl.dat <- urb_dat[urb_dat$complete,]
pixel.grid <- complete.natl.dat[,c('x','y')]
pop.vec <- complete.natl.dat$pop_den
prob.frame <- as.data.frame(mcmc.draws.frame)
adm.vec <- complete.natl.dat$admin2.char

### get results for admin2_16, Mwanza
adm.prob <- split( prob.frame , f = adm.vec )
adm.pop <- split( pop.vec , f = adm.vec )

tmp.adm.prob <- adm.prob[[21]] #19 # 8
tmp.adm.pop <- adm.pop[[21]]

adm.idx <- names(adm.prob)[21]

logit.prob.frame <- logit(tmp.adm.prob)

setwd(paste0(res_dir,country,'/UR/Classification/Calibration',sep=''))
tmp.shift <- readRDS('BART_MCMC_mse_adm_shift.rds')

x=0
frac.samp.ori <- sapply(logit.prob.frame, function(logit.p.vec){sum(tmp.adm.pop*expit(logit.p.vec+x),na.rm=T)/sum(tmp.adm.pop,na.rm=T)})
frac.ori.dat <- data.frame(Frac=frac.samp.ori,method='tau = 0 (uncalibrated)')

x=-0.6
frac.samp.06 <- sapply(logit.prob.frame, function(logit.p.vec){sum(tmp.adm.pop*expit(logit.p.vec+x),na.rm=T)/sum(tmp.adm.pop,na.rm=T)})
frac.06.dat <- data.frame(Frac=frac.samp.06,method='tau = -0.6')


x=-1.2
frac.samp.12 <- sapply(logit.prob.frame, function(logit.p.vec){sum(tmp.adm.pop*expit(logit.p.vec+x),na.rm=T)/sum(tmp.adm.pop,na.rm=T)})
frac.12.dat <- data.frame(Frac=frac.samp.12,method='tau = -1.2')


frac.dat <- rbind(frac.ori.dat,frac.06.dat,frac.12.dat)
frac.dat$method <- factor(frac.dat$method,levels=c('tau = -1.2','tau = -0.6','tau = 0 (uncalibrated)'))

ggplot(frac.dat, aes(x=Frac, fill=method,color=method)) +
  geom_density( alpha=0.4)


Zomba.density.plot <- 
  ggplot(frac.dat,aes(x=Frac,color=method,fill=method)) +
  geom_density( alpha=0.4)+#scale_color_viridis(discrete = TRUE)+scale_fill_viridis(discrete = TRUE)+
  xlab('urban fraction')+#geom_vline(xintercept=0.153054967, linetype="dashed", size=1,color="darkblue")
  ylab('density')+
  geom_vline(xintercept=0.132215889, linetype="dashed", size=1,color="darkblue")+
  #ggtitle("Distributions of urban fractions based on various values of calibration parameters") +
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=15),
         axis.text.x = element_text(size = 15),
         axis.title.x = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         axis.text.y = element_text(size = 15),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=15),
         strip.text.x = element_text(size = 15))+
  annotate(x=0.132215889,y=+Inf,size = 15/.pt,
           label="Target fraction",vjust=1.3,geom="label")
  
setwd(paste0(res_dir,country,'/Visualization/',sep=''))

pdf(paste0('zomba_calibration_example.pdf'), height = 6, 
    width = 9)
{
  Zomba.density.plot
}
dev.off()

  