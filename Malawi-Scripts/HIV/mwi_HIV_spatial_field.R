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
library(RColorBrewer)

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
######### load fitted models
################################################################


setwd(paste0(res_dir,country,'/HIV',sep=''))
#load(file = 'HIV_female_15_49_updated.RData')

load(file = 'HIV_female_15_49_updated.RData')

################################################################
######### extract spatial fields for iid model
################################################################

nSamp <- 1000
admin.ref <- admin2.names
res.inla <- area.iid.strat.int.admin2


# Draw posterior samples
post.sample = inla.posterior.sample(n = nSamp, result = res.inla)

# Initialize storage
n_admin = dim(admin.ref)[1]
iid.U =  matrix(0, nrow = nSamp, ncol = n_admin)
iid.R =  matrix(0, nrow = nSamp, ncol = n_admin)

# Get indicies in the sample
U.iid.Idx  = res.inla$misc$configs$contents$start[2]
R.iid.Idx  = res.inla$misc$configs$contents$start[3]


for(i in 1:nSamp){
  for(j in 1:n_admin){
    
    iid.U[i,j] <-  post.sample[[i]]$latent[U.iid.Idx+j-1] 
    iid.R[i,j] <-  post.sample[[i]]$latent[R.iid.Idx+j-1] 
    
}
}


IID_sp_obj <- poly.adm2
IID_sp_obj$IID_U_re <- apply(iid.U,2,median)
IID_sp_obj$IID_R_re <- apply(iid.R,2,median)

setwd(paste0(res_dir,'Malawi/HIV/Visualization'))

pdf('IID_UR_RE.pdf', height = 10, 
    width = 10.5)
{
grid.arrange(
  spplot(IID_sp_obj, c("IID_U_re"), at= seq(-1.35, 1.35, by = 0.15),
         col.regions = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50),
         main=list(label="IID RE, Urban",cex=1)),
  
  spplot(IID_sp_obj, c("IID_R_re"), at= seq(-1.35, 1.35, by = 0.15),
         col.regions = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50),
         main=list(label="IID RE, Rural",cex=1)),
  ncol=2)
}
dev.off()


### alternative plotting

IID_sp_map_toplot <- data.frame(re = c(IID_sp_obj$IID_U_re, IID_sp_obj$IID_R_re),
                                 method=c(rep('IID RE, Urban',times=length(IID_sp_obj$IID_R_re)),
                                          rep('IID RE, Rural',times=length(IID_sp_obj$IID_R_re))),
                                 admin.name=c(IID_sp_obj$NAME_2,IID_sp_obj$NAME_2))


IID_sp_map_toplot$method <- factor(IID_sp_map_toplot$method,levels=c('IID RE, Urban',
                                                                       'IID RE, Rural'))


ymax <- 1.35
ymin <- (-1.35)

setwd(paste0(res_dir,'Malawi/HIV/Visualization'))

pdf('IID_UR_RE.pdf',width=8, height = 10)
{
  mapPlot(data = IID_sp_map_toplot, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
          values = c("re"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
          is.long = TRUE, ncol = 2)+
    theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
           legend.text=element_text(size=14),
           legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
           strip.text.x = element_text(size = 15))+
    guides(fill = guide_colourbar(title.position = "top",
                                  title.hjust = .5,
                                  title='Area random effects',
                                  label.position = "bottom"))+
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50) )
}

dev.off()



IID_sp_toplot <- data.frame(rural_re = IID_sp_obj$IID_R_re,urban_re= IID_sp_obj$IID_U_re)

range.toplot <- c(min(c(IID_sp_toplot$rural_re,IID_sp_toplot$urban_re)),
                  max(c(IID_sp_toplot$rural_re,IID_sp_toplot$urban_re)))


setwd(paste0(res_dir,'Malawi/HIV/Visualization'))

pdf('IID_UR_RE_scatter.pdf', height = 8, 
    width = 8)
{
ggplot(data = IID_sp_toplot)+geom_point(aes(x = rural_re, y = urban_re),
                                         color=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.7))+
  coord_fixed() +
  xlim(range.toplot)+ylim(range.toplot)+
  geom_abline(slope=1, intercept=0)+
  xlab('Rural IID random effects')+ylab('Urban IID random effects')+
  theme_bw()+
  theme ( strip.text.x = element_text(size = 15),
          axis.title  = element_text(size = 15),
          axis.text  = element_text(size = 12) )
  
}
dev.off()

################################################################
######### extract spatial fields for BYM2 model
################################################################
nSamp <- 1000
admin.ref <- admin2.names
res.inla <- BYM2.strat.int.admin2
  
  
post.sample = inla.posterior.sample(n = nSamp, result = res.inla)


# Initialize storage
n_admin = dim(admin.ref)[1]
sp.U =  matrix(0, nrow = nSamp, ncol = n_admin)
sp.R =  matrix(0, nrow = nSamp, ncol = n_admin)




spat.u.Idx = res.inla$misc$configs$contents$start[2]
spat.r.Idx = res.inla$misc$configs$contents$start[3]
intIdx  = res.inla$misc$configs$contents$start[5]
urbIdx  = res.inla$misc$configs$contents$start[6]



for(i in 1:nSamp){
  
  for(j in 1:n_admin){
    # spatial effects (spatial+iid)
    sp.U[i,j] = post.sample[[i]]$latent[spat.u.Idx + j-1]
    sp.R[i,j] = post.sample[[i]]$latent[spat.r.Idx + j-1]
   }
}

BYM2_sp_obj <- poly.adm2
BYM2_sp_obj$bym2_U_re <- apply(sp.U,2,median)
BYM2_sp_obj$bym2_R_re <- apply(sp.R,2,median)


setwd(paste0(res_dir,'Malawi/HIV/Visualization'))

pdf('BYM2_UR_RE.pdf', height = 10, 
    width = 10.5)
{
grid.arrange(
spplot(BYM2_sp_obj, c("bym2_U_re"), at= seq(-1.25, 1.25, by = 0.25),
       col.regions = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50),
       main=list(label="Spatial+IID RE, BYM2, Urban",cex=1)),

spplot(BYM2_sp_obj, c("bym2_R_re"),at= seq(-1.25, 1.25, by = 0.25),
       col.regions = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50),
       main=list(label="Spatial+IID RE, BYM2, Rural",cex=1)),
 ncol=2)
}

dev.off()



### alternative plotting

BYM2_sp_map_toplot <- data.frame(re = c(BYM2_sp_obj$bym2_R_re, BYM2_sp_obj$bym2_U_re),
                                 method=c(rep('Spatial+IID RE, BYM2, Rural',times=length(BYM2_sp_obj$bym2_R_re)),
                                          rep('Spatial+IID RE, BYM2, Urban',times=length(BYM2_sp_obj$bym2_R_re))),
                                 admin.name=c(BYM2_sp_obj$NAME_2,BYM2_sp_obj$NAME_2))


BYM2_sp_map_toplot$method <- factor(BYM2_sp_map_toplot$method,levels=c('Spatial+IID RE, BYM2, Urban',
                                                                       'Spatial+IID RE, BYM2, Rural'))


ymax <- 1.25
ymin <- (-1.25)

setwd(paste0(res_dir,'Malawi/HIV/Visualization'))

pdf('BYM2_UR_RE.pdf',width=8, height = 10)
{
mapPlot(data = BYM2_sp_map_toplot, geo = poly.adm2, variables = c("method"), size=0.1,ylim= c(ymin,ymax),
                     values = c("re"), by.data = "admin.name", by.geo = "NAME_2",direction = -1, legend.label = "HIV",
                     is.long = TRUE, ncol = 2)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='Area random effects',
                                label.position = "bottom"))+
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50) )
}

dev.off()





BYM2_sp_toplot <- data.frame(rural_re = BYM2_sp_obj$bym2_R_re,urban_re= BYM2_sp_obj$bym2_U_re)

range.toplot <- c(min(c(BYM2_sp_toplot$rural_re,BYM2_sp_toplot$urban_re)),
                  max(c(BYM2_sp_toplot$rural_re,BYM2_sp_toplot$urban_re)))



setwd(paste0(res_dir,'Malawi/HIV/Visualization'))

pdf('BYM2_UR_RE_scatter.pdf', height = 8, 
    width = 8)
{
ggplot(data = BYM2_sp_toplot)+geom_point(aes(x = rural_re, y = urban_re),
                                           color=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.7))+
  coord_fixed() +
  xlim(range.toplot)+ylim(range.toplot)+
  geom_abline(slope=1, intercept=0)+
  xlab('Rural spatial random effects')+ylab('Urban spatial random effects')+
  theme_bw()+
  theme ( strip.text.x = element_text(size = 15),
          axis.title  = element_text(size = 15),
          axis.text  = element_text(size = 12) )
}

dev.off()