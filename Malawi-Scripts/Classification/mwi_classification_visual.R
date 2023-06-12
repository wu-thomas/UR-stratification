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
library(ggmap)
library(data.table)
library(gtable)


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

frame.year<-2008
country.abbrev<-'mwi' #for population rasters
gadm.abbrev = "MWI"

data_dir<-paste0(home_dir,'/Data/')
res_dir<-paste0(home_dir,'/Results/')

# create directories to store results
setwd(paste0(res_dir,country,sep=''))

if(!dir.exists(paths = paste0('Visualization'))){
  dir.create(path = paste0('Visualization'))
}

if(!dir.exists(paths = paste0('Visualization/Maps'))){
  dir.create(path = paste0('Visualization/Maps'))
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

load(file = paste0(poly.path,'/', country, '_Amat.rda'))
load( file = paste0(poly.path, '/', country, '_Amat_Names.rda'))



################################################################
#########   load worldpop
################################################################

worldpop <- raster('worldpop/mwi_ppp_2008_1km_Aggregated_UNadj.tif')


################################################################
#########   load cluster information
################################################################


setwd(paste0(data_dir,country,sep=''))
load('prepared_dat/uncrc_dat_nolikoma.rda')
load('prepared_dat/crc_dat_nolikoma.rda')


################################################################
######### load national grid 
################################################################

setwd(paste0(data_dir,country,sep=''))
load("prepared_dat/natl_grid.rda")

urb_dat$index <- c(1:nrow(urb_dat))

################################################################
#########  load potential locations
################################################################

setwd(paste0(data_dir,country,'/prepared_dat/',sep=''))
load(file='cluster_jitter_frame_nolikoma.rda')
#urban_crc <- readRDS(file='crc_dat.rds')
#urban_crc$light <- urban_crc$noaa


clus_loc_frame <- do.call(rbind, clust_nolikoma)

#rm(clus_loc_list)

# load urban fractions
setwd(paste0(data_dir,country,sep=''))

load('admin2_urb_ref_tab.rda')

################################################################
#########  load predicted probabilities surfaces
################################################################

loss.func='mse'

method='BART_mcmc'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
bart_mcmc_surf<-raster(paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''))
bart_mcmc_uncali_surf<-raster(paste0(method,'_uncalibrated_prob_surf.tif',sep=''))

method='BART_crc'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
bart_crc_surf<-raster(paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''))
bart_crc_uncali_surf <-raster(paste0(method,'_uncalibrated_prob_surf.tif',sep=''))

method='BART_uncrc'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
bart_uncrc_surf<-raster(paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''))
bart_uncrc_uncali_surf<-raster(paste0(method,'_uncalibrated_prob_surf.tif',sep=''))

method='logistic_uncrc'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
logistic_uncrc_surf<-raster(paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''))
logistic_uncrc_uncali_surf<-raster(paste0(method,'_uncalibrated_prob_surf.tif',sep=''))

method='logistic_crc'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
logistic_crc_surf<-raster(paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''))
logistic_crc_uncali_surf<-raster(paste0(method,'_uncalibrated_prob_surf.tif',sep=''))

method='GBT_crc'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
gbt_crc_surf<-raster(paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''))
gbt_crc_uncali_surf<-raster(paste0(method,'_uncalibrated_prob_surf.tif',sep=''))


method='GBT_uncrc'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
gbt_uncrc_surf<-raster(paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''))
gbt_uncrc_uncali_surf<-raster(paste0(method,'_uncalibrated_prob_surf.tif',sep=''))


method='pop'
setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
pop_surf<-raster(paste0(method,'_admin2_thresh_ind_surf.tif',sep=''))


################################################################
#########  helper functions
################################################################

ggmap_extractPoly <- function(ggmap.obj){
  
  # extract box 
  box_range <- attr(ggmap.obj, "bb")
  ymin <- as.numeric(box_range[1])
  ymax <- as.numeric(box_range[3])
  xmin <- as.numeric(box_range[2])
  xmax <- as.numeric(box_range[4])
  
  
  # create rectangle polygon
  
  poly1 <- sp::Polygon(cbind(c(xmin,xmax,xmax,xmin),c(ymin,ymin,ymax, ymax)))
  
  firstPoly <- sp::Polygons(list(poly1), ID = "A")
  
  firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
  
  crs(firstSpatialPoly) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
  
  return(firstSpatialPoly)
  
}


urban_local_map <- function(ras,ggmap.poly,
                              plot.ind =TRUE,
                            p.title=NULL){
  
  
  cropped.ras<-mask(crop(ras,ggmap.poly),
                    ggmap.poly)
  
  ras.dat <- as.data.frame(cropped.ras, xy=TRUE)
  
  if(plot.ind){
    colnames(ras.dat) <- c('x','y','urban_ind')
    ras.dat$urban <- NA
    ras.dat[!is.na(ras.dat$urban_ind)&ras.dat$urban_ind==1,]$urban <- 'urban'
    ras.dat[!is.na(ras.dat$urban_ind)&ras.dat$urban_ind==0,]$urban <- 'rural'
    
  }else {
    colnames(ras.dat) <- c('x','y','urban_prob')
  }
  
  if(plot.ind){
    ras.plot <- ggplot(ras.dat, aes(x = x, y = y, fill = urban)) +
      geom_raster()+
      #scale_fill_manual(values = rev(terrain.colors(2)))+
      scale_fill_manual(values = scales::viridis_pal()(2))+
      coord_fixed()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank())+
      theme(legend.position = "none")+ 
      xlab('lon')+ylab('lat')+
      theme(axis.text.x=element_blank(), #remove x axis labels
            axis.ticks.x=element_blank(), #remove x axis ticks
            axis.text.y=element_blank(),  #remove y axis labels
            axis.ticks.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
      )+
      scale_x_continuous(limits = c(extent(cropped.ras)[1]-1e-10,extent(cropped.ras)[2]+1e-10), expand = c(0, 0)) +
      scale_y_continuous(limits = c(extent(cropped.ras)[3]-1e-10,extent(cropped.ras)[4]+1e-10), expand = c(0, 0)) 
  }else {
    ras.plot <- ggplot(ras.dat, aes(x = x, y = y, fill = urban_prob)) +
      geom_raster()+
      scale_fill_viridis_c(breaks = seq(0,1,by=0.1),direction)+
      coord_fixed()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank())+
      #xlab('lon')+ylab('lat')+
      scale_x_continuous(limits = c(extent(cropped.ras)[1]-1e-10,extent(cropped.ras)[2]+1e-10), expand = c(0, 0)) +
      scale_y_continuous(limits = c(extent(cropped.ras)[3]-1e-10,extent(cropped.ras)[4]+1e-10), expand = c(0, 0)) +
      theme(legend.position = "none")+
      theme(axis.text.x=element_blank(), #remove x axis labels
            axis.ticks.x=element_blank(), #remove x axis ticks
            axis.text.y=element_blank(),  #remove y axis labels
            axis.ticks.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
      )
      
    
  }
  
  if(!is.null(p.title)){
    ras.plot = ras.plot+ggtitle(p.title)+ 
      theme(plot.title = element_text(hjust = 0.5,size=18))
  }
  
  return(ras.plot)
}

prob_to_ind <- function(ras,cutoff=0.5){
  ras.ind <- ras
  values(ras.ind) <- as.numeric(values(ras)>cutoff)
  return(ras.ind)
}


prior_prob_ras <- function(pop.ras,edge.poly,cluster.frame){
  
  ras.tmp <- mask(crop(pop.ras,edge.poly),
                  edge.poly)
  ras.tmp.xy <- coordinates(ras.tmp)
  

  cluster.frame <- as.data.table(cluster.frame)
  
  dist1 <- function(a, b,dt2){
    dt <- data.table((dt2$x-a)^2+(dt2$y-b)^2)
    return(which.min(dt$V1))}
  
  
  results1 <- function() return(cluster.frame[, j = list(Closest =  dist1(x, y,as.data.table(ras.tmp.xy))), 
                                              by = 1:nrow(cluster.frame)])
  
  tmp.pixelID <- results1()
  
  
  cluster.frame$pixelID <- tmp.pixelID$Closest
  
  prior.p.frame <- cluster.frame %>% 
    group_by(pixelID) %>% 
    summarise(prob = sum(norm_w))
  
  
  prior.p.ras <- ras.tmp
  
  values(prior.p.ras)[-prior.p.frame$pixelID] <- NA
  values(prior.p.ras)[prior.p.frame$pixelID] <- prior.p.frame$prob
  
  return(prior.p.ras)
  
  
}
################################################################
#########  Blantyre Prediction
################################################################

blantyre_boundary <- poly.adm2[2,]

### Blantyre overall map

blantyre_bbox <- c(34.7,-16.1, 35.2,-15.3)

blantyre_bmap <-   get_stamenmap(bbox=blantyre_bbox, maptype = "terrain")

setwd(paste0(res_dir,country,'/Visualization/Maps',sep=''))

if(FALSE){
pdf('Blantyre_zoomed_out.pdf', height =6, 
    width = 4)

{

ggmap(blantyre_bmap)+
  geom_polygon(data= fortify(blantyre_boundary),aes(x=long,y=lat),alpha=0.3,colour="black",fill=NA)+
  ggtitle('Blantyre Admin-2 Region Terrain')+ 
    theme(plot.title = element_text(hjust = 0.5))
  
}

dev.off()

}

### Blantyre zoomed-in map

#blantyre_sbox <- c(34.9,-15.9, 35.15,-15.62)
blantyre_sbox <- c(34.95,-15.85, 35.15,-15.65)

blantyre_smap <-   get_map(blantyre_sbox,
                            source="stamen", maptype="terrain", crop=FALSE)

blantyre_terrain_zoomed <- ggmap(blantyre_smap,darken = c(0.05,"black"))+
  ggtitle('Blantyre Terrain')+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()
  )

### prettier blantyre map
if(T){
  blantyre_pmap <- get_googlemap(center = c(35.05,
                                            -15.75), maptype = 'terrain',
                                 zoom=11)
  
  x_coord <- c(34.9,  34.9,  35.15, 35.15, 34.9)
  y_coord <- c(-15.9, -15.62, -15.62,-15.9, -15.9)
  xym <- cbind(x_coord, y_coord)
  xym
  
  p = Polygon(xym)
  ps = Polygons(list(p),1)
  blantyre_poly = SpatialPolygons(list(ps))
  #plot(sps)
  
  blantyre_terrain_zoomed <- ggmap(blantyre_pmap)+
    #geom_polygon(data= fortify(Chikwawa_boundary),aes(x=long,y=lat),alpha=0.3,colour="black",fill=NA)+
    ggtitle('Blantyre Terrain')+ 
    theme(plot.title = element_text(hjust = 0.5,size=18))+
    xlim(c(34.9,35.15))+ylim(c(-15.9,-15.62))+
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()
    )
}

### plot area for rasters
#blantyre_poly <- ggmap_extractPoly (blantyre_smap)
#blantyre_poly <- ggmap_extractPoly (blantyre_pmap)

#plot(blantyre_poly)

### prepare legend 

cropped.mcmc<-mask(crop(bart_mcmc_surf,blantyre_poly),
                   blantyre_poly)

mcmc.dat <- as.data.frame(bart_mcmc_surf, xy=TRUE)
colnames(mcmc.dat) <- c('x','y','urban_prob')

prob.scale.plot <- ggplot(mcmc.dat, aes(x = x, y = y, fill = urban_prob)) +
  geom_raster()+
  scale_fill_viridis_c(breaks = seq(0,1,by=0.1))+
  guides(fill = guide_colourbar(
    title.position = "top"),title.hjust = .5)+
  coord_fixed()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title.align=0.5)+
  labs(fill='urban predicted probability')+
         theme(legend.position="bottom",
         legend.key.width = unit(2.5, 'cm'),
         legend.text=element_text(size=15),
         legend.title=element_text(size=18))
         #legend.text=element_text(angle = 45, vjust = 0, hjust=0.2))

prob.legend <- cowplot::get_legend(prob.scale.plot)
plot(prob.legend)

legend.colour <- gtable_filter(ggplot_gtable(ggplot_build(prob.scale.plot)), "guide-box") 
legend.colour.grob <- grobTree(legend.colour)

#legend.grob(prob.legend)


#scales::viridis_pal()(3)
#barplot(c(2,5), col=c("#440154FF", "#FDE725FF"))



### prepare clipped rasters

blantyre_bart_mcmc_prob <- urban_local_map(ras=bart_mcmc_surf,ggmap.poly=blantyre_poly,
                                     plot.ind =F,p.title = 'BART-MCMC')
blantyre_bart_mcmc_uncali_prob <- urban_local_map(ras=bart_mcmc_uncali_surf,ggmap.poly=blantyre_poly,
                                           plot.ind =F,p.title = 'BART-MCMC')




blantyre_bart_uncrc_prob <-  urban_local_map(ras=bart_uncrc_surf,ggmap.poly=blantyre_poly,
                                             plot.ind =F,p.title = 'BART-uncorrected')
blantyre_bart_uncrc_uncali_prob <- urban_local_map(ras=bart_uncrc_uncali_surf,ggmap.poly=blantyre_poly,
                                                  plot.ind =F,p.title = 'BART-uncorrected')

blantyre_bart_crc_prob <-  urban_local_map(ras=bart_crc_surf,ggmap.poly=blantyre_poly,
                                             plot.ind =F,p.title = 'BART-corrected')
blantyre_bart_crc_uncali_prob <- urban_local_map(ras=bart_crc_uncali_surf,ggmap.poly=blantyre_poly,
                                                   plot.ind =F,p.title = 'BART-corrected')

blantyre_logistic_crc_prob <-  urban_local_map(ras=logistic_crc_surf,ggmap.poly=blantyre_poly,
                                           plot.ind =F,p.title = 'LR-corrected')
blantyre_logistic_crc_uncali_prob <- urban_local_map(ras=logistic_crc_uncali_surf,ggmap.poly=blantyre_poly,
                                                 plot.ind =F,p.title = 'LR-corrected')

blantyre_logistic_uncrc_prob <-  urban_local_map(ras=logistic_uncrc_surf,ggmap.poly=blantyre_poly,
                                               plot.ind =F,p.title = 'LR-uncorrected')
blantyre_logistic_uncrc_uncali_prob <-  urban_local_map(ras=logistic_uncrc_uncali_surf,ggmap.poly=blantyre_poly,
                                                 plot.ind =F,p.title = 'LR-uncorrected')


blantyre_gbt_crc_prob <-  urban_local_map(ras=gbt_crc_surf,ggmap.poly=blantyre_poly,
                                           plot.ind =F,p.title = 'GBT-corrected')
blantyre_gbt_crc_uncali_prob <-  urban_local_map(ras=gbt_crc_uncali_surf,ggmap.poly=blantyre_poly,
                                          plot.ind =F,p.title = 'GBT-corrected')


blantyre_gbt_uncrc_prob <-  urban_local_map(ras=gbt_uncrc_surf,ggmap.poly=blantyre_poly,
                                          plot.ind =F,p.title = 'GBT-uncorrected')

blantyre_gbt_uncrc_uncali_prob <-  urban_local_map(ras=gbt_uncrc_uncali_surf,ggmap.poly=blantyre_poly,
                                            plot.ind =F,p.title = 'GBT-uncorrected')

blantyre_pop_ind <- urban_local_map(ras=pop_surf,ggmap.poly=blantyre_poly,
                                    plot.ind =TRUE,p.title = 'Population Thresholding')


### calibrated predicted probability

#grid.arrange(blantyre_terrain_zoomed, blantyre_mcmc_ind, blantyre_uncrc_ind,blantyre_crc_ind,
#             blantyre_gbt_ind, blantyre_logistic_ind , nrow = 2)

if(FALSE){
setwd(paste0(res_dir,country,'/Visualization/Maps',sep=''))

pdf('Blantyre_compare_calibrated_prob.pdf', height = 12, 
    width = 11)

{
grid.arrange(arrangeGrob(blantyre_bart_uncrc_prob, blantyre_bart_mcmc_prob, blantyre_bart_crc_prob, ncol=3),
             arrangeGrob(blantyre_gbt_uncrc_prob, blantyre_terrain_zoomed, blantyre_gbt_crc_prob, ncol=3),
             arrangeGrob(blantyre_logistic_uncrc_prob, blantyre_pop_ind, blantyre_logistic_crc_prob, ncol=3),
             arrangeGrob(legend.colour.grob),
             nrow=4,
             heights=c(1, 1, 1,1/4))
          

}

dev.off()

setwd(paste0(res_dir,country,'/Visualization/Maps',sep=''))

pdf('Blantyre_compare_uncalibrated_prob.pdf', height = 12, 
    width = 11)

{
  grid.arrange(arrangeGrob(blantyre_bart_uncrc_uncali_prob, blantyre_bart_mcmc_uncali_prob, blantyre_bart_crc_uncali_prob, ncol=3),
               arrangeGrob(blantyre_gbt_uncrc_uncali_prob, blantyre_terrain_zoomed, blantyre_gbt_crc_uncali_prob, ncol=3),
               arrangeGrob(blantyre_logistic_uncrc_uncali_prob, blantyre_pop_ind, blantyre_logistic_crc_uncali_prob, ncol=3),
               arrangeGrob(legend.colour.grob),
               nrow=4,
               heights=c(1, 1, 1,1/4))
  
  
}

dev.off()

### predicted probability
}
