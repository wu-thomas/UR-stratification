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


################################################################
#########   set parameters
################################################################

## set directory

country<-'Malawi'
# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")

# path to AfricaAdmin2Estimates folder 
main_dir<-'E:/Dropbox/AfricaAdmin2Estimates/Data/'
gadm.abbrev = "MWI"

# survey parameters
beg.year <- 2007   ## 8 year prior to the survey 
end.year <- 2015  ## year of the survey

# population related parameter
country.abbrev<-'mwi'
frame.year<- 2008



################################################################
#########  file names (no modification is needed)
################################################################

# directories
data_dir<-paste0(home_dir,'/Data/')
res_dir<-paste0(home_dir,'/Results/')
pop_dir<-paste0(main_dir,'Population/')
cluster_dir<-paste0(main_dir,'countryDataFolders/',country)

# night time light (noaa) file name
noaa.file<-'F162008.v4b_web.stable_lights.avg_vis.tif'

# elevation
elev.folder<-'wc2.1_30s_elev'
elev.file<-paste0(elev.folder,'.tif',sep='')

# average temperature
tavg.folder<-'wc2.1_30s_tavg'
tavg.file<-paste0('global_annual_tavg','.tif',sep='')

# precipitation
prec.folder<-'wc2.1_30s_prec'
prec.file<-paste0('global_annual_prec','.tif',sep='')

# access
access.folder<-'access_50k'
access.file<-paste0('acc_50k','.tif',sep='')

# ndvi
ndvi.file<-paste0('ndvi/tiles/',country,'_annual_avg_ndvi','.tif',sep='')


################################################################
#########   load polygons
################################################################

#### Load polygon files ####
setwd(paste(data_dir,country,sep=''))

poly.path <- paste0("shapeFiles_gadm")
poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_")
poly.layer.adm1 <- paste('gadm36', gadm.abbrev,
                         '1', sep = "_")
poly.layer.adm2 <- paste('gadm36', gadm.abbrev,
                         '2', sep = "_")

poly.adm0 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm0))
poly.adm1 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm1))
poly.adm2 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm2))

proj4string(poly.adm0) <- proj4string(poly.adm1) <- proj4string(poly.adm2)
load(paste0('shapeFiles_gadm/', country, '_Amat.rda'))
load(paste0('shapeFiles_gadm/', country, '_Amat_Names.rda'))


################################################################
#########   load worldpop
################################################################

setwd(paste(pop_dir,country,sep=''))

# UNadjusted population counts
worldpop <- raster(paste0(country.abbrev,
                          '_ppp_',frame.year,
                          '_1km_Aggregated_UNadj.tif',sep=''))


#setwd(paste0(data_dir,country,sep=''))
#writeRaster(worldpop, paste0('prepared_raster/',country,'_frame_total_pop','.tif',sep=''))



################################################################
#########   prepare rasters
################################################################

# create directory
setwd(paste0(data_dir,country,sep=''))

if(!dir.exists(paths = paste0('prepared_raster'))){
  dir.create(path = paste0('prepared_raster'))
}

buffered_adm0 <- buffer(poly.adm0, width=0.25)

# clip global rasters to individual country

plot(poly.adm0)
plot(buffered_adm0, add=TRUE)

################################################################
#########   night time light
################################################################
setwd(paste(data_dir,'/Global_raster/',sep=''))

# load global noaa
noaa_g <- raster(paste0('noaa/',noaa.file))

#light <- projectRaster(noaa_g, crs=crs(worldpop))

light_buffered<-mask(crop(noaa_g,buffered_adm0),
                     buffered_adm0)

setwd(paste0(data_dir,country,sep=''))
writeRaster(light_buffered, paste0('prepared_raster/',country.abbrev,'_light_buffered','.tif',sep=''))

#plot(light_buffered)
  
  

################################################################
#########   average temperature
################################################################

setwd(paste(data_dir,'/Global_raster/',sep=''))

tavg_global <- raster(paste0(tavg.folder,'/',tavg.file,sep=''))
#tavg <- projectRaster(tavg_global, crs=crs(worldpop))

tavg_buffered<-mask(crop(tavg_global,buffered_adm0),
                    buffered_adm0)

#plot(tavg_buffered)
setwd(paste0(data_dir,country,sep=''))
writeRaster(tavg_buffered, paste0('prepared_raster/',country.abbrev,'_tavg_buffered','.tif',sep=''))

#plot(tavg_buffered)

################################################################
#########   elevation
################################################################

setwd(paste(data_dir,'/Global_raster/',sep=''))

elev_global <- raster(paste0(elev.folder,'/',elev.file,sep=''))

elev_buffered<-mask(crop(elev_global,buffered_adm0),
                    buffered_adm0)

#plot(elev_buffered)
setwd(paste0(data_dir,country,sep=''))
writeRaster(elev_buffered, paste0('prepared_raster/',country.abbrev,'_elev_buffered','.tif',sep=''))

#plot(elev_buffered)

################################################################
#########  precipitation
################################################################

setwd(paste(data_dir,'/Global_raster/',sep=''))

prec_global <- raster(paste0(prec.folder,'/',prec.file,sep=''))

prec_buffered<-mask(crop(prec_global,buffered_adm0),
                    buffered_adm0)

#plot(prec_buffered)
setwd(paste0(data_dir,country,sep=''))
writeRaster(prec_buffered, paste0('prepared_raster/',country.abbrev,'_prec_buffered','.tif',sep=''))

#plot(prec)


################################################################
#########  distance to nearest city
################################################################

setwd(paste(data_dir,'/Global_raster/',sep=''))

access_global <- raster(paste0(access.folder,'/',access.file,sep=''))

access_buffered<-mask(crop(access_global,buffered_adm0),
                      buffered_adm0)

#plot(access_buffered)
setwd(paste0(data_dir,country,sep=''))
writeRaster(access_buffered, paste0('prepared_raster/',country.abbrev,'_access_buffered','.tif',sep=''))

#plot(access)


################################################################
#########   ndvi
################################################################

setwd(paste0(data_dir,country,sep=''))

ndvi<- raster(paste0(ndvi.file))

ndvi_buffered<-mask(crop(ndvi,buffered_adm0),
           buffered_adm0)

setwd(paste0(data_dir,country,sep=''))
writeRaster(ndvi, paste0('prepared_raster/',country.abbrev,'_ndvi_buffered','.tif',sep=''))

#plot(ndvi_buffered)





### The following is for aggreagation model, could be ran latter

################################################################
#########   process U5 pop
################################################################

years <- c(beg.year:end.year)

for ( t in 1:9){

print(t)
year <- years[t]

setwd(paste(pop_dir,country,sep=''))
pop_surf<-raster(paste0(country.abbrev,'_u5_',year,'_100m.tif'))


pop_u5_aggregate <- aggregate(pop_surf, fact=10,sum)

pop_grid<-as.data.frame(coordinates(worldpop))
colnames(pop_grid)<-c('x','y')

pop_grid$u5_pop <- raster::extract(pop_u5_aggregate, 
                                   pop_grid[c('x','y')])

u5_pop<-worldpop
values(u5_pop)<-pop_grid$u5_pop 

writeRaster(u5_pop, overwrite=TRUE,
            paste0(country.abbrev,'_u5_',year,'_1km.tif'))

}