################################################################
#########   load libraries
################################################################

library(gdalUtils)
library(raster)
library(sf)
################################################################
#########   set parameters
################################################################



# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")

# path to AfricaAdmin2Estimates folder 
country<-'Malawi'
main_dir<-'E:/Dropbox/AfricaAdmin2Estimates/Data/'
gadm.abbrev = "MWI"

# survey parameters
beg.year <- 2007   ## 8 year prior to the survey 
end.year <- 2015  ## year of the survey

# population related parameter
country.abbrev<-'mwi'
frame.year<- 2008

# directories
data_dir<-paste0(home_dir,'/Data/')
res_dir<-paste0(home_dir,'/Results/')


setwd(paste(data_dir,country,'/ndvi/tiles/',sep=''))


################################################################
#########   convert hdf to tif 
################################################################

files <- dir(pattern = ".hdf")

filename <- substr(files,1, nchar(files)[1]-4)
filename <- paste0( filename, ".tif")


for (i in 1:length(filename)){
  print(i)
  sds <- get_subdatasets(files[i])
  gdal_translate(sds[1], dst_dataset = filename[i])
}

################################################################
#########  aggregate tiles (across time and space)
################################################################

ndvi_tif_list <- dir(pattern = ".tif")


setwd(paste0(data_dir,country,sep=''))

# UNadjusted population counts
worldpop <- raster(paste0('worldpop/',country.abbrev,
                          '_ppp_',frame.year,
                          '_1km_Aggregated_UNadj.tif',sep=''))



# group tifs from same time
setwd(paste(data_dir,country,'/ndvi/tiles/',sep=''))
ndvi_code<-(substr(ndvi_tif_list,10, 16))
ndvi_group<-unique(ndvi_code)

##### aggregate across space for each time point
for (time_idx in 1:length(ndvi_group)){
print(time_idx)
  
t_ndvi_idx<-which(ndvi_code==ndvi_group[time_idx])
t_ndvi_ras<-ndvi_tif_list[t_ndvi_idx]

# create an empty raster
ext_grid<-extent(worldpop)
template <- raster(ext_grid)
projection(template) <- crs(worldpop)
writeRaster(template, file=paste0(country,'_',ndvi_group[time_idx],"_ndvi.tif"), format="GTiff",
            overwrite=TRUE)

# aggregate ndvi files from the same time point
mosaic_rasters(gdalfile=t_ndvi_ras,dst_dataset=paste0(country,'_',ndvi_group[time_idx],"_ndvi.tif"),of="GTiff")
gdalinfo(paste0(country,'_',ndvi_group[time_idx],"_ndvi.tif"))

temp_ndvi<-raster(paste0(country,'_',ndvi_group[time_idx],"_ndvi.tif"))
proj_ndvi <- projectRaster(temp_ndvi, crs=crs(worldpop))

writeRaster(proj_ndvi, file=paste0(country,'_',ndvi_group[time_idx],"_ndvi.tif"), format="GTiff",
            overwrite=TRUE)

}

#### average across time
ndvi_ras_natl_list<-list()

for (time_idx in 1:length(ndvi_group)){
  
temp_ndvi<-raster(paste0(country,'_',ndvi_group[time_idx],"_ndvi.tif"))
ndvi_ras_natl_list[[time_idx]]<-temp_ndvi
}


stack_ndvi<-stack(ndvi_ras_natl_list)
ndvi_annual_avg <- calc(stack_ndvi, fun = mean, na.rm = T)

setwd(paste(data_dir,country,'/ndvi/',sep=''))

writeRaster(ndvi_annual_avg, paste0(country,'_ndvi','.tif',sep=''))
writeRaster(ndvi_annual_avg, paste0('tiles/',country,'_annual_avg_ndvi','.tif',sep=''))
