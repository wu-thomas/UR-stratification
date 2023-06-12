##############################################################################
######### calibrate probabilities based on all samples, national
##############################################################################

natl.cali.samp <- function(prob.frame,pop.vec,bench.frac=NULL,loss.func='mse'){
  
  
  logit.frame <- logit(prob.frame)
  
  if(loss.func=='mse'){
    bench.loss <- function(x, pop=NULL,logit.prob.frame=NULL,bench.frac=NULL)
    {
      
      frac.samp <- sapply(logit.prob.frame, function(logit.p.vec){sum(pop*expit(logit.p.vec+x),na.rm=T)/sum(pop,na.rm=T)
      } )
      
      loss <- mean((frac.samp-bench.frac)^2)
      loss
    }
  }else{
    
    if(loss.func=='mean'){
      
      
      bench.loss <- function(x, pop=NULL,logit.prob.frame=NULL,bench.frac=NULL)
      {
        
        frac.samp <- sapply(logit.prob.frame, function(logit.p.vec){sum(pop*expit(logit.p.vec+x),na.rm=T)/sum(pop,na.rm=T)
        } )
        
        loss <- mean(frac.samp-bench.frac)
        loss
      }
    }else{
      
      if(loss.func=='abs'){
        
        bench.loss <- function(x, pop=NULL,logit.prob.frame=NULL,bench.frac=NULL)
        {
          
          frac.samp <- sapply(logit.prob.frame, function(logit.p.vec){sum(pop*expit(logit.p.vec+x),na.rm=T)/sum(pop,na.rm=T)
          } )
          
          loss <- mean(abs(frac.samp-bench.frac))
          loss
        }
        
      }else{
        
        if(loss.func=='median'){
          
          
          bench.loss <- function(x, pop=NULL,logit.prob.frame=NULL,bench.frac=NULL)
          {
            
            frac.samp <- sapply(logit.prob.frame, function(logit.p.vec){sum(pop*expit(logit.p.vec+x),na.rm=T)/sum(pop,na.rm=T)
            } )
            
            loss <- median(frac.samp-bench.frac)
            loss
          }
        
        
      }
      }
    }
    
  }
  
  
  xstart <- c(0)
  
  bench.sol <- nleqslv(xstart, bench.loss, control=list(btol=.0001),method="Newton",
                       pop=pop.vec,logit.prob.frame=logit.frame,bench.frac=bench.frac)
  
  
  calibrated.prob <- expit(logit.frame+bench.sol$x)
  
  
  return.obj <- list()
  return.obj[[1]] <- bench.sol$x
  return.obj[[2]] <- calibrated.prob
  names(return.obj) <- c('tau.tab','calibrated.prob')
  
  return(return.obj)
  
}




##############################################################################
######### calibrate probabilities based on all samples, subnational
##############################################################################

subnatl.cali.samp <- function(prob.frame,pop.vec,adm.vec,
                              loss.func='mse',ref.tab){
  
  adm.prob <- split( prob.frame , f = adm.vec )
  adm.pop <- split( pop.vec , f = adm.vec )
  
  
  adm.find.tau <- function(idx){
    tmp.adm.prob <- adm.prob[[idx]]
    tmp.adm.pop <- adm.pop[[idx]]
    
    adm.idx <- names(adm.prob)[idx]
    bench.frac <- ref.tab[ref.tab$Internal==adm.idx,]$urb_frac
    
    #return(names(adm.prob)[idx])
    p.shift <- natl.cali.samp(prob.frame=tmp.adm.prob,
                              pop.vec=tmp.adm.pop,
                              bench.frac =bench.frac,
                              loss.func=loss.func)
    
    return(p.shift[[1]])
  }
  
  urb.list<-lapply(seq_along(adm.prob), FUN=adm.find.tau)
  
  urb.class <- do.call("rbind", urb.list)
  
  adm.calibrate.prob <- function(idx){
    tmp.adm.prob <- adm.prob[[idx]]
    tmp.adm.tau <- urb.class[idx]
    
    calibrate.prob <- expit(logit(tmp.adm.prob)+tmp.adm.tau)
    
    return(calibrate.prob)
    
  }
  
  calibrated.prob.list<-lapply(seq_along(adm.prob), FUN=adm.calibrate.prob)
  
  calibrated.prob <- do.call("rbind", calibrated.prob.list)
  calibrated.prob <- calibrated.prob[ order(as.numeric(row.names(calibrated.prob))), ]
  
  
  return.obj <- list()
  return.obj[[1]] <- data.frame(admin=names(adm.prob),tau=urb.class)
  return.obj[[2]] <- calibrated.prob
  names(return.obj) <- c('tau.tab','calibrated.prob')
  
  return(return.obj)
  
}






##############################################################################
######### create probability and indicator surface raster
##############################################################################



get.urb.surf <- function(natl.grid,pred.vec=NULL,ras.template){
  
  # assign classification based on 0.5 cutoff
  natl.grid$urban <- NA
  natl.grid[natl.grid$complete,]$urban <- (pred.vec >0.5)
  natl.grid[!natl.grid$complete,]$pop_den <- NA    # omit pixels without classification
  
  # prepare rasters for classification surfaces
  ind.surf <- ras.template
  values(ind.surf) <- natl.grid$urban
  
  # prepare return object
  return.obj <- list()
  return.obj[[1]] <- natl.grid
  return.obj[[2]] <- ind.surf
  names(return.obj) <- c('grid_dat','indicator_surface')
  
  # prepare prob surf
  natl.grid$pred_prob <- NA
  natl.grid[natl.grid$complete,]$pred_prob <- pred.vec
  
  pred.p.surf <- ras.template
  values(pred.p.surf) <- natl.grid$pred_prob
  
  return.obj[[3]] <- pred.p.surf
  names(return.obj) <- c('grid_dat','indicator_surface','probability_surface')
  
  
  
  return(return.obj)
}




##############################################################################
######### calculate fractions based on samples 
##############################################################################



urb.frac.sample <- function(pixel.grid,pop.ras,urb.prob.frame,
                            admin.poly,admin.level,admin.name){
  # assign population density
  pixel.grid$pop_den<-raster::extract(pop.ras, pixel.grid[c('x','y')])
  
  pixel.grid <- pixel.grid %>% mutate_at(vars('pop_den'), ~replace(., is.na(.), 0))
  #pixel.grid[is.na(pixel.grid$pop_den),]$pop_den <- 0
  
  # prob.to.frac <- function(prob.vec,pop.vec){
  
  #   return(sum(pop.vec*prob.vec,na.rm=T)/sum(pop.vec,na.rm=T))
  
  # }
  
  # calculate national fraction
  if(is.null(admin.poly)){
    
    natl.frac.samp <- t(urb.prob.frame)%*%complete_urb_dat$pop_den/sum(complete_urb_dat$pop_den,na.rm=T)
    
    # prepare return object
    return.obj <- list()
    return.obj[[1]] <- as.vector(natl.frac.samp)
    return.obj[[2]] <- median(natl.frac.samp)
    names(return.obj) <- c('post.sample','post.median')
    
    
    return(return.obj)  
  }
  
  
  # assign admin regions
  points.frame <- as.data.frame(pixel.grid[c('x','y')])
  points.frame <- SpatialPoints(points.frame)
  proj4string(points.frame) = proj4string(admin.poly)
  
  poly.over.adm <- SpatialPolygons(admin.poly@polygons)
  proj4string(poly.over.adm) = proj4string(admin.poly)
  
  admin = over(points.frame, poly.over.adm)
  
  pixel.grid$adm.num <- admin
  pixel.grid$adm.char <-NA
  pixel.grid[!is.na(pixel.grid$adm.num),]$adm.char<- paste0(admin.level,'_',pixel.grid[!is.na(pixel.grid$adm.num),]$adm.num,sep='')
  
  
  adm.prob <- split( urb.prob.frame , f = pixel.grid$adm.char )
  adm.pop <- split( pixel.grid$pop_den , f = pixel.grid$adm.char )
  
  find.urb.frac <- function(idx){
    tmp.adm.prob <- adm.prob[[idx]]
    tmp.adm.pop <- adm.pop[[idx]]

    adm.frac.samp <- t(tmp.adm.prob)%*%tmp.adm.pop/sum(tmp.adm.pop,na.rm=T)
    
    return(adm.frac.samp)
    
  }
  
  frac.samp.list<-lapply(seq_along(adm.prob), FUN=find.urb.frac)
  
  frac.samp <- do.call("cbind", frac.samp.list)
  rownames(frac.samp) <- NULL
  colnames(frac.samp) <- names(adm.prob)
  
  
  frac.median <- apply(frac.samp,2,median)
  frac.median.frame <- data.frame(admin=names(adm.prob),Frac=as.vector(frac.median))
  
  # prepare return object
  return.obj <- list()
  return.obj[[1]] <- frac.samp
  return.obj[[2]] <- frac.median.frame
  names(return.obj) <- c('post.sample','post.median')
  
  
  return(return.obj)  
  
  
}



##############################################################################
######### create surface for uncalibrated probabilities 
##############################################################################


get.uncali.surface <- function(method,
                               prob.frame,
                               natl.dat,
                               ras.template,
                               pop.ras
){
  
  uncalibrated.prob.samp <- as.data.frame(prob.frame)
    summarize.uncali.prob <- apply(uncalibrated.prob.samp, 1, median, na.rm=T)
    
    surf.and.dat <- get.urb.surf(natl.grid=natl.dat,
                                 pred.vec=summarize.uncali.prob,
                                 ras.template=ras.template)
    
    setwd(paste0(res_dir,country,'/UR/Classification/',sep=''))
    if(!dir.exists(paths = paste0(method))){
      dir.create(path = paste0(method))
    }
    setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
    
    
    writeRaster(surf.and.dat$indicator_surface, file=paste0(method,'_','uncalibrated_ind_surf.tif',sep=''), 
                format="GTiff",overwrite=TRUE)
    writeRaster(surf.and.dat$probability_surface, file=paste0(method,'_','uncalibrated_prob_surf.tif',sep=''),
                format="GTiff",overwrite=TRUE)

}

  
##############################################################################
######### pipeline for national benchmark based on samples
##############################################################################

natl.cali.samp.pipeline <- function(method,
                                    prob.frame,
                                    bench.frac,
                                    loss.func='mse',
                                    natl.dat,
                                    ras.template,
                                    pop.ras,
                                    admin.poly=NULL,
                                    admin.level=NULL,
                                    admin.name=NULL
){
  
  complete.natl.dat <- natl.dat[natl.dat$complete,]
  pixel.grid <- complete.natl.dat[,c('x','y')]
  pop.vec <- complete.natl.dat$pop_den
  prob.frame <- as.data.frame(prob.frame)
  
  natl.cali.res <- natl.cali.samp(prob.frame=prob.frame,
                                  pop.vec=pop.vec,
                                  bench.frac=bench.frac,
                                  loss.func=loss.func)
  
  
  # create directory
  
  
  calibrated.prob.samp <- natl.cali.res[[2]]
  
  setwd(paste0(res_dir,country,'/UR/Classification/Calibration',sep=''))
  
  saveRDS(natl.cali.res[[1]],file=paste0(method,'_',loss.func,
                                         '_natl_shift.rds',sep=''))
  
  calibrated.prob.samp <- as.data.frame(calibrated.prob.samp)
  summarize.cali.prob <- apply(calibrated.prob.samp, 1, median, na.rm=T)
  
  surf.and.dat <- get.urb.surf(natl.grid=natl.dat,
                               pred.vec=summarize.cali.prob,
                               ras.template=ras.template)
  
  setwd(paste0(res_dir,country,'/UR/Classification/',sep=''))
  
  if(!dir.exists(paths = paste0(method))){
    dir.create(path = paste0(method))
  }
  setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
  
  
  writeRaster(surf.and.dat$indicator_surface, file=paste0(method,'_',loss.func,'_natl_bench_ind_surf.tif',sep=''), 
              format="GTiff",overwrite=TRUE)
  writeRaster(surf.and.dat$probability_surface, file=paste0(method,'_',loss.func,'_natl_bench_prob_surf.tif',sep=''),
              format="GTiff",overwrite=TRUE)
  save(calibrated.prob.samp,file=paste0(method,'_',loss.func,'_natl_bench_calibrated_psamp.rda',sep=''))
  
  
  frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                              pop.ras = pop.ras,
                              urb.prob.frame = calibrated.prob.samp,
                              admin.poly = admin.poly ,
                              admin.level = admin.level,
                              admin.name=admin.name)
  
  
  setwd(paste0(res_dir,country,'/UR/Fractions',sep=''))
  
  if(!dir.exists(paths = paste0(method))){
    dir.create(path = paste0(method))
  }
  setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))
  
  
  
  save(frac.adm,file=paste0(method,'_',loss.func,'_natl_bench_total_',admin.level,'_frac.rda'))
  
  
  
}











##############################################################################
######### pipeline for subnational benchmark based on samples
##############################################################################

adm.cali.samp.pipeline <- function(method,
                                   prob.frame,
                                   bench.frac,
                                   loss.func='mse',
                                   natl.dat,
                                   ref.tab,
                                   ras.template,
                                   pop.ras,
                                   admin.poly=NULL,
                                   admin.level=NULL,
                                   admin.name=NULL){
  
  complete.natl.dat <- natl.dat[natl.dat$complete,]
  pixel.grid <- complete.natl.dat[,c('x','y')]
  pop.vec <- complete.natl.dat$pop_den
  adm.vec <- complete.natl.dat$bench_adm
  prob.frame <- as.data.frame(prob.frame)
  rownames(prob.frame) <- c(1:nrow(prob.frame))
    
  adm.cali.res <- subnatl.cali.samp(prob.frame=prob.frame,
                                    pop.vec=pop.vec,
                                    adm.vec=adm.vec,
                                    loss.func=loss.func,
                                    ref.tab=ref.tab)
  
  setwd(paste0(res_dir,country,'/UR/Classification/Calibration',sep=''))
  
  calibrated.prob.samp <- adm.cali.res[[2]]
  calibrated.prob.samp <- as.data.frame(calibrated.prob.samp)
  
  saveRDS(adm.cali.res[[1]],file=paste0(method,'_',loss.func,
                                        '_adm_shift.rds',sep=''))
  
  summarize.cali.prob <- apply(calibrated.prob.samp, 1, median, na.rm=T)
  
  surf.and.dat <- get.urb.surf(natl.grid=natl.dat,
                               pred.vec=summarize.cali.prob,
                               ras.template=ras.template)
  
  setwd(paste0(res_dir,country,'/UR/Classification/',sep=''))
  
  if(!dir.exists(paths = paste0(method))){
    dir.create(path = paste0(method))
  }
  setwd(paste0(res_dir,country,'/UR/Classification/',method,sep=''))
  
  writeRaster(surf.and.dat$indicator_surface, file=paste0(method,'_',loss.func,'_adm_bench_ind_surf.tif',sep=''), 
              format="GTiff",overwrite=TRUE)
  writeRaster(surf.and.dat$probability_surface, file=paste0(method,'_',loss.func,'_adm_bench_prob_surf.tif',sep=''),
              format="GTiff",overwrite=TRUE)
  save(calibrated.prob.samp,file=paste0(method,'_',loss.func,'_adm_bench_calibrated_psamp.rda',sep=''))
  
  
  frac.adm <- urb.frac.sample(pixel.grid = pixel.grid, ## must match the order in the urban probability frame!
                              pop.ras = pop.ras,
                              urb.prob.frame = calibrated.prob.samp,
                              admin.poly = admin.poly ,
                              admin.level = admin.level,
                              admin.name=admin.name)
  
  setwd(paste0(res_dir,country,'/UR/Fractions',sep=''))
  
  if(!dir.exists(paths = paste0(method))){
    dir.create(path = paste0(method))
  }
  setwd(paste0(res_dir,country,'/UR/Fractions/',method,sep=''))
  
  save(frac.adm,file=paste0(method,'_',loss.func,'_adm_bench_total_',admin.level,'_frac.rda'))
  
  
  
}








