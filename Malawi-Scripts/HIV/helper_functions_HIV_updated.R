###############################################################################
# functions.R                                                                 #
#    Help functions for EThiopia vaccination case study.                      #
# adapted from: Geir-Arne Fuglstad <geir-arne.fuglstad@ntnu.no>               #
###############################################################################

# Required libraries
library(survey)
library(INLA)
library(ggplot2)
library(raster)
library(readxl)


#setwd(paste0(code_dir,'/',country,'/HIV'))
#source('helper_functions_HIV_updated.R')

################################################################
#########  Direct estimates admin1
################################################################

# Compute national direct estimates

# Compute admin1 direct estimates
getHIVDirect = function(myData, byUrban = FALSE){
  
  myData$adm1_ur <- as.factor(myData$adm1_ur )
  # Get direct stratum-specific direct estimates
  my.svydesign <- svydesign(id= ~clusterIdx ,
                            strata=~stratum, nest=T, 
                            weights=~weight, data=myData)
  
  if(byUrban){
    direct.res = svyglm(HIV~(-1)+adm1_ur,
                        design = my.svydesign,
                        family = quasibinomial)
    strat.admin1 = rep(sort(unique(myData$admin1.char)), each = 2)
    strat.admin1 = strat.admin1
    urb.admin1 = rep(c("R", "U"), dim(admin1.names)[1])
    urb.admin1 = urb.admin1
    direct.est = data.frame(logitP = direct.res$coefficients,
                            se = sqrt(diag(vcov(direct.res))),
                            admin1 = strat.admin1,
                            urb = urb.admin1)
  }else{
    direct.res = svyglm(HIV~(-1)+admin1Fac,
                        design = my.svydesign,
                        family = quasibinomial)
    strat.admin1 = levels(myData$admin1Fac)
    direct.est = data.frame(logitP = direct.res$coefficients,
                            se = sqrt(diag(vcov(direct.res))), 
                            admin1 = strat.admin1)
  }
  
  # Get confidence intervals
  k = qnorm(0.975)
  direct.est$logitP_Upp = direct.est$logitP + k*direct.est$se
  direct.est$logitP_Low = direct.est$logitP - k*direct.est$se
  
  # Map to probability
  direct.est$p_Low = 1/(1+exp(-direct.est$logitP_Low))
  direct.est$p_Med = 1/(1+exp(-direct.est$logitP))
  direct.est$p_Upp = 1/(1+exp(-direct.est$logitP_Upp))
  
  # Sort admin1
  if(byUrban){
    idx = order(direct.est$admin1, direct.est$urb)
    row.names(direct.est) <- paste0(direct.est$admin1,':',direct.est$urb,sep='')
  }else{
    idx = order(direct.est$admin1)
  }
  direct.est = direct.est[idx,]
  
  return(direct.est)
}

aggreDirect = function(urb.res,urb.frac){
  
  aggre.est <- vector()
  
  for (i in 1:dim(urb.frac)[1]){
    
    adm.name <- urb.frac$admin.name[i]
    frac <- urb.frac$Frac[i]
    
    if(frac == 1){
      urban.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='U',]$p_Med
      aggre.est[i] <- urban.est
      
    }else{
      
      if(frac == 0){
        rural.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='R',]$p_Med
        aggre.est[i] <- rural.est
        
      }else{
        
        urban.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='U',]$p_Med     
        rural.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='R',]$p_Med
        
        aggre.est[i] <- urban.est*frac+ rural.est*(1-frac)
      }
      
    }
    
  }
  
  return (data.frame(admin1 = urb.frac$admin.name,
                     aggre.est = aggre.est, 
                     urb.frac =  urb.frac$Frac))
}


################################################################
#########  Direct estimates admin2
################################################################
getHIVDirectAdm2 = function(myData, byUrban = FALSE,save.draws = FALSE){
  
  myData$adm2_ur <- as.factor(myData$adm2_ur )
  # Get direct stratum-specific direct estimates
  my.svydesign <- svydesign(id= ~clusterIdx ,
                            strata=~stratum, nest=T, 
                            weights=~weight, data=myData)
  
  if(byUrban){
    direct.res = svyglm(HIV~(-1)+adm2_ur,
                        design = my.svydesign,
                        family = quasibinomial)
    strat.admin2 = rep(sort(unique(myData$admin2.char)), each = 2)
    strat.admin2 = strat.admin2
    urb.admin2 = rep(c("R", "U"), dim(admin2.names)[1])
    urb.admin2 = urb.admin2
    direct.est = data.frame(logitP = direct.res$coefficients,
                            se = sqrt(diag(vcov(direct.res))),
                            admin2 = strat.admin2,
                            urb = urb.admin2)
  }else{
    direct.res = svyglm(HIV~(-1)+admin2.char,
                        design = my.svydesign,
                        family = quasibinomial)
    strat.admin2 = sort(unique(myData$admin2.char))
    direct.est = data.frame(logitP = direct.res$coefficients,
                            se = sqrt(diag(vcov(direct.res))), 
                            admin2 = strat.admin2)
  }
  
  # Get confidence intervals
  k = qnorm(0.975)
  direct.est$logitP_Upp = direct.est$logitP + k*direct.est$se
  direct.est$logitP_Low = direct.est$logitP - k*direct.est$se
  
  # Map to probability
  direct.est$p_Low = 1/(1+exp(-direct.est$logitP_Low))
  direct.est$p_Med = 1/(1+exp(-direct.est$logitP))
  direct.est$p_Upp = 1/(1+exp(-direct.est$logitP_Upp))
  
  # Sort admin2
  if(byUrban){
    idx = order(direct.est$admin2, direct.est$urb)
    row.names(direct.est) <- paste0(direct.est$admin2,':',direct.est$urb,sep='')
  }else{
    idx = order(direct.est$admin2)
  }
  direct.est = direct.est[idx,]
  
  if(!save.draws){
  return(direct.est)
  }
  
  draws.out <- matrix(,nrow=1000,ncol=dim(direct.est)[1])
  #colnames(draws.out) <- row.names(direct.est)
  for(i in 1:dim(direct.est)[1]){
    normal.samp <- rnorm(1000)
    adm.samp = expit(direct.est$logitP[i] + normal.samp*direct.est$se[i])
    draws.out[,i]<-adm.samp
  }
  draws.out <- as.data.frame(draws.out)
  colnames(draws.out) <- row.names(direct.est)
  
  return(list(overall=direct.est,draws=draws.out))
}

aggreDirectAdm2  = function(urb.res,urb.frac){
  
  aggre.est <- vector()
  
  for (i in 1:dim(urb.frac)[1]){
    
    adm.name <- urb.frac$admin.name[i]
    frac <- urb.frac$Frac[i]
    
    if(frac == 1){
      urban.est <- urb.res[urb.res$admin2==adm.name & urb.res$urb =='U',]$p_Med
      aggre.est[i] <- urban.est
      
    }else{
      
      if(frac == 0){
        rural.est <- urb.res[urb.res$admin2==adm.name & urb.res$urb =='R',]$p_Med
        aggre.est[i] <- rural.est
        
      }else{
        
        urban.est <- urb.res[urb.res$admin2==adm.name & urb.res$urb =='U',]$p_Med     
        rural.est <- urb.res[urb.res$admin2==adm.name & urb.res$urb =='R',]$p_Med
        
        aggre.est[i] <- urban.est*frac+ rural.est*(1-frac)
      }
      
    }
    
  }
  
  return (data.frame(admin2 = urb.frac$admin.name,
                     aggre.est = aggre.est, 
                     urb.frac =  urb.frac$Frac))
}



#' @ur.draws UR specific draws, dim: n_samp*n_admin, colname names: 'admin_name:U/R' (ex. Blantyre:R) 
#' @urb.frac UR fractions, if fixed, n_admin rows, if samples, n_samp*n_admin
#' @fixed.urb Indicator whether is fixed 
#' @nSamp number of samples 
#' @admin.ref reference table for admin names and admin index
aggreDirectAdm2_CI  = function(ur.draws,urb.frac,fixed.urb=T,nSamp=1000,
                               admin.ref=NULL){
  
  ## fixed urban proportions
  if(fixed.urb){
    
  n_admin <- dim(urb.frac)[1]  
  aggre.samp <- matrix(0, nrow = nSamp, ncol = n_admin)
  
  for (i in 1:n_admin){
    
    adm.name <- urb.frac$admin.name[i]
    frac <- urb.frac$Frac[i]
    
    if(frac == 1){
      urban.samp <- ur.draws[,paste0(adm.name,':U')]
      aggre.samp[,i] <- urban.samp
      
    }else{
      
      if(frac == 0){
        rural.samp <- ur.draws[,paste0(adm.name,':R')]
        aggre.samp[,i]  <- rural.samp
        
      }else{
        
        urban.samp <- ur.draws[,paste0(adm.name,':U')]
        rural.samp <- ur.draws[,paste0(adm.name,':R')]
        
        aggre.samp[,i]  <- urban.samp*frac+ rural.samp*(1-frac)
      }
      
    }
    
  }
  
  all.q = matrix(0, nrow = n_admin, ncol = 3)
  
  for(i in 1:n_admin){
    
    all.q[i,] = quantile(aggre.samp[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    
  }
  
  overall.res = data.frame(admin.char = urb.frac$adm.char,
                           admin.name = urb.frac$admin.name,
                           p_Low = all.q[,1],
                           p_Med = all.q[,2],
                           p_Upp = all.q[,3],
                           urb.frac =  urb.frac$Frac)
  
  }else{
  
  n_admin <- dim(admin.ref)[1]  
  urb.frac  <- urb.frac[,admin.ref$Internal]
    
  aggre.samp <- matrix(0, nrow = nSamp, ncol = n_admin)
  
  for (i in 1:n_admin){
    
    adm.name <- admin.ref$GADM[i]
    frac <- urb.frac[,i]
    
    
        
        urban.samp <- ur.draws[,paste0(adm.name,':U')]
        rural.samp <- ur.draws[,paste0(adm.name,':R')]
        
        aggre.samp[,i]  <- urban.samp*frac+ rural.samp*(1-frac)
      
      
    }
    
  all.q = matrix(0, nrow = n_admin, ncol = 3)
  frac.vec <- vector()
  for(i in 1:n_admin){
    
    all.q[i,] = quantile(aggre.samp[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    frac.vec[i] = median(urb.frac[,i])
  }
  
  overall.res = data.frame(admin.char = admin.ref$Internal,
                           admin.name = admin.ref$GADM,
                           p_Low = all.q[,1],
                           p_Med = all.q[,2],
                           p_Upp = all.q[,3],
                           urb.frac =  frac.vec)
  
  }
  

  
  row.names(overall.res) <- NULL
  
  colnames(aggre.samp) <- admin.ref$Internal

  return(list(Combined.est = overall.res,
              Combined.draws = aggre.samp ))
  
  

 # return (overall.res)
  
  
  
  
}




################################################################
#########  Aggregate draws of estimates, fixed/sampled fractions 
################################################################

#' @u.draws Urban specific draws at native scale, dim: n_samp*n_admin, ordered by admin index
#' @r.draws Rural specific draws at native scale, dim: n_samp*n_admin, ordered by admin index
#' @fixed.urb Indicator whether is fixed 
#' @nSamp number of samples 
#' @admin.ref reference table for admin names and admin index
aggre_Adm_samp  = function(u.draws,r.draws,urb.frac,fixed.urb=T,nSamp=1000,
                         admin.ref){
  
  n_admin <- dim(admin.ref)[1]  
  
  ## fixed urban proportions
  if(fixed.urb){
    
    aggre.samp <- matrix(0, nrow = nSamp, ncol = n_admin)
    
    for (i in 1:n_admin){
      
      adm.name <- admin.ref$GADM[i]
      frac <- urb.frac$Frac[i]
      
      if(frac == 1){
        urban.samp <- u.draws[,i]
        aggre.samp[,i] <- urban.samp
        
      }else{
        
        if(frac == 0){
          rural.samp <- r.draws[,i]
          aggre.samp[,i]  <- rural.samp
          
        }else{
          
          urban.samp <- u.draws[,i]
          rural.samp <- r.draws[,i]
          
          aggre.samp[,i]  <- urban.samp*frac+ rural.samp*(1-frac)
        }
        
      }
      
    }
    
    all.q = matrix(0, nrow = n_admin, ncol = 3)
    
    for(i in 1:n_admin){
      
      all.q[i,] = quantile(aggre.samp[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      
    }
    
    overall.res = data.frame(admin.char = urb.frac$adm.char,
                             admin.name = urb.frac$admin.name,
                             p_Low = all.q[,1],
                             p_Med = all.q[,2],
                             p_Upp = all.q[,3],
                             urb.frac =  urb.frac$Frac)
    
  }else{
    
    urb.frac  <- urb.frac[,admin.ref$Internal]
    
    aggre.samp <- matrix(0, nrow = nSamp, ncol = n_admin)
    
    for (i in 1:n_admin){
      
      adm.name <- admin.ref$GADM[i]
      frac <- urb.frac[,i]
      
      
      
      urban.samp <- u.draws[,i]
      rural.samp <- r.draws[,i]
      
      aggre.samp[,i]  <- urban.samp*frac+ rural.samp*(1-frac)
      
      
    }
    
    all.q = matrix(0, nrow = n_admin, ncol = 3)
    frac.vec <- vector()
    for(i in 1:n_admin){
      
      all.q[i,] = quantile(aggre.samp[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      frac.vec[i] = median(urb.frac[,i])
    }
    
    overall.res = data.frame(admin.char = admin.ref$Internal,
                             admin.name = admin.ref$GADM,
                             p_Low = all.q[,1],
                             p_Med = all.q[,2],
                             p_Upp = all.q[,3],
                             urb.frac =  frac.vec)
    
  }
  
  
  # Calculate quantiles
  U.q = matrix(0, nrow = n_admin, ncol = 3)
  R.q = matrix(0, nrow = n_admin, ncol = 3)
  
  for(i in 1:n_admin){
    U.q[i,] = quantile(u.draws[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    R.q[i,] = quantile(r.draws[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }

  Urban.res = data.frame(admin.char = admin.ref$Internal,
                         admin.name = admin.ref$GADM,
                         p_Low = U.q[,1],
                         p_Med = U.q[,2],
                         p_Upp = U.q[,3])
  
  Rural.res = data.frame(admin.char = admin.ref$Internal,
                         admin.name = admin.ref$GADM,
                         p_Low = R.q[,1],
                         p_Med = R.q[,2],
                         p_Upp = R.q[,3])
  
  row.names(Urban.res) <- NULL
  row.names(Rural.res) <- NULL
  row.names(overall.res) <- NULL
  
  colnames(aggre.samp) <- admin.ref$Internal
  colnames(u.draws) <- admin.ref$Internal
  colnames(r.draws) <- admin.ref$Internal
  
  return(list(Combined.est = overall.res,
              Urban.res = Urban.res,
              Rural.res = Rural.res,
              Combined.draws = aggre.samp,
              Urban.draws = u.draws,
              Rural.draws = r.draws ))
  
  
  
  #row.names(overall.res) <- NULL
  
  #return (overall.res)
  
  
  
  
}


################################################################
#########  Smoothed Direct estimates 
################################################################
# Compute smoothed direct


################################################################
#########  national fixed effects cluster level model
################################################################

# Compute national nonstratified cluster level model
FitNatlNonstrat = function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                               param = c(1, 0.05)))){
  
  
  formula = HIV ~ 1+ f(clusterIdx,
                           model = "iid",
                           hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}


# Compute national stratified cluster level model
FitNatlStrat = function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                            param = c(1, 0.05)))){
  
  
  formula = HIV ~ urban + f(clusterIdx,
                                model = "iid",
                                hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}


# Aggregate national nonstratified cluster level model
aggNatlNonstrat = function(res.inla, nSamp = 1000){
  
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  res.vec = vector()
  
  # Get indicies in the sample
  intIdx  = res.inla$misc$configs$contents$start[3]
  
  # Integrate out random effects
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){
    cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
    res.vec[i] = logit(mean(expit(post.sample[[i]]$latent[intIdx] + nugSimStd*cSD)))
    
  }
  
  
  res.vec.q = quantile(res.vec, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  res.vec.q = 1/(1+exp(-res.vec.q))
  
  etaNatl.overD = data.frame(p_Low = res.vec.q[1],
                             p_Med = res.vec.q[2],
                             p_Upp = res.vec.q[3])
  
  row.names(etaNatl.overD) <- NULL
  
  return(list(Natl.direct= etaNatl.overD,
              samples = res.vec))
}

aggNatlStrat = function(res.inla, nSamp = 1000,natl.frac){
  
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  res.U = vector()
  res.R = vector()
  res.all = vector()
  
  # Get indicies in the sample
  intIdx  = res.inla$misc$configs$contents$start[3]
  urbIdx  = res.inla$misc$configs$contents$start[4]
  
  # Integrate out random effects
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){
    
    etaUrb.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[urbIdx]
    etaRur.tmp = post.sample[[i]]$latent[intIdx]
    
    # Nugget is overdispersion
    cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
    res.U[i] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
    res.R[i] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))
    
  }
  
  # Calculate combined samples
  res.all = logit(natl.frac*expit(res.U) + (1-natl.frac)*expit(res.R))
  
  
  res.U.q = quantile(res.U, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  res.R.q = quantile(res.R, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  res.all.q = quantile(res.all, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  
  res.U.q = 1/(1+exp(-res.U.q))
  res.R.q = 1/(1+exp(-res.R.q))
  res.all.q = 1/(1+exp(-res.all.q))
  
  Urban.res = data.frame(p_Low = res.U.q[1],
                         p_Med = res.U.q[2],
                         p_Upp = res.U.q[3])
  
  Rural.res = data.frame(p_Low = res.R.q[1],
                         p_Med = res.R.q[2],
                         p_Upp = res.R.q[3])
  
  combined.res = data.frame(p_Low = res.all.q[1],
                            p_Med = res.all.q[2],
                            p_Upp = res.all.q[3])
  
  row.names(Urban.res) <- NULL
  row.names(Rural.res) <- NULL
  row.names(combined.res) <- NULL
  
  return(list(Combined.est = combined.res,
              Urban.res = Urban.res,
              Rural.res = Rural.res,
              samples = list(urban.samp = 1/(1+exp(-res.U)),
                             rural.samp = 1/(1+exp(-res.R)),
                             combined.samp = 1/(1+exp(-res.all)))
              
  ))
}



################################################################
#########  admin-2 fixed effects cluster level model
################################################################

### nonstratified
FitAdmin2Nonstrat = function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                                 param = c(1, 0.05)))){
  
  
  formula = HIV ~  -1 + admin2.char +
    f(clusterIdx,
      model = "iid",
      hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  #control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}


aggAdmin2Nonstrat = function(res.inla, admin.ref, nSamp = 1000){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  areaIdx  = res.inla$misc$configs$contents$start[seq(3,2+n_admin)]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){
    for(j in 1:n_admin){
      # Nugget is measurement error
      est.area.temp = post.sample[[i]]$latent[areaIdx[j]] 
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      est.admin[i,j] = logit(mean(expit(est.area.temp + nugSimStd*cSD)))
    }
  }
  
  
  # Calculate quantiles
  etaQuant = matrix(0, nrow = n_admin, ncol = 3)
  
  for(i in 1:n_admin){
    etaQuant[i,] = quantile(est.admin[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  
  pAdmin2 = 1/(1+exp(-etaQuant))
  
  est.admin <- expit(est.admin)
  colnames(est.admin) <- admin.ref$Internal
  admin2.res = list(Combined.est=data.frame(admin.ref$Internal,
                          admin.name = admin.ref$GADM,
                          p_Low = pAdmin2[, 1],
                          p_Med = pAdmin2[, 2],
                          p_Upp = pAdmin2[, 3]),Combined.draws = est.admin)
  
  
  return(admin2.res)
}

### stratified without interaction
FitAdmin2Strat = function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                              param = c(1, 0.05)))){
  
  
  formula = HIV ~  -1 + admin2.char + urban +
    f(clusterIdx,
      model = "iid",
      hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  #control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}


aggAdmin2Strat = function(res.inla, nSamp = 1000,urb.frac,
                          fixed.urb=T,admin.ref=NULL){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.R =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  areaIdx  = res.inla$misc$configs$contents$start[seq(3,2+n_admin)]
  urbIdx  = res.inla$misc$configs$contents$start[3+n_admin]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){
    for(j in 1:n_admin){
      # Nugget is measurement error
      etaUrb.tmp = post.sample[[i]]$latent[areaIdx[j]] + post.sample[[i]]$latent[urbIdx]
      etaRur.tmp = post.sample[[i]]$latent[areaIdx[j]]
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
  }
  
  # calculate combined results 
  combined.res <- aggre_Adm_samp(u.draws=expit(est.U),r.draws=expit(est.R),urb.frac=urb.frac,
                                 fixed.urb=fixed.urb,nSamp=1000,admin.ref=admin.ref)
  
  return(combined.res)
  
 
  
  
}

### stratified with interaction

FitAdmin2StratInt = function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                              param = c(1, 0.05)))){
  
  
  formula = HIV ~  -1 + adm2_ur  +
    f(clusterIdx,
      model = "iid",
      hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  #control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE),
                  control.inla = list(int.strategy = "grid"))
  
  return(res.inla)
}





aggAdmin2StratInt = function(res.inla, nSamp = 1000,urb.frac,
                             fixed.urb=T,admin.ref=NULL){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.R =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  strata.Idx  = res.inla$misc$configs$contents$start[seq(3,2+2*n_admin)]

  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){

    for(j in 1:n_admin){
      # Nugget is measurement error
      etaRur.tmp = post.sample[[i]]$latent[strata.Idx[j*2-1]] 
      etaUrb.tmp = post.sample[[i]]$latent[strata.Idx[j*2]]
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
  }
  
  # calculate combined results 
  combined.res <- aggre_Adm_samp(u.draws=expit(est.U),r.draws=expit(est.R),urb.frac=urb.frac,
                                 fixed.urb=fixed.urb,nSamp=1000,admin.ref=admin.ref)
  
  return(combined.res)

  
  
}




################################################################
#########  admin-2 spatial iid random effects
################################################################


### nonstratified 
FitAdmin2IIDNonstrat = function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                                 param = c(1, 0.05)))){
  
  
  formula = HIV ~  f(admin2.char,
                       model = "iid",
                       hyper = clustPrior) +
    f(clusterIdx,
      model = "iid",
      hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  #control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE),
                  control.inla = list(int.strategy = "grid"))
  
  return(res.inla)
}


aggAdmin2IIDNonstrat = function(res.inla, nSamp = 1000,admin.ref=NULL){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  areaIdx  = res.inla$misc$configs$contents$start[2]
  IntIdx  = res.inla$misc$configs$contents$start[4]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){
    for(j in 1:n_admin){
      # Nugget is measurement error
      est.area.temp = post.sample[[i]]$latent[areaIdx+j-1] +post.sample[[i]]$latent[IntIdx]
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      est.admin[i,j] = logit(mean(expit(est.area.temp + nugSimStd*cSD)))
    }
  }
  
  
  # Calculate quantiles
  etaQuant = matrix(0, nrow = n_admin, ncol = 3)
  
  for(i in 1:n_admin){
    etaQuant[i,] = quantile(est.admin[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  
  pAdmin2 = 1/(1+exp(-etaQuant))
  
  est.admin <- expit(est.admin)
  colnames(est.admin) <- admin.ref$Internal
  admin2.res = list(Combined.est=data.frame(admin.ref$Internal,
                                            admin.name = admin.ref$GADM,
                                            p_Low = pAdmin2[, 1],
                                            p_Med = pAdmin2[, 2],
                                            p_Upp = pAdmin2[, 3]),Combined.draws = est.admin)
  
  
  
  return(admin2.res)
}


### stratified, same iid random effects for urban/rural

FitAdmin2IIDStrat= function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                                    param = c(1, 0.05)))){
  

  formula = HIV ~ urban +  f(admin2.char,
                               model = "iid",
                               hyper = clustPrior) +
    f(clusterIdx,
      model = "iid",
      hyper = clustPrior) 

  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  #control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE),
                  control.inla = list(int.strategy = "grid"))
  
  return(res.inla)
}


aggAdmin2IIDStrat = function(res.inla, nSamp = 1000,urb.frac,
                             fixed.urb=T,admin.ref=NULL){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.R =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  iid.Idx  = res.inla$misc$configs$contents$start[2]
  int.Idx  = res.inla$misc$configs$contents$start[4]
  urb.Idx  = res.inla$misc$configs$contents$start[5]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){
    for(j in 1:n_admin){
      # Nugget is measurement error
      etaRur.tmp = post.sample[[i]]$latent[iid.Idx+j-1] + post.sample[[i]]$latent[int.Idx]
      etaUrb.tmp = post.sample[[i]]$latent[iid.Idx+j-1] +  post.sample[[i]]$latent[int.Idx] +
        post.sample[[i]]$latent[urb.Idx]
      
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
  }
  
  # calculate combined results 
  combined.res <- aggre_Adm_samp(u.draws=expit(est.U),r.draws=expit(est.R),urb.frac=urb.frac,
                                 fixed.urb=fixed.urb,nSamp=1000,admin.ref=admin.ref)
  
  return(combined.res)
  
  
  
}


### stratified with strata random effects

FitAdmin2IIDStratInt = function(myData, clustPrior=list(prec = list(prior = "pc.prec",
                                                                    param = c(1, 0.05)))){
  
  myData$urban.admin2 <- NA
  myData[myData$urban=='U',]$urban.admin2 <- myData[myData$urban=='U',]$admin2.char
  
  myData$rural.admin2 <- NA
  myData[myData$urban=='R',]$rural.admin2 <- myData[myData$urban=='R',]$admin2.char
  
  formula = HIV ~ urban +  f(urban.admin2,
                       model = "iid",
                       hyper = clustPrior) +
    f(rural.admin2,
      model = "iid",
      hyper = clustPrior) +
    f(clusterIdx,
      model = "iid",
      hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  #control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE),
                  control.inla = list(int.strategy = "grid"))
  
  return(res.inla)
}


aggAdmin2IIDStratInt = function(res.inla, nSamp = 1000,urb.frac,
                                fixed.urb=T,admin.ref=NULL){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.R =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  U.iid.Idx  = res.inla$misc$configs$contents$start[2]
  R.iid.Idx  = res.inla$misc$configs$contents$start[3]
  
  int.Idx  = res.inla$misc$configs$contents$start[5]
  urb.Idx  = res.inla$misc$configs$contents$start[6]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  for(i in 1:nSamp){

    for(j in 1:n_admin){
      # Nugget is measurement error
      etaRur.tmp = post.sample[[i]]$latent[R.iid.Idx+j-1] + post.sample[[i]]$latent[int.Idx]
      etaUrb.tmp = post.sample[[i]]$latent[U.iid.Idx+j-1] +  post.sample[[i]]$latent[int.Idx] +
                   post.sample[[i]]$latent[urb.Idx]
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
  }
  
  # calculate combined results 
  combined.res <- aggre_Adm_samp(u.draws=expit(est.U),r.draws=expit(est.R),urb.frac=urb.frac,
                                 fixed.urb=fixed.urb,nSamp=1000,admin.ref=admin.ref)
  
  return(combined.res)
  
  
  
  
}


################################################################
#########  admin-2 BYM model
################################################################


# Compute areal BYM model
getBYMAdmin2 = function(myData, Admin2Graph, bym2prior, clustPrior, admin2 = FALSE,
                        strata=TRUE,
                        int=FALSE){
  
  if(strata==FALSE){
    formula = HIV ~  f(as.numeric(admin2.char), 
                                model = "bym2",
                                graph = Admin2Graph,
                                scale.model = TRUE) +
      f(clusterIdx,
        model = "iid",
        hyper = clustPrior)
  }else{
  if(int==T){
  # Formula
  myData$urban.admin2 <- NA
  myData[myData$urban=='U',]$urban.admin2 <- myData[myData$urban=='U',]$admin2.char
  
  myData$rural.admin2 <- NA
  myData[myData$urban=='R',]$rural.admin2 <- myData[myData$urban=='R',]$admin2.char
  

  formula = HIV ~ urban + f(as.numeric(urban.admin2), 
                                  model = "bym2",
                                  graph = Admin2Graph,
                                  scale.model = TRUE) +
    f(as.numeric(rural.admin2), 
      model = "bym2",
      graph = Admin2Graph,
      scale.model = TRUE)+
      f(clusterIdx,
        model = "iid",
        hyper = clustPrior)
  }else{
    
    
    formula = HIV ~ urban + f(as.numeric(admin2.char), 
                                model = "bym2",
                                graph = Admin2Graph,
                                scale.model = TRUE) +
      f(clusterIdx,
        model = "iid",
        hyper = clustPrior)
    
  }
  }
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}



aggBYMAdmin2 = function(res.inla, nSamp = 1000,urb.frac,strata=TRUE,int=FALSE,
                        fixed.urb=T,admin.ref=NULL){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.R =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  if(strata==FALSE){
    spat.Idx = res.inla$misc$configs$contents$start[2]
    intIdx  = res.inla$misc$configs$contents$start[4]

    
    nugSimStd = rnorm(1e4, mean = 0, sd = 1)
    
    for(i in 1:nSamp){

      for(j in 1:n_admin){
        # Nugget is measurement error
        eta.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spat.Idx + j-1]

        # Nugget is overdispersion
        cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
        est.admin[i,j] = logit(mean(expit(eta.tmp + nugSimStd*cSD)))
    }
    }
    
    # Calculate quantiles
    all.q = matrix(0, nrow = n_admin, ncol = 3)
    
    for(i in 1:n_admin){
        all.q[i,] = quantile(est.admin[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      
    }
    

    all.q = 1/(1+exp(-all.q))
    
    
    est.admin <- expit(est.admin)
    colnames(est.admin) <- admin.ref$Internal
    admin2.res = list(Combined.est=data.frame(admin.ref$Internal,
                                              admin.name = admin.ref$GADM,
                                              p_Low = all.q[, 1],
                                              p_Med = all.q[, 2],
                                              p_Upp = all.q[, 3]),Combined.draws = est.admin)
    
    
    return(admin2.res)
    

  }else{
  if(int==TRUE){
  spat.u.Idx = res.inla$misc$configs$contents$start[2]
  spat.r.Idx = res.inla$misc$configs$contents$start[3]
  intIdx  = res.inla$misc$configs$contents$start[5]
  urbIdx  = res.inla$misc$configs$contents$start[6]
  }else{
    spatIdx = res.inla$misc$configs$contents$start[2]
    intIdx  = res.inla$misc$configs$contents$start[4]
    urbIdx  = res.inla$misc$configs$contents$start[5]
  }
  }
  
  
  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  
  if(int==TRUE){
  for(i in 1:nSamp){

    for(j in 1:n_admin){
      # Nugget is measurement error
      etaUrb.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spat.u.Idx + j-1]+
        post.sample[[i]]$latent[urbIdx] 
      etaRur.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spat.r.Idx + j-1]
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
  }
  }else {
    nugSimStd = rnorm(1e4, mean = 0, sd = 1)
    for(i in 1:nSamp){

      for(j in 1:n_admin){
        # Nugget is measurement error
        etaUrb.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spatIdx + j-1]+
          post.sample[[i]]$latent[urbIdx] 
        etaRur.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spatIdx + j-1]
        
        # Nugget is overdispersion
        cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
        est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
        est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
    }
    
  }
  
  # calculate combined results 
  combined.res <- aggre_Adm_samp(u.draws=expit(est.U),r.draws=expit(est.R),urb.frac=urb.frac,
                                 fixed.urb=fixed.urb,nSamp=1000,admin.ref=admin.ref)
  
  return(combined.res)
  
  
}









################################################################
#########  admin-3 BYM2 model
################################################################








# Compute areal BYM model
getBYMAdmin3 = function(myData, Admin3Graph, bym2prior, clustPrior, 
                        strata=TRUE,
                        int=FALSE){
  
  if(strata==FALSE){
    formula = HIV ~  f(as.numeric(admin3Fac), 
                       model = "bym2",
                       graph = Admin3Graph,
                       scale.model = TRUE) +
      f(clusterIdx,
        model = "iid",
        hyper = clustPrior)
  }else{
    if(int==T){
      # Formula
      myData$urban.admin3 <- NA
      myData[myData$urban=='U',]$urban.admin3 <- myData[myData$urban=='U',]$admin3Fac
      
      myData$rural.admin3 <- NA
      myData[myData$urban=='R',]$rural.admin3 <- myData[myData$urban=='R',]$admin3Fac
      
      
      formula = HIV ~ urban + f(as.numeric(urban.admin3), 
                                model = "bym2",
                                graph = Admin3Graph,
                                scale.model = TRUE) +
        f(as.numeric(rural.admin3), 
          model = "bym2",
          graph = Admin3Graph,
          scale.model = TRUE)+
        f(clusterIdx,
          model = "iid",
          hyper = clustPrior)
    }else{
      
      
      formula = HIV ~ urban + f(as.numeric(admin3Fac), 
                                model = "bym2",
                                graph = Admin3Graph,
                                scale.model = TRUE) +
        f(clusterIdx,
          model = "iid",
          hyper = clustPrior)
      
    }
  }
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE,dic=TRUE,cpo=TRUE,waic=TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}




aggBYMAdmin3 = function(res.inla, nSamp = 1000,urb.frac,strata=TRUE,int=FALSE,
                        fixed.urb=T,admin.ref=NULL){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_admin = dim(admin.ref)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.R =  matrix(0, nrow = nSamp, ncol = n_admin)
  est.admin = matrix(0, nrow = nSamp, ncol = n_admin)
  
  # Get indicies in the sample
  if(strata==FALSE){
    spat.Idx = res.inla$misc$configs$contents$start[2]
    intIdx  = res.inla$misc$configs$contents$start[4]
    
    
    nugSimStd = rnorm(1e4, mean = 0, sd = 1)
    
    for(i in 1:nSamp){
      
      for(j in 1:n_admin){
        # Nugget is measurement error
        eta.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spat.Idx + j-1]
        
        # Nugget is overdispersion
        cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
        est.admin[i,j] = logit(mean(expit(eta.tmp + nugSimStd*cSD)))
      }
    }
    
    # Calculate quantiles
    all.q = matrix(0, nrow = n_admin, ncol = 3)
    
    for(i in 1:n_admin){
      all.q[i,] = quantile(est.admin[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      
    }
    
    
    all.q = 1/(1+exp(-all.q))
    
    all.res = data.frame(admin.char = admin.ref$Internal,
                         admin.name = admin.ref$GADM,
                         p_Low = all.q[,1],
                         p_Med = all.q[,2],
                         p_Upp = all.q[,3])
    
    return(all.res)
    
  }else{
    if(int==TRUE){
      spat.u.Idx = res.inla$misc$configs$contents$start[2]
      spat.r.Idx = res.inla$misc$configs$contents$start[3]
      intIdx  = res.inla$misc$configs$contents$start[5]
      urbIdx  = res.inla$misc$configs$contents$start[6]
    }else{
      spatIdx = res.inla$misc$configs$contents$start[2]
      intIdx  = res.inla$misc$configs$contents$start[4]
      urbIdx  = res.inla$misc$configs$contents$start[5]
    }
  }
  
  
  # Sample urban and rural
  nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  
  if(int==TRUE){
    for(i in 1:nSamp){
      
      for(j in 1:n_admin){
        # Nugget is measurement error
        etaUrb.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spat.u.Idx + j-1]+
          post.sample[[i]]$latent[urbIdx] 
        etaRur.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spat.r.Idx + j-1]
        
        # Nugget is overdispersion
        cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
        est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
        est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
    }
  }else {
    nugSimStd = rnorm(1e4, mean = 0, sd = 1)
    for(i in 1:nSamp){
      
      for(j in 1:n_admin){
        # Nugget is measurement error
        etaUrb.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spatIdx + j-1]+
          post.sample[[i]]$latent[urbIdx] 
        etaRur.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spatIdx + j-1]
        
        # Nugget is overdispersion
        cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
        est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
        est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
    }
    
  }
  
  # calculate combined results 
  combined.res <- aggre_Adm_samp(u.draws=expit(est.U),r.draws=expit(est.R),urb.frac=urb.frac,
                                 fixed.urb=fixed.urb,nSamp=1000,admin.ref=admin.ref)
  
  return(combined.res)
  
  
}

