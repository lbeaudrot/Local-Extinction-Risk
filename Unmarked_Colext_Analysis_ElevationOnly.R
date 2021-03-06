######### R Code for Dynamic Occupancy Modeling #########
######### Model Set 1: Elevation                #########
######### Beaudrot, Acevedo, Lessard et al.     #########

rm(list=ls())
library(unmarked)
library(plyr)
library(AICcmodavg)


######################################
# Define helper functions
######################################

isEmpty <- function(x) {
  return(length(x)==0)
}


CondNum <- function(model){
  max(eigen(hessian(model))$values)/min(eigen(hessian(model))$values)
}


######################################
# Load formatted data
######################################

load("Species7sites_Include.RData")
load("All_covs.RData")

Species_data <- Species7sites_Include

nms=names(Species_data)

# Empty objects for loop
results.all=list()
mods.all=list()
results.table.ma=list()
results.table.aic=list()

############################################################################
################ Begin loop to run all models and extract results ##########
############################################################################


for(k in 1:length(nms)){
  print(k)
  
  # Define species for analysis and site using index value for list of all species
  sp.name <- names(Species_data)[k]
  species <- Species_data[[k]]
  site <- substr(names(Species_data)[[k]],1,3)
  
  
  # Define site covariates
  site_covs <- paste(site, "covs", sep="_")
  covs <- All_covs[names(All_covs)==site_covs]
  Elevation <- unlist(as.matrix(sapply(covs, "[", 1)))
  site.covs<-data.frame(Elevation)
  
  # Define number of primary periods
  yrs <- as.data.frame(sapply(covs, "[", 4))
  to=dim(yrs)[2]
  
  # Create object with data formatted for unmarked
  umf<-unmarkedMultFrame(y=species, siteCovs=site.covs, numPrimary=dim(yrs)[2])
  
  # Create list to hold model set for model selection
  mods=list()
  
  
  ######################################
  # Define & Run Models
  ######################################
  
  # A model is only included in the model set (i.e., mods) if convergence occurs and the condition number is < 2000 (i.e., CondNum < 2000).
  
  # Null model (no covariates) ###################################################################
  try((fm0=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~1,
                  pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm0")) {
    if(CondNum(fm0)<5000){
      if(CondNum(fm0)>0){mods=c(mods,fm0)}
    } 
  }
  
#  try((fm0.1=colext(psiformula=~Elevation,
#                  gammaformula=~1,
#                  epsilonformula=~1,
#                  pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
#  if(exists("fm0.1")) {
#    if(CondNum(fm0.1)<5000){
#      if(CondNum(fm0.1)>0){mods=c(mods,fm0.1)}
#    } 
#  }
  
#  try((fm0.12=colext(psiformula=~Elevation^2,
#                  gammaformula=~1,
#                  epsilonformula=~1,
#                  pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
#  if(exists("fm0.12")) {
#    if(CondNum(fm0.12)<5000){
#      if(CondNum(fm0.12)>0){mods=c(mods,fm0.12)}
#    } 
#  }
  
   # NO COVARIATES FOR INITIAL OCCUPANCY
   # Elevation as a covariate of colonization and extinction ######################################
  try((fm1=colext(psiformula=~1,
                  gammaformula=~Elevation,
                  epsilonformula=~Elevation,
                  pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm1")) {
    if(CondNum(fm1)<5000){
      if(CondNum(fm1)>0){mods=c(mods,fm1)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm1.1=colext(psiformula=~1,
                    gammaformula=~Elevation,
                    epsilonformula=~1,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm1.1")) {
    if(CondNum(fm1.1)<5000){
      if(CondNum(fm1.1)>0){mods=c(mods,fm1.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm1.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Elevation,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm1.2")) {
    if(CondNum(fm1.2)<5000){
      if(CondNum(fm1.2)>0){mods=c(mods,fm1.2)}
    } 
  }
  
  
  # ELEVATION AS A COVARIATE FOR INITIAL OCCUPANCY
  # Elevation as a covariate of colonization and extinction ######################################
  try((fm2=colext(psiformula=~Elevation,
                  gammaformula=~Elevation,
                  epsilonformula=~Elevation,
                  pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2")) {
    if(CondNum(fm2)<5000){
      if(CondNum(fm2)>0){mods=c(mods,fm2)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm2.1=colext(psiformula=~Elevation,
                    gammaformula=~Elevation,
                    epsilonformula=~1,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.1")) {
    if(CondNum(fm2.1)<5000){
      if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm2.2=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Elevation,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.2")) {
    if(CondNum(fm2.2)<5000){
      if(CondNum(fm2.2)>0){mods=c(mods,fm2.2)}
    } 
  }
  
  ################# WITH QUADRATICS ##############################################################

  # Elevation as a covariate of colonization and extinction ######################################
  try((fm12.0=colext(psiformula=~1,
                  gammaformula=~Elevation^2,
                  epsilonformula=~Elevation^2,
                  pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.0")) {
    if(CondNum(fm12.0)<5000){
      if(CondNum(fm12.0)>0){mods=c(mods,fm12.0)}
    } 
  }
  
  try((fm12=colext(psiformula=~Elevation^2,
                  gammaformula=~Elevation^2,
                  epsilonformula=~Elevation^2,
                  pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12")) {
    if(CondNum(fm12)<5000){
      if(CondNum(fm12)>0){mods=c(mods,fm12)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm12.10=colext(psiformula=~1,
                    gammaformula=~Elevation^2,
                    epsilonformula=~1,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.10")) {
    if(CondNum(fm12.10)<5000){
      if(CondNum(fm12.10)>0){mods=c(mods,fm12.10)}
    } 
  }
  
  try((fm12.1=colext(psiformula=~Elevation^2,
                    gammaformula=~Elevation^2,
                    epsilonformula=~1,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.1")) {
    if(CondNum(fm12.1)<5000){
      if(CondNum(fm12.1)>0){mods=c(mods,fm12.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm12.20=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Elevation^2,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.20")) {
    if(CondNum(fm12.20)<5000){
      if(CondNum(fm12.20)>0){mods=c(mods,fm12.20)}
    } 
  }

  try((fm12.2=colext(psiformula=~Elevation^2,
                    gammaformula=~1,
                    epsilonformula=~Elevation^2,
                    pformula=~Elevation,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.2")) {
    if(CondNum(fm12.2)<5000){
      if(CondNum(fm12.2)>0){mods=c(mods,fm12.2)}
    } 
  }
    
  ######################################
  # Run Model Selection
  ######################################
  
  
  # If zero models for a species converged:
  
  if(isEmpty(mods)==TRUE){
    results.all[[k]]=NA  
    results.table.ma[[k]]=NA
    names(results.table.ma)[k] <- nms[k]
    results.table.aic[[k]]=NA
    names(results.table.aic)[k] <- nms[k]
  }else{
    models=fitList(fits=mods)
    (ms <- modSel(models))
    results.all[[k]]=ms    
    
    
    # For species with models that converged, run model selection:
    
    models=fitList(fits=mods)
    
    (ms <- modSel(models))
    
    results.all[[k]] <- ms
    mods.all[[k]] <- mods
    
    
    
    ######################################
    # Extract Results
    ######################################
    
    toExport <- as(ms,"data.frame")
    
    null0.aic <- toExport$delta[toExport$formula=="~1 ~ 1 ~ 1 ~ 1"]
    #null01.aic <- toExport$delta[toExport$formula=="~Elevation ~ 1 ~ 1 ~ 1"]
    #null012.aic <- toExport$delta[toExport$formula=="~Elevation^2 ~ 1 ~ 1 ~ 1"]

    
    #if null didn't converge
    if(isEmpty(null0.aic)==TRUE){
    #if(isEmpty(null0.aic)==TRUE && isEmpty(null01.aic)==TRUE && isEmpty(null012.aic)==TRUE){
      nulls <- NA
    }else{
      nulls <- rbind(toExport[toExport$formula=="~1 ~ 1 ~ 1 ~ 1",])
              #rbind(toExport[toExport$formula=="~1 ~ 1 ~ 1 ~ 1",],
                     #toExport[toExport$formula=="~Elevation ~ 1 ~ 1 ~ 1",],
                     #toExport[toExport$formula=="~Elevation^2 ~ 1 ~ 1 ~ 1",])
    }
    
    
    if((null0.aic==0) || (isEmpty(null0.aic)==TRUE)){
    #if((null0.aic==0) || (null01.aic==0) || (null012.aic==0) || (isEmpty(null0.aic)==TRUE && isEmpty(null01.aic)==TRUE && isEmpty(null012.aic)==TRUE) || (isEmpty(null01.aic)==TRUE && isEmpty(null012.aic)==TRUE)){
      results.table.ma[[k]] <- rbind(nulls)
      results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
      names(results.table.ma)[k] <- nms[k]
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[k]] <- rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
                                #rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",],
                                      #temp[temp$formula=="~Elevation ~ 1 ~ 1 ~ 1",],
                                      #temp[temp$formula=="~Elevation^2 ~ 1 ~ 1 ~ 1",])
      names(results.table.aic)[k] <- nms[k]
      
    }else{
      results.table.ma[[k]] <- rbind(toExport[1,], nulls)
      results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
      names(results.table.ma)[k] <- nms[k]
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[k]] <- rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
                                #rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",],
                                               #temp[temp$formula=="~Elevation ~ 1 ~ 1 ~ 1",],
                                               #temp[temp$formula=="~Elevation^2 ~ 1 ~ 1 ~ 1",])
      results.table.aic[[k]] <- cbind(nms[k], results.table.aic[[k]])
      names(results.table.aic)[k] <- nms[k]
    }
  }    
    
  # Remove all models and results
  rm(fm0, #fm0.1, fm0.12, 
     fm1, fm1.1, fm1.2,
     fm2, fm2.1, fm2.2, 
     fm12, fm12.1, fm12.2,
     fm12.0, fm12.10, fm12.20,
     #fm20, fm22, fm22.1, fm22.2,
     #fm30, fm32, fm32.1, fm32.2,
     mods, ms, temp, toExport,
     null0.aic, nulls)
}

############################################################################
# End loop
############################################################################




######################################
# Write Results
######################################

# Coerce the lists of results into dataframes and write to files
results.table.ma.df <- ldply(results.table.ma, data.frame)
results.table.aic.df <- ldply(results.table.aic, data.frame)
#colext.transformed.df <- ldply(colext.transformed, data.frame)

write.csv(results.table.ma.df, file="results.table.ma_Elevation.csv")
write.csv(results.table.aic.df, file="results.table.aic_Elevation.csv")
#write.csv(colext.transformed.df, file="colext.transformed.csv")


for(i in 1:length(nms)) {
   outputname <- paste(nms[i], "colextAIC.elevation", "csv", sep=".")
   if(is.na(results.all[[i]])==TRUE){
     output <- NA
   }else{
     output <- results.all[[i]]@Full
   } 
write.csv(output, file=outputname)
}