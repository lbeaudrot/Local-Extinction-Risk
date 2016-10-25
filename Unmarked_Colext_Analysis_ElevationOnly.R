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
colext.transformed=list()

# Create "year" object for each site and repeat in a list as many times as we have populations from the site

VB_.year <- list(year=matrix(c("08", "09", "10", "11", "12", "13", "14", "15"), 60, 8, byrow=TRUE)) # VB 7 populations from 2008-2015 60 cameras c("08", "09", "10", "11", "12", "13", "14", "15")
UDZ.year <- list(year=matrix(c("09", "10", "11", "12", "13", "14"), 61, 6, byrow=TRUE)) # UDZ 13 populations from 2009-2014 61 cameras c("09", "10", "11", "12", "13", "14")
BIF.year <- list(year=matrix(c("10", "11", "12", "13", "14"), 60, 5, byrow=TRUE)) # BIF 9 populations from 2010-2013 60 cameras c("10", "11", "12", "13", "14)
PSH.year <- list(year=matrix(c("11", "12", "13", "14", "15"), 30, 5, byrow=TRUE)) # PSH 10 populations from 2011-2015 30 cameras c("11", "12", "13", "14", "15")
YAN.year <- list(year=matrix(c("11", "12", "13", "14"), 61, 4, byrow=TRUE)) # YAN 10 populations from 2011-2014 61 cameras c("11", "12", "13, "14)
NAK.year <- list(year=matrix(c("10", "11", "12", "13", "14", "15"), 60, 6, byrow=TRUE)) # NAK 5 populations from 2010-2015 60 cameras c("10", "11", "12", "13", "14", "15")
RNF.year <- list(year=matrix(c("10", "11", "12", "13", "14"), 60, 5, byrow=TRUE)) # RNF 8 populations from 2010-2014 60 cameras c("10", "11", "12", "13", "14)

years <- c(rep(list(VB_.year), 7), 
          rep(list(UDZ.year), 13),
          rep(list(BIF.year), 9),
          rep(list(PSH.year), 10),
          rep(list(YAN.year), 10),
          rep(list(NAK.year), 5),
          rep(list(RNF.year), 8))


############################################################################
################ Begin loop to run all models and extract results ##########
############################################################################


for(k in 1:length(nms)){
  print(k)
  
  # Define species for analysis and site using index value for list of all species
  sp.name <- names(Species_data)[k]
  species <- Species_data[[k]]
  site <- substr(names(Species_data)[[k]],1,3)
  
  # Define yearly site covariates (i.e. year)
  year <- years[[k]]
  
  # Define site covariates
  site_covs <- paste(site, "covs", sep="_")
  covs <- All_covs[names(All_covs)==site_covs]
  Elevation <- unlist(as.matrix(sapply(covs, "[", 1)))
  site.covs<-data.frame(Elevation)
  
  # Define number of primary periods
  yrs <- as.data.frame(sapply(covs, "[", 4))
  to=dim(yrs)[2]
  
  # Create object with data formatted for unmarked
  umf<-unmarkedMultFrame(y=species, siteCovs=site.covs, numPrimary=dim(yrs)[2], yearlySiteCovs=year)
  
  # Create list to hold model set for model selection
  mods=list()
  
  
  ######################################
  # Define & Run Models
  ######################################
  
  # A model is only included in the model set (i.e., mods) if convergence occurs and the condition number is < 2000 (i.e., CondNum < 2000).
  
  # Null model (no covariates) ###################################################################
  try((fm0=colext(psiformula=~Elevation,
                  gammaformula=~1,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm0")) {
    if(CondNum(fm0)<2000){
      if(CondNum(fm0)>0){mods=c(mods,fm0)}
    } 
  }
  
  # Elevation as a covariate of colonization and extinction ######################################
  try((fm2=colext(psiformula=~Elevation,
                  gammaformula=~Elevation,
                  epsilonformula=~Elevation,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2")) {
    if(CondNum(fm2)<2000){
      if(CondNum(fm2)>0){mods=c(mods,fm2)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm2.1=colext(psiformula=~Elevation,
                    gammaformula=~Elevation,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.1")) {
    if(CondNum(fm2.1)<2000){
      if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm2.2=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Elevation,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2.2")) {
    if(CondNum(fm2.2)<2000){
      if(CondNum(fm2.2)>0){mods=c(mods,fm2.2)}
    } 
  }
  
  ################# WITH QUADRATICS ##############################################################
  # Null model (no covariates) ###################################################################
  try((fm10=colext(psiformula=~Elevation^2,
                  gammaformula=~1,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm10")) {
    if(CondNum(fm10)<2000){
      if(CondNum(fm10)>0){mods=c(mods,fm10)}
    } 
  }
  
  # Elevation as a covariate of colonization and extinction ######################################
  try((fm12=colext(psiformula=~Elevation^2,
                  gammaformula=~Elevation^2,
                  epsilonformula=~Elevation^2,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12")) {
    if(CondNum(fm12)<2000){
      if(CondNum(fm12)>0){mods=c(mods,fm12)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm12.1=colext(psiformula=~Elevation^2,
                    gammaformula=~Elevation^2,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.1")) {
    if(CondNum(fm12.1)<2000){
      if(CondNum(fm12.1)>0){mods=c(mods,fm12.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm12.2=colext(psiformula=~Elevation^2,
                    gammaformula=~1,
                    epsilonformula=~Elevation^2,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm12.2")) {
    if(CondNum(fm12.2)<2000){
      if(CondNum(fm12.2)>0){mods=c(mods,fm12.2)}
    } 
  }
  
  ################# WITH ANNUAL COL/EXT ESTIMATES ##############################################################
  # Null model (no covariates) ###################################################################
  try((fm20=colext(psiformula=~Elevation,
                   gammaformula=~year-1,
                   epsilonformula=~year-1,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm20")) {
    if(CondNum(fm20)<2000){
      if(CondNum(fm20)>0){mods=c(mods,fm20)}
    } 
  }
  
  # Elevation as a covariate of colonization and extinction ######################################
  try((fm22=colext(psiformula=~Elevation,
                   gammaformula=~year-1 + Elevation,
                   epsilonformula=~year-1 + Elevation,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22")) {
    if(CondNum(fm22)<2000){
      if(CondNum(fm22)>0){mods=c(mods,fm22)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm22.1=colext(psiformula=~Elevation,
                     gammaformula=~year-1 + Elevation,
                     epsilonformula=~year-1,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22.1")) {
    if(CondNum(fm22.1)<2000){
      if(CondNum(fm22.1)>0){mods=c(mods,fm22.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm22.2=colext(psiformula=~Elevation,
                     gammaformula=~year-1,
                     epsilonformula=~year-1 + Elevation,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm22.2")) {
    if(CondNum(fm22.2)<2000){
      if(CondNum(fm22.2)>0){mods=c(mods,fm22.2)}
    } 
  }
  

  ################# WITH QUADRATICS AND ANNUAL COL/EXT ESTIMATES ##############################################################
  # Null model (no covariates) ###################################################################
  try((fm30=colext(psiformula=~Elevation^2,
                   gammaformula=~year-1,
                   epsilonformula=~year-1,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm30")) {
    if(CondNum(fm30)<2000){
      if(CondNum(fm30)>0){mods=c(mods,fm30)}
    } 
  }
  
  # Elevation as a covariate of colonization and extinction ######################################
  try((fm32=colext(psiformula=~Elevation^2,
                   gammaformula=~year-1 + Elevation^2,
                   epsilonformula=~year-1 + Elevation^2,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32")) {
    if(CondNum(fm32)<2000){
      if(CondNum(fm32)>0){mods=c(mods,fm32)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm32.1=colext(psiformula=~Elevation^2,
                     gammaformula=~year-1 + Elevation^2,
                     epsilonformula=~year-1,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32.1")) {
    if(CondNum(fm32.1)<2000){
      if(CondNum(fm32.1)>0){mods=c(mods,fm32.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm32.2=colext(psiformula=~Elevation^2,
                     gammaformula=~year-1,
                     epsilonformula=~year-1 + Elevation^2,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm32.2")) {
    if(CondNum(fm32.2)<2000){
      if(CondNum(fm32.2)>0){mods=c(mods,fm32.2)}
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
    colext.transformed[[k]]=NA
    names(colext.transformed)[k] <- nms[k]
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
    
    null.aic <- toExport$delta[toExport$formula=="~Elevation ~ 1 ~ 1 ~ 1"]
    
    #if null didn't converge
    if(isEmpty(null.aic)==TRUE){
      null <- NA
    }else{
      null <- toExport[toExport$formula=="~Elevation ~ 1 ~ 1 ~ 1",]
    }
    
    
    if((null.aic==0) ||isEmpty(null.aic)==TRUE){
      results.table.ma[[k]] <- rbind(null)
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[k]] <- rbind(temp[temp$formula=="~Elevation ~ 1 ~ 1 ~ 1",])
      
    }else{
      results.table.ma[[k]] <- rbind(toExport[1,],null)
      results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
      names(results.table.ma)[k] <- nms[k]
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[k]] <- rbind(temp[1,],temp[temp$formula=="~Elevation ~ 1 ~ 1 ~ 1",])
      results.table.aic[[k]] <- cbind(nms[k], results.table.aic[[k]])
      names(results.table.aic)[k] <- nms[k]
    }
    
    test <- seq(3,length(toExport)-10,by=2)
    tmp <- as.numeric(toExport[1,test])
    colext.transformed[[k]] <- exp(tmp)
    colext.transformed[[k]] <- cbind(nms[k], toExport[1,2], colext.transformed[[k]])
    names(colext.transformed)[k] <- nms[k]
  }
  
  # Remove all models and results
  rm(fm0, fm2, fm2.1, fm2.2, 
     fm10, fm12, fm12.1, fm12.2,
     fm20, fm22, fm22.1, fm22.2,
     fm30, fm32, fm32.1, fm32.2,
     mods, ms, tmp, temp, toExport)
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
colext.transformed.df <- ldply(colext.transformed, data.frame)

write.csv(results.table.ma.df, file="results.table.ma.csv")
write.csv(results.table.aic.df, file="results.table.aic.csv")
write.csv(colext.transformed.df, file="colext.transformed.csv")

length(nms)
for(i in 1:length(nms)) {
  outputname <- paste(nms[i], "colextAIC", "csv", sep=".")
  output <- results.all[[i]]@Full
  write.csv(output, file=outputname)
} 