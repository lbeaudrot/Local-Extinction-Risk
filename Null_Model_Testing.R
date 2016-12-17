######### Null Model Testing for 9 "elevation shift" populations in Nature submission #########


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



#load("reshuffled.col.RData") #reshuffling columns (i.e. time)
#load("reshuffled.all.RData") #reshuffling all values in the matrix
load("reshuffled.row.RData") #reshuffling rows (i.e. space)


#load("Species7sites_Include.RData")
#load("All_covs.RData")

#Species_data <- Species7sites_Include
#nms=names(Species_data)
  
#null_species <- c("BIF.Cephalophus_nigrifrons",#22
#                  "BIF.Pan_troglodytes",#27
#                  "NAK.Tragulus_kanchil",#54
#                  "UDZ.Cricetomys_gambianus",#13
#                  "UDZ.Potamochoerus_larvatus",#19
#                  "VB_.Dasypus_novemcinctus",#3
#                  "VB_.Pecari_tajacu",#6
#                  "YAN.Eira_barbara",#43
#                  "YAN.Mazama_americana")#45

Species_use <- reshuffled.row

#names(Species_use) <- null_species

null_species <- names(Species_use)


# Empty objects for loop
results.all=list()
mods.all=list()
results.table.ma=list()
results.table.ma.k=list()
results.table.aic=list()
results.table.aic.k=list()
colext.transformed=list()


############################################################################
################ Begin loop to run all models and extract results ##########
############################################################################


for(k in 1:length(Species_use)){
  print(k)
  
  for(j in 1:length(Species_use[[k]])){
    #for(j in 1:3){
    print(j)
  
  # Define species for analysis and site using index value for list of all species
  sp.name <- names(Species_use)[k]
  species <- as.matrix(data.frame(Species_use[[k]][j]))
  site <- substr(names(Species_use)[[k]],1,3)
  

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
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm0")) {
    if(CondNum(fm0)<2000){
      if(CondNum(fm0)>0){mods=c(mods,fm0)}
    } 
  }
  
  
  # Elevation as a covariate of colonization and extinction ######################################
  try((fm2=colext(psiformula=~1,
                  gammaformula=~Elevation,
                  epsilonformula=~Elevation,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm2")) {
    if(CondNum(fm2)<2000){
      if(CondNum(fm2)>0){mods=c(mods,fm2)}
    } 
  }
  
  # Elevation as a covariate of colonization only ######################################
  try((fm3=colext(psiformula=~1,
                  gammaformula=~Elevation,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm3")) {
    if(CondNum(fm3)<2000){
      if(CondNum(fm3)>0){mods=c(mods,fm3)}
    } 
  }
  
  # Elevation as a covariate of extinction only ######################################
  try((fm4=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Elevation,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm4")) {
    if(CondNum(fm4)<2000){
      if(CondNum(fm4)>0){mods=c(mods,fm4)}
    } 
  }
  
  
  # ELEVATION AS A COVARIATE FOR INITIAL OCCUPANCY
  # Elevation as a covariate of colonization and extinction ######################################
  try((fm5=colext(psiformula=~Elevation,
                  gammaformula=~Elevation,
                  epsilonformula=~Elevation,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5")) {
    if(CondNum(fm5)<2000){
      if(CondNum(fm5)>0){mods=c(mods,fm5)}
    } 
  }
  
  # Elevation as a covariate of colonization only ################################################
  try((fm5.1=colext(psiformula=~Elevation,
                    gammaformula=~Elevation,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.1")) {
    if(CondNum(fm5.1)<2000){
      if(CondNum(fm5.1)>0){mods=c(mods,fm5.1)}
    } 
  }
  
  # Elevation as a covariate of extinction only ##################################################
  try((fm5.2=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Elevation,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
  
  if(exists("fm5.2")) {
    if(CondNum(fm5.2)<2000){
      if(CondNum(fm5.2)>0){mods=c(mods,fm5.2)}
    } 
  }
  
  ######################################
  # Run Model Selection
  ######################################
  
  if(isEmpty(mods)==TRUE){  # If zero models for a species converged:
    results.all[[j]]=NA  
    results.table.ma[[j]]=NA
    names(results.table.ma)[[j]] <- null_species[k]
    results.table.aic[[j]]=NA
    names(results.table.aic)[[j]] <- null_species[k]
  }else{   # For species with models that converged, run model selection:
    models=fitList(fits=mods)
    (ms <- modSel(models))
    results.all[[j]] <- ms
    mods.all[[j]] <- mods
    
    ######################################
    # Extract Results
    ######################################
    toExport <- as(ms,"data.frame")
    null.aic <- toExport$delta[toExport$formula=="~1 ~ 1 ~ 1 ~ 1"]
    
    #if null didn't converge
    if(isEmpty(null.aic)==TRUE){
      null <- NA
    }else{
      null <- toExport[toExport$formula=="~1 ~ 1 ~ 1 ~ 1",]
    }
    
    if((null.aic==0) ||isEmpty(null.aic)==TRUE){
      results.table.ma[[j]] <- rbind(null)
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[j]] <- rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
      
    }else{
      results.table.ma[[j]] <- rbind(toExport[1,],null)
      results.table.ma[[j]] <- cbind(null_species[k], results.table.ma[[j]])
      #names(results.table.ma)[j] <- null_species[k]
      temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
      names(temp) <- c("formula","delta","AICwt")
      results.table.aic[[j]] <- rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
      results.table.aic[[j]] <- cbind(null_species[k], results.table.aic[[j]])
      #names(results.table.aic)[j] <- null_species[k]
    }
  }
  results.table.ma.k[[k]] <- results.table.ma
  results.table.aic.k[[k]] <- results.table.aic
  
    # Remove all models and results
  rm(fm0, fm2, fm3, fm4, 
     fm5, fm5.1, fm5.2,
     models, mods, ms, temp, toExport, null.aic, null)
  }
}


############################################################################
# End loop
############################################################################

### NEED TO FIGURE OUT HOW TO WRITE RESULTS EFFICIENTLY FOR 62 POPULATIONS
# WITH 250 SIMULATIONS EACH


######################################
# Write Results
######################################

# Coerce the lists of results into dataframes and write to files
results.table.ma.df.k1 <- ldply(results.table.ma.k[[1]], data.frame)
results.table.ma.df.k2 <- ldply(results.table.ma.k[[2]], data.frame)
results.table.ma.df.k3 <- ldply(results.table.ma.k[[3]], data.frame)
results.table.ma.df.k4 <- ldply(results.table.ma.k[[4]], data.frame)
results.table.ma.df.k5 <- ldply(results.table.ma.k[[5]], data.frame)
results.table.ma.df.k6 <- ldply(results.table.ma.k[[6]], data.frame)
results.table.ma.df.k7 <- ldply(results.table.ma.k[[7]], data.frame)
results.table.ma.df.k8 <- ldply(results.table.ma.k[[8]], data.frame)
results.table.ma.df.k9 <- ldply(results.table.ma.k[[9]], data.frame)

results.table.aic.df.k1 <- ldply(results.table.aic.k[[1]], data.frame)
results.table.aic.df.k2 <- ldply(results.table.aic.k[[2]], data.frame)
results.table.aic.df.k3 <- ldply(results.table.aic.k[[3]], data.frame)
results.table.aic.df.k4 <- ldply(results.table.aic.k[[4]], data.frame)
results.table.aic.df.k5 <- ldply(results.table.aic.k[[5]], data.frame)
results.table.aic.df.k6 <- ldply(results.table.aic.k[[6]], data.frame)
results.table.aic.df.k7 <- ldply(results.table.aic.k[[7]], data.frame)
results.table.aic.df.k8 <- ldply(results.table.aic.k[[8]], data.frame)
results.table.aic.df.k9 <- ldply(results.table.aic.k[[9]], data.frame)

write.csv(results.table.ma.df.k1, file="results.table.ma_BIF.Cephalophus_nigrifrons.csv")
write.csv(results.table.ma.df.k2, file="results.table.ma_BIF.Pan_troglodytes.csv")
write.csv(results.table.ma.df.k3, file="results.table.ma_NAK.Tragulus_kanchil.csv")
write.csv(results.table.ma.df.k4, file="results.table.ma_UDZ.Cricetomys_gambianus.csv")
write.csv(results.table.ma.df.k5, file="results.table.ma_UDZ.Potamochoerus_larvatus.csv")
write.csv(results.table.ma.df.k6, file="results.table.ma_VB_.Dasypus_novemcinctus.csv")
write.csv(results.table.ma.df.k7, file="results.table.ma_VB_.Pecari_tajacu.csv")
write.csv(results.table.ma.df.k8, file="results.table.ma_YAN.Eira_barbara.csv")
write.csv(results.table.ma.df.k9, file="results.table.ma_YAN.Mazama_americana.csv")

write.csv(results.table.aic.df.k1, file="results.table.aic_BIF.Cephalophus_nigrifrons.csv")
write.csv(results.table.aic.df.k2, file="results.table.aic_BIF.Pan_troglodytes.csv")
write.csv(results.table.aic.df.k3, file="results.table.aic_NAK.Tragulus_kanchil.csv")
write.csv(results.table.aic.df.k4, file="results.table.aic_UDZ.Cricetomys_gambianus.csv")
write.csv(results.table.aic.df.k5, file="results.table.aic_UDZ.Potamochoerus_larvatus.csv")
write.csv(results.table.aic.df.k6, file="results.table.aic_VB_.Dasypus_novemcinctus.csv")
write.csv(results.table.aic.df.k7, file="results.table.aic_VB_.Pecari_tajacu.csv")
write.csv(results.table.aic.df.k8, file="results.table.aic_YAN.Eira_barbara.csv")
write.csv(results.table.aic.df.k9, file="results.table.aic_YAN.Mazama_americana.csv")

#for(i in 1:length(nms)) {
#  outputname <- paste(nms[i], "colextAIC.elevation", "csv", sep=".")
#  if(is.na(results.all[[i]])==TRUE){
#    output <- NA
#  }else{
#    output <- results.all[[i]]@Full
#  } 
#  write.csv(output, file=outputname)
#}