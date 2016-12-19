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


load("reshuffled.row.RData") #reshuffling rows (i.e. space)

Species_use <- reshuffled.row

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

# FIRST SAVE LOOP OUTPUT AS Rdata OBJECT

save(results.table.ma.k, file="results.table.ma_250shuffles.RData")
save(results.table.aic.k, file="results.table.aic_250shuffles.RData")

### THEN FIGURE OUT HOW TO WRITE SUMMARY RESULTS EFFICIENTLY FOR 62 POPULATIONS WITH 250 SIMULATIONS EACH
# Automate this process for all 62 populations
# Extract only top row formula and second row delta (if applicable) from each element in the list, convert into a dataframe
# Then select only the rows where actual top model occurs and if null=0 or if elevation model delta >2
# take the length/dimensions of the output


pops <- list()
hold1 <- list()
hold2 <- data.frame()
outs <- vector()
null.df <- matrix(NA, nrow=250, ncol=2, dimnames=list(NULL,cn))
ifempty <- as.vector(c(NA, NA))
pops.df.list <- list()

for(p in 1:length(results.table.aic.k)){
  for(s in 1:length(results.table.aic.k[[p]])){
    hold1 <- results.table.aic.k[[p]][s]
    hold2 <- ldply(hold1)
    outs <- cbind(as.character(hold2$formula[1]), ifelse(hold2$formula[1]=="~1 ~ 1 ~ 1 ~ 1", 0, hold2$delta[2]))
    outs <- cbind(ifelse(isEmpty(outs)==TRUE, ifempty, outs[1]), ifelse(isEmpty(outs)==TRUE, ifempty, outs[2]))
    #outs <- ifelse(dim(outs)[1]==0, ifempty, outs)
    null.df[s,] <- outs
  }
  pops[[p]] <- null.df
  pops.df.list[[p]] <- as.matrix(pops[[p]], nrow=250, ncol=2)
  pops.df.list[[p]] <- as.data.frame(pops.df.list[[p]])
  pops.df.list[[p]]$delta <- as.numeric(as.character(pops.df.list[[p]]$delta)) # Creates a data fame of 250 results for each population
}


# Then select only the rows where actual top model occurs and if null=0 or if elevation model delta >2
# take the length/dimensions of the output


pops.significant <- list()
pops.sig.table <- list()

for(p in 1:length(pops.df.list)){
  #pops.table[[p]] <- table(pops.df.list[[p]]$formula)
  pops.significant[[p]] <- pops.df.list[[p]][pops.df.list[[p]]$delta>2,]
  pops.sig.table[[p]] <- table(pops.significant[[p]]$formula)
}
names(pops.sig.table) <- nms

pops.sig.table.df <- ldply(pops.sig.table, data.frame)
sig.per <- (pops.sig.table.df$Freq/250)*100
sig.flag <- ifelse(sig.per>5, 1, 0)
pops.sig.table.df <- data.frame(pops.sig.table.df, sig.per, sig.flag)

null.flagged <- pops.sig.table.df[pops.sig.table.df$sig.flag==1,]
write.csv(null.flagged, file="null.flagged.csv")





