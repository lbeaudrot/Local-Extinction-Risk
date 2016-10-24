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

# Elevation as a covariate of colonization only ################################################
try((fm2.1=colext(psiformula=~1,
                  gammaformula=~Elevation,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.1")) {
  if(CondNum(fm2.1)<2000){
    if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
} 
}

# Elevation as a covariate of extinction only ##################################################
try((fm2.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Elevation,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.2")) {
  if(CondNum(fm2.2)<2000){
    if(CondNum(fm2.2)>0){mods=c(mods,fm2.2)}
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

null.aic <- toExport$delta[toExport$formula=="~1 ~ 1 ~ 1 ~ 1"]

#if null didn't converge
if(isEmpty(null.aic)==TRUE){
	null <- NA
}else{
null <- toExport[toExport$formula=="~1 ~ 1 ~ 1 ~ 1",]
}


if((null.aic==0) ||isEmpty(null.aic)==TRUE){
results.table.ma[[k]] <- rbind(null)
temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
names(temp) <- c("formula","delta","AICwt")
results.table.aic[[k]] <- rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])

}else{
results.table.ma[[k]] <- rbind(toExport[1,],null)
results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
names(results.table.ma)[k] <- nms[k]
temp <- data.frame(toExport$formula,toExport$delta,toExport$AICwt)
names(temp) <- c("formula","delta","AICwt")
results.table.aic[[k]] <- rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
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
rm(fm0, fm2, fm2.1, fm2.2, mods, ms, tmp, temp, toExport)
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

