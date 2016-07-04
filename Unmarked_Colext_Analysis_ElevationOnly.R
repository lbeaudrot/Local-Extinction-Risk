# Unmarked analysis of 62 TEAM populations with >5 detections/year at sites with >500 m elevation change
library(unmarked)
library(plyr)
library(AICcmodavg)

rm(list=ls())
load("All_covs.RData")
load("PCA1_covs.RData")
load("All500m_covariate_species.RData") # 32 populations (less with birds excluded)
load("All_species7sites.RData") # 166 populations
load("Species7sites_Include.RData") # 62 populations (excludes binomial cases)

load("BIOTIC_all.RData") # For 166 populations
#load("BIOTIC_ALL_YEARS.RData") # For 166 populations

load("BIOTIC_Include.RData") # For 62 populations (excludes binomial cases)
load("BIOTIC_ALL_YEARS_Include.RData") # For 62 populations (excludes binomial cases)

#All500m_covariate_species <- All_species7sites
All500m_covariate_species <- Species7sites_Include
BIOTIC_166 <- BIOTIC_Include
#BIOTIC_ALL_YEARS <- BIOTIC_ALL_YEARS_Include

# Matrices for each population are contained in the object "All500m_covariate_species"
nms=names(All500m_covariate_species)

results.all=list()
mods.all=list()

# add 3 new objects from Miguel's script here 
results.table.ma=list() #add at the begining
results.table.aic=list() #add at the begining
colext.transformed=list() #add at the beginning

isEmpty <- function(x) {
    return(length(x)==0)
}


CondNum <- function(model){
  max(eigen(hessian(model))$values)/min(eigen(hessian(model))$values)
}
# produces identical condition number to function "extractCN" from the AICcmodavg package

####
for(k in 45:45){
#for(k in 1:length(nms)){
print(k)

# DEFINE SPECIES for analysis and site USING INDEX VALUE for list of all species (see previous call for list of species names)
index <- k
sp.name <- names(All500m_covariate_species)[index]
sp.name
species <- All500m_covariate_species[[index]]
site <- substr(names(All500m_covariate_species)[[index]],1,3)
site

# Define covariate object based on site

#siteCovs
site_covs <- paste(site, "covs", sep="_")
covs <- All_covs[names(All_covs)==site_covs]
Elevation <- unlist(as.matrix(sapply(covs, "[", 1)))
ForestLossCT <- unlist(as.matrix(sapply(covs, "[", 2)))
ForestGainCT <- unlist(as.matrix(sapply(covs, "[", 3)))
Biotic <- BIOTIC_166[[index]]
pca_site <- paste(site, "pca1", sep="_")
pca_covs <- PCA1_covs[names(PCA1_covs)==pca_site]

#yearlySiteCovs
Tmin <- as.data.frame(sapply(covs, "[", 4))
Tmax <- as.data.frame(sapply(covs, "[", 5))
Tvar <- as.data.frame(sapply(covs, "[", 6))
Tsd <- as.data.frame(sapply(covs, "[", 7))
Tmean <- as.data.frame(sapply(covs, "[", 8))
#BioticYearly <- BIOTIC_ALL_YEARS[[index]]
PCA1 <- data.frame(pca_covs)


to=dim(Tmin)[2]
years=as.character(1:to)
years=matrix(years,nrow(species),to,byrow=TRUE)

# ADD BIOTIC TO UMF COVARIATES
site.covs<-data.frame(Elevation, ForestLossCT, ForestGainCT, Biotic)

umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(Tmin=Tmin,Tmax=Tmax,Tvar=Tvar,Tsd=Tsd,Tmean=Tmean, PCA1=PCA1), siteCovs=site.covs, numPrimary=dim(Tmin)[2])
#umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(year=years,Tmin=Tmin,Tmax=Tmax,Tvar=Tvar,Tsd=Tsd,Tmean=Tmean), siteCovs=site.covs, numPrimary=dim(Tmin)[2])

mods=list()

# Null ##################################################################################
try((fm0=colext(psiformula=~1,
                gammaformula=~1,
                epsilonformula=~1,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm0")) {
  if(CondNum(fm0)<2000){
    if(CondNum(fm0)>0){mods=c(mods,fm0)}
} 
}

# Elevation only ##################################################################################
try((fm2=colext(psiformula=~1,
                gammaformula=~Elevation,
                epsilonformula=~Elevation,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2")) {
  if(CondNum(fm2)<2000){
    if(CondNum(fm2)>0){mods=c(mods,fm2)}
} 
}

try((fm2.1=colext(psiformula=~1,
                  gammaformula=~Elevation,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.1")) {
  if(CondNum(fm2.1)<2000){
    if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
} 
}

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
#Model Selection
######################################

#ifelse(isEmpty(mods)==TRUE, NA, mods)

#if no models converged

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


####

models=fitList(fits=mods)

(ms <- modSel(models))

results.all[[k]]=ms
mods.all[[k]]=mods



#### add to Export to the end (including the bracket)

toExport<-as(ms,"data.frame") #add after ms object

#null.elev.aic=toExport$delta[toExport$formula=="~Elevation ~ 1 ~ 1 ~ 1"]
null.aic=toExport$delta[toExport$formula=="~1 ~ 1 ~ 1 ~ 1"]

#if null.elev didn't converge
#if(isEmpty(null.elev.aic)==TRUE){
#  null.elev=NA  
#}else{
#null.elev=toExport[toExport$formula=="~Elevation ~ 1 ~ 1 ~ 1",]		
#}

#if null didn't converge
if(isEmpty(null.aic)==TRUE){
	null=NA	
}else{
null=toExport[toExport$formula=="~1 ~ 1 ~ 1 ~ 1",]
}


#if((null.elev.aic==0 || isEmpty(null.elev.aic)==TRUE) || (null.aic==0) ||isEmpty(null.aic)==TRUE){	
if((null.aic==0) ||isEmpty(null.aic)==TRUE){  

	
#results.table.ma[[k]]=rbind(null,null.elev)	
results.table.ma[[k]]=rbind(null)  


temp=data.frame(toExport$formula,toExport$delta,toExport$AICwt)
names(temp)=c("formula","delta","AICwt")
#results.table.aic[[k]]=rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",],temp[temp$formula=="~Elevation ~ 1 ~ 1 ~ 1",])}else{
results.table.aic[[k]]=rbind(temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])}else{

#results.table.ma[[k]]=rbind(toExport[1,],null,null.elev)
results.table.ma[[k]]=rbind(toExport[1,],null)
results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])
names(results.table.ma)[k] <- nms[k]

temp=data.frame(toExport$formula,toExport$delta,toExport$AICwt)
names(temp)=c("formula","delta","AICwt")
#results.table.aic[[k]]=rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",],temp[temp$formula=="~Elevation ~ 1 ~ 1 ~ 1",])
results.table.aic[[k]]=rbind(temp[1,],temp[temp$formula=="~1 ~ 1 ~ 1 ~ 1",])
results.table.aic[[k]] <- cbind(nms[k], results.table.aic[[k]])
names(results.table.aic)[k] <- nms[k]
}

test=seq(3,length(toExport)-10,by=2)
tmp=as.numeric(toExport[1,test])

colext.transformed[[k]]=exp(tmp)
colext.transformed[[k]] <- cbind(nms[k], toExport[1,2], colext.transformed[[k]])
names(colext.transformed)[k] <- nms[k]
}

#add tmp, temp and toExport to the rm object 
rm(fm0,fm0.1,fm1,fm1.1,fm1.2,fm2,fm2.1,fm2.2,fm3,fm3.1,fm3.2,
   mods,ms, tmp, temp, toExport)
}

# Need to coerce the lists into dataframes before writing to a files
results.table.ma.df <- ldply(results.table.ma, data.frame)
results.table.aic.df <- ldply(results.table.aic, data.frame)
colext.transformed.df <- ldply(colext.transformed, data.frame)

write.csv(results.table.ma.df, file="results.table.ma.csv")
write.csv(results.table.aic.df, file="results.table.aic.csv")
write.csv(colext.transformed.df, file="colext.transformed.csv")


for(i in 1:length(nms)) {
  outputname <- paste(nms[i], "colextAIC", "csv", sep=".")
  if(is.na(results.all[[i]])==TRUE){
    output <- NA
  }else{
    output <- results.all[[i]]@Full
  } 
  write.csv(output, file=outputname)
} 

