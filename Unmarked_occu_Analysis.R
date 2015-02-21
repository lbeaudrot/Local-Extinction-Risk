# Unmarked analysis of 33 TEAM populations with >8% detection rates at sites with >500 m elevation change
library(unmarked)
library(plyr)

#load('/Volumes/SCIENCEWORK/Working_folder/UPR_Prof/Collaborations/TEAM/33spp/All_covs_scaled.RData')
#load('/Volumes/SCIENCEWORK/Working_folder/UPR_Prof/Collaborations/TEAM/32spp/All500m_covariate_species.RData')

rm(list=ls())
load("All_covs.RData")
load("PCA1_covs.RData")
load("All500m_covariate_species.RData") # 32 populations (less with birds excluded)
load("All_species7sites.RData") # 166 populations
load("Species7sites_Include.RData") # 62 populations (excludes binomial cases)

load("BIOTIC_all.RData") # For 166 populations
load("BIOTIC_ALL_YEARS.RData") # For 166 populations

load("BIOTIC_Include.RData") # For 62 populations (excludes binomial cases)
load("BIOTIC_ALL_YEARS_Include.RData") # For 62 populations (excludes binomial cases)

#All500m_covariate_species <- All_species7sites
All500m_covariate_species <- Species7sites_Include
BIOTIC_166 <- BIOTIC_Include
BIOTIC_ALL_YEARS <- BIOTIC_ALL_YEARS_Include

# Matrices for each population are contained in the object "All500m_covariate_species"
nms=names(All500m_covariate_species)

results.all=list()
mods.all=list()

# add 3 new objects from Miguel's script here 
results.table.ma=list() #add at the begining
results.table.aic=list() #add at the begining
occu.transformed=list() #add at the beginning

isEmpty <- function(x) {
    return(length(x)==0)
}

CondNum <- function(model){
  max(eigen(hessian(model))$values)/min(eigen(hessian(model))$values)
}


####
#for(k in 50:54){
for(k in 1:length(nms)){
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
#pca1 <- pca_covs[names(pca_covs)==pca_site]


#yearlySiteCovs # Subset appropriate year
Tmin <- as.data.frame(sapply(covs, "[", 4))
Tmax <- as.data.frame(sapply(covs, "[", 5))
Tvar <- as.data.frame(sapply(covs, "[", 6))
Tsd <- as.data.frame(sapply(covs, "[", 7))
Tmean <- as.data.frame(sapply(covs, "[", 8))
BioticYearly <- BIOTIC_ALL_YEARS[[index]]



# Subset appropriate year of species data
species <- species[,1:15]

# Subset appropriate year of temperature data
Tmin <- Tmin[ ,1]
Tmax <- Tmax[ ,1]
Tvar <- Tvar[ ,1]
Tsd <- Tsd[ ,1]
Tmean <- Tmean[ ,1]
BioticYearly <- BioticYearly[ ,1]
PCA1 <- unlist(as.matrix(sapply(pca_covs, "[", 1)))

# ADD BIOTIC TO UMF COVARIATES
site.covs<-data.frame(Elevation, ForestLossCT, ForestGainCT, Biotic, BioticYearly, Tmin, Tmax, Tvar, Tsd, Tmean, PCA1)

umf<-unmarkedFrameOccu(y=species, siteCovs=site.covs, obsCovs=NULL)
#umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(year=years,Tmin=Tmin,Tmax=Tmax,Tvar=Tvar,Tsd=Tsd,Tmean=Tmean), siteCovs=site.covs, numPrimary=dim(Tmin)[3])

mods=list()

#Null##################################################################################
try((fm0=occu(~ 1 ~ 1, 
              data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm0")){
  if(CondNum(fm0)<2000){
    if(CondNum(fm0)>0){mods=c(mods,fm0)}
}
}

#Elevation##################################################################################
try((fm2=occu(~ 1 ~ Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2")){
  if(CondNum(fm2)<2000){
    if(CondNum(fm2)>0){mods=c(mods,fm2)}
}
}

#PCA1##################################################################################
try((fm3=occu(~ 1 ~ PCA1, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3")){
  if(CondNum(fm3)<2000){
    if(CondNum(fm3)>0){mods=c(mods,fm3)}
}
}

#Biotic##################################################################################
try((fm4=occu(~ 1 ~ Biotic, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
if(exists("fm4")){
  if(CondNum(fm4)<2000){
    if(CondNum(fm4)>0){mods=c(mods,fm4)}
}
}

#BioticYearly##################################################################################
try((fm5=occu(~ 1 ~ BioticYearly, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5")){
  if(CondNum(fm5)<2000){
    if(CondNum(fm5)>0){mods=c(mods,fm5)}
}
}

#PCA1 + Elevation ##################################################################################
try((fm5.1=occu(~ 1 ~ PCA1 + Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5.1")){
  if(CondNum(fm5.1)<2000){
    if(CondNum(fm5.1)>0){mods=c(mods,fm5.1)}
}
}

#PCA1 * Elevation ##################################################################################
try((fm6=occu(~ 1 ~ PCA1 * Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6")){
  if(CondNum(fm6)<2000){
    if(CondNum(fm6)>0){mods=c(mods,fm6)}
}
}

#Biotic + Elevation##################################################################################
try((fm6.1=occu(~ 1 ~ Biotic + Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6.1")){
  if(CondNum(fm6.1)<2000){
    if(CondNum(fm6.1)>0){mods=c(mods,fm6.1)}
}
}

#Biotic * Elevation##################################################################################
try((fm7=occu(~ 1 ~ Biotic * Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7")){
  if(CondNum(fm7)<2000){
    if(CondNum(fm7)>0){mods=c(mods,fm7)}
}
}

#BioticYearly + Elevation##################################################################################

try((fm7.1=occu(~ 1 ~ BioticYearly + Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.1")){
  if(CondNum(fm7.1)<2000){
    if(CondNum(fm7.1)>0){mods=c(mods,fm7.1)}
}
}

#BioticYearly * Elevation############################################################################
try((fm8=occu(~ 1 ~ BioticYearly * Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8")){
  if(CondNum(fm8)<2000){
    if(CondNum(fm8)>0){mods=c(mods,fm8)}
}
}

#Biotic + PCA1##########################################################################
try((fm9=occu(~ 1 ~ Biotic + PCA1, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9")){
  if(CondNum(fm9)<2000){
    if(CondNum(fm9)>0){mods=c(mods,fm9)}
}
}

#Biotic * PCA1############################################################################
try((fm10=occu(~ 1 ~ Biotic * PCA1, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10")){
  if(CondNum(fm10)<2000){
    if(CondNum(fm10)>0){mods=c(mods,fm10)}
}
}

#BioticYearly + PCA1##########################################################################
try((fm11=occu(~ 1 ~ BioticYearly + PCA1, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11")){
  if(CondNum(fm11)<2000){
    if(CondNum(fm11)>0){mods=c(mods,fm11)}
}
}

#BioticYearly * PCA1############################################################################
try((fm12=occu(~ 1 ~ BioticYearly * PCA1, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm12")){
  if(CondNum(fm12)<2000){
    if(CondNum(fm12)>0){mods=c(mods,fm12)}
}
}

#Biotic + PCA1 + Elevation##########################################################################
try((fm13=occu(~ 1 ~ Biotic + PCA1 + Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm13")){
  if(CondNum(fm13)<2000){
    if(CondNum(fm13)>0){mods=c(mods,fm13)}
}
}

#BioticYearly + PCA1 + Elevation############################################################################
try((fm14=occu(~ 1 ~ BioticYearly + PCA1 + Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm14")){
  if(CondNum(fm14)<2000){
    if(CondNum(fm14)>0){mods=c(mods,fm14)}
}
}

#Biotic * PCA1 + Elevation##########################################################################
try((fm15=occu(~ 1 ~ Biotic * PCA1 + Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm15")){
  if(CondNum(fm15)<2000){
    if(CondNum(fm15)>0){mods=c(mods,fm15)}
}
}

#Biotic + PCA1 * Elevation############################################################################
try((fm16=occu(~ 1 ~ Biotic + PCA1 * Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm16")){
  if(CondNum(fm16)<2000){
    if(CondNum(fm16)>0){mods=c(mods,fm16)}
}
}

#BioticAnnual * PCA1 + Elevation##########################################################################
try((fm17=occu(~ 1 ~ BioticAnnual * PCA1 + Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm17")){
  if(CondNum(fm17)<2000){
    if(CondNum(fm17)>0){mods=c(mods,fm17)}
}
}

#BioticAnnual + PCA1 * Elevation############################################################################
try((fm18=occu(~ 1 ~ BioticAnnual + PCA1 * Elevation, data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm18")){
  if(CondNum(fm18)<2000){
    if(CondNum(fm18)>0){mods=c(mods,fm18)}
}
}





######################################
#Model Selection
######################################

#if no models converged

if(isEmpty(mods)==TRUE){
  results.all[[k]]=NA  
  results.table.ma[[k]]=NA
  names(results.table.ma)[k] <- nms[k]
  results.table.aic[[k]]=NA
  names(results.table.aic)[k] <- nms[k]
  occu.transformed[[k]]=NA
  names(occu.transformed)[k] <- nms[k]
}else{
  models=fitList(fits=mods)
  (ms <- modSel(models))
  results.all[[k]]=ms    

models=fitList(fits=mods)

(ms <- modSel(models))

results.all[[k]]=ms
mods.all[[k]]=mods

#### add to Export to the end (including the bracket)

toExport<-as(ms,"data.frame") #add after ms object

#null.elev.aic=toExport$delta[toExport$formula=="~1 ~ Elevation"]
null.aic=toExport$delta[toExport$formula=="~1 ~ 1"]

#if null.elev didn't converge
#if(isEmpty(null.elev.aic)==TRUE){
#  null.elev=NA  
#}else{
#null.elev=toExport[toExport$formula=="~1 ~ Elevation",]		
#}

#if null didn't converge
if(isEmpty(null.aic)==TRUE){
	null=NA	
}else{
null=toExport[toExport$formula=="~1 ~ 1",]
}


#if((null.elev.aic==0 || isEmpty(null.elev.aic)==TRUE) || (null.aic==0) ||isEmpty(null.aic)==TRUE){	
	if(((null.aic==0) ||isEmpty(null.aic)==TRUE)){  

#results.table.ma[[k]]=rbind(null,null.elev)  
results.table.ma[[k]]=null	

temp=data.frame(toExport$formula,toExport$delta,toExport$AICwt)
names(temp)=c("formula","delta","AICwt")
#results.table.aic[[k]]=rbind(temp[temp$formula=="~1 ~ 1",],temp[temp$formula=="~ 1 ~ Elevation",])}else{
results.table.aic[[k]]=rbind(temp[temp$formula=="~1 ~ 1",])}else{
#results.table.ma[[k]]=rbind(toExport[1,],null,null.elev)
results.table.ma[[k]]=rbind(toExport[1,],null)
results.table.ma[[k]] <- cbind(nms[k], results.table.ma[[k]])

temp=data.frame(toExport$formula,toExport$delta,toExport$AICwt)
names(temp)=c("formula","delta","AICwt")
#results.table.aic[[k]]=rbind(temp[1,],temp[temp$formula=="~1 ~ 1",],temp[temp$formula=="~ 1 ~ Elevation",])
results.table.aic[[k]]=rbind(temp[1,],temp[temp$formula=="~1 ~ 1",])
results.table.aic[[k]] <- cbind(nms[k], results.table.aic[[k]])
}

test=seq(3,length(toExport)-2,by=2)
tmp=as.numeric(toExport[1,test])

occu.transformed[[k]]=exp(tmp)
occu.transformed[[k]] <- cbind(nms[k], occu.transformed[[k]])
}

#add tmp, temp and toExport to the rm object 
rm(fm0,fm0.1,fm1,fm1.1,fm1.2,fm2,fm2.1,fm2.2,fm3,fm3.1,fm3.2,
   fm4,fm4.1,fm4.2,fm5,fm5.1,fm5.2,fm6,fm6.1,fm6.2,fm7,fm7.1,
   fm7.2,fm7.3,fm7.4,fm7.5,fm8,fm8.1,fm8.2,fm9,fm9.1,fm9.2,
   fm10,fm10.1,fm10.2,fm11,fm11.1,fm11.2,fm12,fm12.1,fm12.2,
   fm13,fm13.1,fm13.2,fm14,fm14.1,fm14.2,fm15,fm15.1,fm15.2,
   fm16,fm16.1,fm16.2,fm17,fm17.1,fm17.2,fm18,fm18.1,fm18.2,
   fm19,fm19.1,fm19.2,fm20,fm20.1,fm20.2,fm21,fm21.1,fm21.2,
   fm22,fm22.1,fm22.2,fm23,fm23.1,fm23.2,fm24,fm24.1,fm24.2,
   fm25,fm25.1,fm25.2,fm26,fm26.1,fm26.2,fm27,fm27.1,fm27.2,
   fm28,fm28.1,fm28.2,fm29,fm29.1,fm29.2,fm30,fm30.1,fm30.2,
   fm31,fm31.1,fm31.2,fm32,fm32.1,fm32.2,fm33,fm33.1,fm33.2,
   fm34,fm34.1,fm34.2,fm35,fm35.1,fm35.2,
   fm36,fm36.1,fm36.2,
   fm37,fm37.1,fm37.2,
   fm38,fm38.1,fm38.2,
   fm39,fm39.1,fm39.2,
   fm40,fm40.1,fm40.2,
   fm41,fm41.1,fm41.2,
   fm42,fm42.1,fm42.2,
   mods,ms, tmp, temp, toExport, test)
}

# Need to coerce the lists into dataframes before writing to a files
results.table.ma.df <- ldply(results.table.ma, data.frame)
results.table.aic.df <- ldply(results.table.aic, data.frame)
occu.transformed.df <- ldply(occu.transformed, data.frame)

write.csv(results.table.ma.df, file="results.table.ma.csv")
write.csv(results.table.aic.df, file="results.table.aic.csv")
write.csv(occu.transformed.df, file="occu.transformed.csv")

length(nms)
for(i in 1:length(nms)) {
  outputname <- paste(nms[i], "occuAIC", "csv", sep=".")
  if(is.na(results.all[[i]])==TRUE){
    output <- NA
  }else{
    output <- results.all[[i]]@Full
  } 
  write.csv(output, file=outputname)
} 





