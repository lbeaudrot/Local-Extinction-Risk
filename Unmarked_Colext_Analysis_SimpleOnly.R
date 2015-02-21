# Unmarked analysis of 32 TEAM populations with >8% detection rates at sites with >500 m elevation change
library(unmarked)
library(plyr)
library(AICcmodavg)

#load('/Volumes/SCIENCEWORK/Working_folder/UPR_Prof/Collaborations/TEAM/32spp/All_covs_scaled.RData')
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
colext.transformed=list() #add at the beginning

isEmpty <- function(x) {
    return(length(x)==0)
}


CondNum <- function(model){
  max(eigen(hessian(model))$values)/min(eigen(hessian(model))$values)
}
# produces identical condition number to function "extractCN" from the AICcmodavg package

####
#for(k in 57:62){
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

#yearlySiteCovs
Tmin <- as.data.frame(sapply(covs, "[", 4))
Tmax <- as.data.frame(sapply(covs, "[", 5))
Tvar <- as.data.frame(sapply(covs, "[", 6))
Tsd <- as.data.frame(sapply(covs, "[", 7))
Tmean <- as.data.frame(sapply(covs, "[", 8))
BioticYearly <- BIOTIC_ALL_YEARS[[index]]
PCA1 <- data.frame(pca_covs)


to=dim(Tmin)[2]
years=as.character(1:to)
years=matrix(years,nrow(species),to,byrow=TRUE)

# ADD BIOTIC TO UMF COVARIATES
site.covs<-data.frame(Elevation, ForestLossCT, ForestGainCT, Biotic)

umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(Tmin=Tmin,Tmax=Tmax,Tvar=Tvar,Tsd=Tsd,Tmean=Tmean, BioticYearly=BioticYearly, PCA1=PCA1), siteCovs=site.covs, numPrimary=dim(Tmin)[2])
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

# ELEVATION only ##################################################################################
#try((fm2=colext(psiformula=~1,
#                gammaformula=~Elevation,
#                epsilonformula=~Elevation,
#                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm2")) {
#  if(CondNum(fm2)<2000){
#    if(CondNum(fm2)>0){mods=c(mods,fm2)}
#} 
#}

#try((fm2.1=colext(psiformula=~1,
#                  gammaformula=~Elevation,
#                  epsilonformula=~1,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm2.1")) {
#  if(CondNum(fm2.1)<2000){
#    if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
#} 
#}

#try((fm2.2=colext(psiformula=~1,
#                  gammaformula=~1,
#                  epsilonformula=~Elevation,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm2.2")) {
#  if(CondNum(fm2.2)<2000){
#    if(CondNum(fm2.2)>0){mods=c(mods,fm2.2)}
#} 
#}


# Temperature only ##################################################################################
try((fm3=colext(psiformula=~1,
                gammaformula=~Tmean + Tsd,
                epsilonformula=~Tmean + Tsd,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3")) {
  if(CondNum(fm3)<2000){
    if(CondNum(fm3)>0){mods=c(mods,fm3)}
} 
}

try((fm3.1=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.1")) {
  if(CondNum(fm3.1)<2000){
    if(CondNum(fm3.1)>0){mods=c(mods,fm3.1)}
} 
}

try((fm3.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmean + Tsd,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.2")) {
  if(CondNum(fm3.2)<2000){
    if(CondNum(fm3.2)>0){mods=c(mods,fm3.2)}
} 
}

# Biotic only ##################################################################################
try((fm4=colext(psiformula=~1,
                gammaformula=~Biotic,
                epsilonformula=~Biotic,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
if(exists("fm4")) {
  if(CondNum(fm4)<2000){
    if(CondNum(fm4)>0){mods=c(mods,fm4)}
} 
}   
   
try((fm4.1=colext(psiformula=~1,
                  gammaformula=~Biotic,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
if(exists("fm4.1")) {
  if(CondNum(fm4.1)<2000){
    if(CondNum(fm4.1)>0){mods=c(mods,fm4.1)}
} 
}   
   
try((fm4.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm4.2")) {
  if(CondNum(fm4.2)<2000){
    if(CondNum(fm4.2)>0){mods=c(mods,fm4.2)}
} 
}

# BioticYearly only ##################################################################################
try((fm5=colext(psiformula=~1,
                gammaformula=~BioticYearly,
                epsilonformula=~BioticYearly,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5")){
  if(CondNum(fm5)<2000){
    if(CondNum(fm5)>0){mods=c(mods,fm5)}
} 
}

try((fm5.1=colext(psiformula=~1,
                  gammaformula=~BioticYearly,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5.1")) {
  if(CondNum(fm5.1)<2000){
    if(CondNum(fm5.1)>0){mods=c(mods,fm5.1)}
} 
}

try((fm5.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~BioticYearly,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5.2")) {
  if(CondNum(fm5.2)<2000){
    if(CondNum(fm5.2)>0){mods=c(mods,fm5.2)}
} 
}

# PCA1 + Elevation ##################################################################################
#try((fm6=colext(psiformula=~1,
#                  gammaformula=~PCA1 + Elevation,
#                  epsilonformula=~PCA1 + Elevation,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm6")) {
#  if(CondNum(fm6)<2000){
#    if(CondNum(fm6)>0){mods=c(mods,fm6)}
#} 
#}

#try((fm6.1=colext(psiformula=~1,
#                  gammaformula=~PCA1 + Elevation,
#                  epsilonformula=~1,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm6.1")) {
# if(CondNum(fm6.1)<2000){
#    if(CondNum(fm6.1)>0){mods=c(mods,fm6.1)} 
#}
#}


#try((fm6.2=colext(psiformula=~1,
#                  gammaformula=~1,
#                  epsilonformula=~PCA1 + Elevation,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm6.2")) {
#  if(CondNum(fm6.2)<2000){
#    if(CondNum(fm6.2)>0){mods=c(mods,fm6.2)}
#} 
#}

# PCA1 * Elevation #############################################################################

#try((fm7=colext(psiformula=~1,
#                  gammaformula=~PCA1 * Elevation,
#                  epsilonformula=~PCA1 * Elevation,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm7")) {
#  if(CondNum(fm7)<2000){
#    if(CondNum(fm7)>0){mods=c(mods,fm7)}
#} 
#}

#try((fm7.1=colext(psiformula=~1,
#                  gammaformula=~PCA1 * Elevation,
#                  epsilonformula=~1,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm7.1")) {
#  if(CondNum(fm7.1)<2000){
#    if(CondNum(fm7.1)>0){mods=c(mods,fm7.1)}
#} 
#}

#try((fm7.2=colext(psiformula=~1,
#                  gammaformula=~1,
#                  epsilonformula=~PCA1 * Elevation,
#                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm7.2")) {
#  if(CondNum(fm7.2)<2000){
#    if(CondNum(fm7.2)>0){mods=c(mods,fm7.2)}
#} 
#}


# Tmean + Tsd + BioticYearly ##################################################################################
try((fm8=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd + BioticYearly,
                  epsilonformula=~Tmean + Tsd + BioticYearly,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8")) {
  if(CondNum(fm8)<2000){
    if(CondNum(fm8)>0){mods=c(mods,fm8)}
} 
}

try((fm8.1=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd + BioticYearly,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8.1")) {
 if(CondNum(fm8.1)<2000){
    if(CondNum(fm8.1)>0){mods=c(mods,fm8.1)} 
}
}


try((fm8.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmean + Tsd + BioticYearly,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8.2")) {
  if(CondNum(fm8.2)<2000){
    if(CondNum(fm8.2)>0){mods=c(mods,fm8.2)}
} 
}

# Tmean + Tsd * BioticYearly #############################################################################

try((fm9=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd * BioticYearly,
                  epsilonformula=~Tmean + Tsd * BioticYearly,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9")) {
  if(CondNum(fm9)<2000){
    if(CondNum(fm9)>0){mods=c(mods,fm9)}
} 
}

try((fm9.1=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd * BioticYearly,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9.1")) {
  if(CondNum(fm9.1)<2000){
    if(CondNum(fm9.1)>0){mods=c(mods,fm9.1)}
} 
}

try((fm9.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmean + Tsd * BioticYearly,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9.2")) {
  if(CondNum(fm9.2)<2000){
    if(CondNum(fm9.2)>0){mods=c(mods,fm9.2)}
} 
}

# Tmean + Tsd + Biotic ##################################################################################
try((fm10=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd + Biotic,
                  epsilonformula=~Tmean + Tsd + Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10")) {
  if(CondNum(fm10)<2000){
    if(CondNum(fm10)>0){mods=c(mods,fm10)}
} 
}

try((fm10.1=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd + Biotic,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10.1")) {
 if(CondNum(fm10.1)<2000){
    if(CondNum(fm10.1)>0){mods=c(mods,fm10.1)} 
}
}


try((fm10.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmean + Tsd + Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10.2")) {
  if(CondNum(fm10.2)<2000){
    if(CondNum(fm10.2)>0){mods=c(mods,fm10.2)}
} 
}

# Tmean + Tsd * Biotic #############################################################################

try((fm11=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd * Biotic,
                  epsilonformula=~Tmean + Tsd * Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11")) {
  if(CondNum(fm11)<2000){
    if(CondNum(fm11)>0){mods=c(mods,fm11)}
} 
}

try((fm11.1=colext(psiformula=~1,
                  gammaformula=~Tmean + Tsd * Biotic,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11.1")) {
  if(CondNum(fm11.1)<2000){
    if(CondNum(fm11.1)>0){mods=c(mods,fm11.1)}
} 
}

try((fm11.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmean + Tsd * Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11.2")) {
  if(CondNum(fm11.2)<2000){
    if(CondNum(fm11.2)>0){mods=c(mods,fm11.2)}
} 
}



# Biotic + Elevation ##########################################################################
#try((fm28=colext(psiformula=~1,
#                 gammaformula=~Biotic + Elevation,
#                 epsilonformula=~Biotic + Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm28")) {
#  if(CondNum(fm28)<2000){
#    if(CondNum(fm28)>0){mods=c(mods,fm28)}
#} 
#}

#try((fm28.1=colext(psiformula=~1,
#                   gammaformula=~Biotic + Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm28.1")) {
#  if(CondNum(fm28.1)<2000){
#    if(CondNum(fm28.1)>0){mods=c(mods,fm28.1)}
#} 
#}

#try((fm28.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~Biotic + Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm28.2")) {
#  if(CondNum(fm28.2)<2000){
#    if(CondNum(fm28.2)>0){mods=c(mods,fm28.2)}
#} 
#}


#Biotic * Elevation##########################################################################
#try((fm30=colext(psiformula=~1,
#                 gammaformula=~Biotic * Elevation,
#                 epsilonformula=~Biotic * Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm30")) {
#  if(CondNum(fm30)<2000){
#    if(CondNum(fm30)>0){mods=c(mods,fm30)}
#} 
#}

#try((fm30.1=colext(psiformula=~1,
#                   gammaformula=~Biotic * Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm30.1")) {
#  if(CondNum(fm30.1)<2000){
#    if(CondNum(fm30.1)>0){mods=c(mods,fm30.1)}
#} 
#}

#try((fm30.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~Biotic * Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm30.2")) {
#  if(CondNum(fm30.2)<2000){
#    if(CondNum(fm30.2)>0){mods=c(mods,fm30.2)}
#} 
#}



# BioticYearly + Elevation ##########################################################################
#try((fm32=colext(psiformula=~1,
#                 gammaformula=~BioticYearly + Elevation,
#                 epsilonformula=~BioticYearly + Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm32")) {
#  if(CondNum(fm32)<2000){
#    if(CondNum(fm32)>0){mods=c(mods,fm32)}
#} 
#}

#try((fm32.1=colext(psiformula=~1,
#                   gammaformula=~BioticYearly + Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm32.1")) {
#  if(CondNum(fm32.1)<2000){
#    if(CondNum(fm32.1)>0){mods=c(mods,fm32.1)}
#} 
#}

#try((fm32.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~BioticYearly + Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm32.2")) {
#  if(CondNum(fm32.2)<2000){
#    if(CondNum(fm32.2)>0){mods=c(mods,fm32.2)}
#} 
#}



# BioticYearly * Elevation ##########################################################################
#try((fm35=colext(psiformula=~1,
#                 gammaformula=~BioticYearly * Elevation,
#                 epsilonformula=~BioticYearly * Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm35")) {
#  if(CondNum(fm35)<2000){
#    if(CondNum(fm35)>0){mods=c(mods,fm35)}
#} 
#}

#try((fm35.1=colext(psiformula=~1,
#                   gammaformula=~BioticYearly * Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm35.1")){
#  if(CondNum(fm35.1)<2000){
#    if(CondNum(fm35.1)>0){mods=c(mods,fm35.1)}
#} 
#}

#try((fm35.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~BioticYearly * Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm35.2")) {
#  if(CondNum(fm35.2)<2000){
#    if(CondNum(fm35.2)>0){mods=c(mods,fm35.2)}
#} 
#}


# PCA1 + Biotic + Elevation ##########################################################################
#try((fm36=colext(psiformula=~1,
#                 gammaformula=~PCA1 + Biotic + Elevation,
#                 epsilonformula=~PCA1 + Biotic + Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm36")) {
#  if(CondNum(fm36)<2000){
#    if(CondNum(fm36)>0){mods=c(mods,fm36)}
#} 
#}

#try((fm36.1=colext(psiformula=~1,
#                   gammaformula=~PCA1 + Biotic + Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm36.1")){
#  if(CondNum(fm36.1)<2000){
#    if(CondNum(fm36.1)>0){mods=c(mods,fm36.1)}
#} 
#}

#try((fm36.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~PCA1 + Biotic + Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm36.2")) {
#  if(CondNum(fm36.2)<2000){
#    if(CondNum(fm36.2)>0){mods=c(mods,fm36.2)}
#} 
#}


# PCA1 + BioticYearly + Elevation ##########################################################################
#try((fm38=colext(psiformula=~1,
#                 gammaformula=~PCA1 + BioticYearly + Elevation,
#                 epsilonformula=~PCA1 + BioticYearly + Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm38")) {
#  if(CondNum(fm38)<2000){
#    if(CondNum(fm38)>0){mods=c(mods,fm38)}
#} 
#}

#try((fm38.1=colext(psiformula=~1,
#                   gammaformula=~PCA1 + BioticYearly + Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm38.1")) {
#  if(CondNum(fm38.1)<2000){
#    if(CondNum(fm38.1)>0){mods=c(mods,fm38.1)}
#} 
#}

#try((fm38.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~PCA1 + BioticYearly + Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm38.2")) {
#  if(CondNum(fm38.2)<2000){
#    if(CondNum(fm38.2)>0){mods=c(mods,fm38.2)}
#} 
#}


# PCA1 * BioticYearly + Elevation ##########################################################################
#try((fm37=colext(psiformula=~1,
#                 gammaformula=~PCA1 * BioticYearly + Elevation,
#                 epsilonformula=~PCA1 * BioticYearly + Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm37")) {
#  if(CondNum(fm37)<2000){
#    if(CondNum(fm37)>0){mods=c(mods,fm37)}
#} 
#}

#try((fm37.1=colext(psiformula=~1,
#                   gammaformula=~PCA1 * BioticYearly + Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm37.1")) {
#  if(CondNum(fm37.1)<2000){
#    if(CondNum(fm37.1)>0){mods=c(mods,fm37.1)}
#} 
#}

#try((fm37.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~PCA1 * BioticYearly + Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm37.2")) {
#  if(CondNum(fm37.2)<2000){
#    if(CondNum(fm37.2)>0){mods=c(mods,fm37.2)}
#} 
#}



#PCA1 + BioticYearly * Elevation ##########################################################################
#try((fm39=colext(psiformula=~1,
#                 gammaformula=~PCA1 + BioticYearly * Elevation,
#                 epsilonformula=~PCA1 + BioticYearly * Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm39")) {
#  if(CondNum(fm39)<2000){
#    if(CondNum(fm39)>0){mods=c(mods,fm39)}
#} 
#}

#try((fm39.1=colext(psiformula=~1,
#                   gammaformula=~PCA1 + BioticYearly * Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm39.1")) {
#  if(CondNum(fm39.1)<2000){
#    if(CondNum(fm39.1)>0){mods=c(mods,fm39.1)}
#} 
#}

#try((fm39.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~PCA1 + BioticYearly * Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm39.2")) {
#  if(CondNum(fm39.2)<2000){
#    if(CondNum(fm39.2)>0){mods=c(mods,fm39.2)}
#} 
#}




# PCA1 * Biotic + Elevation ##########################################################################
#try((fm40=colext(psiformula=~1,
#                 gammaformula=~PCA1 * Biotic + Elevation,
#                 epsilonformula=~PCA1 * Biotic + Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm40")) {
#  if(CondNum(fm40)<2000){
#    if(CondNum(fm40)>0){mods=c(mods,fm40)}
#} 
#}

#try((fm40.1=colext(psiformula=~1,
#                   gammaformula=~PCA1 * Biotic + Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm40.1")) {
#  if(CondNum(fm40.1)<2000){
#    if(CondNum(fm40.1)>0){mods=c(mods,fm40.1)}
#} 
#}

#try((fm40.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~PCA1 * Biotic + Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm40.2")) {
#  if(CondNum(fm40.2)<2000){
#    if(CondNum(fm40.2)>0){mods=c(mods,fm40.2)}
#} 
#}




#PCA1 + Biotic * Elevation ##########################################################################
#try((fm42=colext(psiformula=~1,
#                 gammaformula=~PCA1 + Biotic * Elevation,
#                 epsilonformula=~PCA1 + Biotic * Elevation,
#                 pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm42")) {
#  if(CondNum(fm42)<2000){
#    if(CondNum(fm42)>0){mods=c(mods,fm42)}
#} 
#}

#try((fm42.1=colext(psiformula=~1,
#                   gammaformula=~PCA1 + Biotic * Elevation,
#                   epsilonformula=~1,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm42.1")) {
#  if(CondNum(fm42.1)<2000){
#    if(CondNum(fm42.1)>0){mods=c(mods,fm42.1)}
#} 
#}

#try((fm42.2=colext(psiformula=~1,
#                   gammaformula=~1,
#                   epsilonformula=~PCA1 + Biotic * Elevation,
#                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm42.2")) {
#  if(CondNum(fm42.2)<2000){
#    if(CondNum(fm42.2)>0){mods=c(mods,fm42.2)}
#} 
#}




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
   fm4,fm4.1,fm4.2,fm5,fm5.1,fm5.2,fm6, fm6.1,fm6.2,fm7,fm7.1,
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
