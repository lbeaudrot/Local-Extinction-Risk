# Unmarked analysis of 32 TEAM populations with >8% detection rates at sites with >500 m elevation change
#load('/Volumes/SCIENCEWORK/Working_folder/UPR_Prof/Collaborations/TEAM/32spp/All_covs_scaled.RData')
#load('/Volumes/SCIENCEWORK/Working_folder/UPR_Prof/Collaborations/TEAM/32spp/All500m_covariate_species.RData')
rm(list=ls())
load("All_covs.RData")
#load("All500m_covariate_species.RData")
load("All_species7sites.RData")
All500m_covariate_species <- All_species7sites

library(unmarked)

# Matrices for each population are contained in the object "All500m_covariate_species"
nms=names(All500m_covariate_species)

results.all=list()
mods.all=list()

for(k in 27:30){
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

#yearlySiteCovs
Tmin <- as.data.frame(sapply(covs, "[", 4))
Tmax <- as.data.frame(sapply(covs, "[", 5))
Tvar <- as.data.frame(sapply(covs, "[", 6))
Tsd <- as.data.frame(sapply(covs, "[", 7))
Tmean <- as.data.frame(sapply(covs, "[", 8))



to=dim(Tmin)[2]
years=as.character(1:to)
years=matrix(years,nrow(species),to,byrow=TRUE)


site.covs<-data.frame(Elevation, ForestLossCT, ForestGainCT)

umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(Tmin=Tmin,Tmax=Tmax,Tvar=Tvar,Tsd=Tsd,Tmean=Tmean), siteCovs=site.covs, numPrimary=dim(Tmin)[2])
#umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(year=years,Tmin=Tmin,Tmax=Tmax,Tvar=Tvar,Tsd=Tsd,Tmean=Tmean), siteCovs=site.covs, numPrimary=dim(Tmin)[2])

mods=list()

#Null##################################################################################
try((fm0=colext(psiformula=~1,gammaformula=~1,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm0")){mods=c(mods,fm0)}

#NULL + ELEVATION for initial occupancy
try((fm0.1=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm0.1")){mods=c(mods,fm0.1)}


#Tmin##################################################################################
try((fm2=colext(psiformula=~Elevation,gammaformula=~Tmin,epsilonformula=~Tmin,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2")){mods=c(mods,fm2)}


try((fm2.1=colext(psiformula=~Elevation,gammaformula=~Tmin,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.1")){mods=c(mods,fm2.1)}


try((fm2.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.2")){mods=c(mods,fm2.2)}



#Tmax##################################################################################
try((fm3=colext(psiformula=~Elevation,gammaformula=~Tmax,epsilonformula=~Tmax,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3")){mods=c(mods,fm3)}



try((fm3.1=colext(psiformula=~Elevation,gammaformula=~Tmax,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.1")){mods=c(mods,fm3.1)}



try((fm3.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.2")){mods=c(mods,fm3.2)}


#Tvar##################################################################################
   try((fm4=colext(psiformula=~Elevation,gammaformula=~Tvar,epsilonformula=~Tvar,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
   if(exists("fm4")){mods=c(mods,fm4)}
   
   
   
   try((fm4.1=colext(psiformula=~Elevation,gammaformula=~Tvar,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
   if(exists("fm4.1")){mods=c(mods,fm4.1)}
   
   
   
try((fm4.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm4.2")){mods=c(mods,fm4.2)}


#ForestLossCT##################################################################################
try((fm5=colext(psiformula=~Elevation,gammaformula=~ForestLossCT,epsilonformula=~ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5")){mods=c(mods,fm5)}



try((fm5.1=colext(psiformula=~Elevation,gammaformula=~ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5.1")){mods=c(mods,fm5.1)}



try((fm5.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5.2")){mods=c(mods,fm5.2)}

#ForestGainCT##################################################################################
try((fm6=colext(psiformula=~Elevation,gammaformula=~ForestGainCT,epsilonformula=~ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6")){mods=c(mods,fm6)}



try((fm6.1=colext(psiformula=~Elevation,gammaformula=~ForestGainCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6.1")){mods=c(mods,fm6)}



try((fm6.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6.2")){mods=c(mods,fm6.2)}

#Tmean*Tsd##################################################################################
try((fm7=colext(psiformula=~Elevation,gammaformula=~Tmean*Tsd,epsilonformula=~Tmean*Tsd,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7")){mods=c(mods,fm7)}

try((fm7.1=colext(psiformula=~Elevation,gammaformula=~Tmean*Tsd,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.1")){mods=c(mods,fm7.1)}

try((fm7.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmean*Tsd,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.2")){mods=c(mods,fm7.2)}

try((fm7.3=colext(psiformula=~Elevation,gammaformula=~Tmean+Tsd,epsilonformula=~Tmean+Tsd,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.3")){mods=c(mods,fm7.3)}

try((fm7.4=colext(psiformula=~Elevation,gammaformula=~Tmean+Tsd,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.4")){mods=c(mods,fm7.4)}

try((fm7.5=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmean+Tsd,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.5")){mods=c(mods,fm7.5)}

#Tmin+ForestLossCT############################################################################
try((fm8=colext(psiformula=~Elevation,gammaformula=~Tmin+ForestLossCT,epsilonformula=~Tmin+ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8")){mods=c(mods,fm8)}

try((fm8.1=colext(psiformula=~Elevation,gammaformula=~Tmin+ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8.1")){mods=c(mods,fm8.1)}

try((fm8.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin+ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8.2")){mods=c(mods,fm8.2)}

#Tmin * ForestLossCT##########################################################################
try((fm9=colext(psiformula=~Elevation,gammaformula=~Tmin*ForestLossCT,epsilonformula=~Tmin*ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9")){mods=c(mods,fm9)}

try((fm9.1=colext(psiformula=~Elevation,gammaformula=~Tmin*ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9.1")){mods=c(mods,fm9.1)}

try((fm9.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin*ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9.2")){mods=c(mods,fm9.2)}


#Tmin+ForestGainCT############################################################################
try((fm10=colext(psiformula=~Elevation,gammaformula=~Tmin+ForestGainCT,epsilonformula=~Tmin+ForestLossCA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10")){mods=c(mods,fm10)}

try((fm10.1=colext(psiformula=~Elevation,gammaformula=~Tmin+ForestGainCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10.1")){mods=c(mods,fm10.1)}

try((fm10.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin+ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10.2")){mods=c(mods,fm10.2)}

#Tmin * ForestGainCT##########################################################################
try((fm11=colext(psiformula=~Elevation,gammaformula=~Tmin*ForestGainCT,epsilonformula=~Tmin*ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11")){mods=c(mods,fm11)}

try((fm11.1=colext(psiformula=~Elevation,gammaformula=~Tmin*ForestGainCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11.1")){mods=c(mods,fm11.1)}

try((fm11.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin*ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11.2")){mods=c(mods,fm11.2)}


#Tvar+ForestGainCT############################################################################
try((fm12=colext(psiformula=~Elevation,gammaformula=~Tvar+ForestGainCT,epsilonformula=~Tvar+ForestLossCA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm12")){mods=c(mods,fm12)}

try((fm12.1=colext(psiformula=~Elevation,gammaformula=~Tvar+ForestGainCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm12.1")){mods=c(mods,fm12.1)}

try((fm12.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar+ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm12.2")){mods=c(mods,fm12.2)}

#Tvar * ForestGainCT##########################################################################
try((fm13=colext(psiformula=~Elevation,gammaformula=~Tvar*ForestGainCT,epsilonformula=~Tvar*ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm13")){mods=c(mods,fm13)}

try((fm13.1=colext(psiformula=~Elevation,gammaformula=~Tvar*ForestGainCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm13.1")){mods=c(mods,fm13.1)}

try((fm13.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar*ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm13.2")){mods=c(mods,fm13.2)}

#Tvar+ForestLossCT############################################################################
try((fm14=colext(psiformula=~Elevation,gammaformula=~Tvar+ForestLossCT,epsilonformula=~Tvar+ForestLossCA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm14")){mods=c(mods,fm14)}

try((fm14.1=colext(psiformula=~Elevation,gammaformula=~Tvar+ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm14.1")){mods=c(mods,fm14.1)}

try((fm14.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar+ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm14.2")){mods=c(mods,fm14.2)}

#Tvar * ForestLossCT##########################################################################
try((fm15=colext(psiformula=~Elevation,gammaformula=~Tvar*ForestLossCT,epsilonformula=~Tvar*ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm15")){mods=c(mods,fm15)}

try((fm15.1=colext(psiformula=~Elevation,gammaformula=~Tvar*ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm15.1")){mods=c(mods,fm15.1)}

try((fm15.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar*ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm15.2")){mods=c(mods,fm15.2)}

#Tmax+ForestGainCT############################################################################
try((fm16=colext(psiformula=~Elevation,gammaformula=~Tmax+ForestGainCT,epsilonformula=~Tmax+ForestLossCA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm16")){mods=c(mods,fm16)}

try((fm16.1=colext(psiformula=~Elevation,gammaformula=~Tmax+ForestGainCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm16.1")){mods=c(mods,fm16.1)}

try((fm16.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax+ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm16.2")){mods=c(mods,fm16.2)}

#Tmax * ForestGainCT##########################################################################
try((fm17=colext(psiformula=~Elevation,gammaformula=~Tmax*ForestGainCT,epsilonformula=~Tmax*ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm17")){mods=c(mods,fm17)}

try((fm17.1=colext(psiformula=~Elevation,gammaformula=~Tmax*ForestGainCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm17.1")){mods=c(mods,fm17.1)}

try((fm17.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax*ForestGainCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm17.2")){mods=c(mods,fm17.2)}

#Tmax+ForestLossCT############################################################################
try((fm18=colext(psiformula=~Elevation,gammaformula=~Tmax+ForestLossCT,epsilonformula=~Tmax+ForestLossCA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm18")){mods=c(mods,fm18)}

try((fm18.1=colext(psiformula=~Elevation,gammaformula=~Tmax+ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm18.1")){mods=c(mods,fm18.1)}

try((fm18.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax+ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm18.2")){mods=c(mods,fm18.2)}

#Tmax * ForestLossCT##########################################################################
try((fm19=colext(psiformula=~Elevation,gammaformula=~Tmax*ForestLossCT,epsilonformula=~Tmax*ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm19")){mods=c(mods,fm19)}

try((fm19.1=colext(psiformula=~Elevation,gammaformula=~Tmax*ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm19.1")){mods=c(mods,fm19.1)}

try((fm19.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax*ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm19.2")){mods=c(mods,fm19.2)}

### NEW FROM LB ####

#Tmax * Elevation##########################################################################
try((fm20=colext(psiformula=~Elevation,gammaformula=~Tmax*Elevation,epsilonformula=~Tmax*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm20")){mods=c(mods,fm20)}

try((fm20.1=colext(psiformula=~Elevation,gammaformula=~Tmax*Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm20.1")){mods=c(mods,fm20.1)}

try((fm20.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm20.2")){mods=c(mods,fm20.2)}

#Tmax + Elevation##########################################################################
try((fm24=colext(psiformula=~Elevation,gammaformula=~Tmax+Elevation,epsilonformula=~Tmax+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm24")){mods=c(mods,fm24)}

try((fm24.1=colext(psiformula=~Elevation,gammaformula=~Tmax+Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm24.1")){mods=c(mods,fm24.1)}

try((fm24.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm24.2")){mods=c(mods,fm24.2)}


#Tmin * Elevation##########################################################################
try((fm21=colext(psiformula=~Elevation,gammaformula=~Tmin*Elevation,epsilonformula=~Tmin*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm21")){mods=c(mods,fm21)}

try((fm21.1=colext(psiformula=~Elevation,gammaformula=~Tmin*Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm21.1")){mods=c(mods,fm21.1)}

try((fm21.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm21.2")){mods=c(mods,fm21.2)}

#Tmin + Elevation##########################################################################
try((fm25=colext(psiformula=~Elevation,gammaformula=~Tmin+Elevation,epsilonformula=~Tmin+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm25")){mods=c(mods,fm25)}

try((fm25.1=colext(psiformula=~Elevation,gammaformula=~Tmin+Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm25.1")){mods=c(mods,fm25.1)}

try((fm25.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm25.2")){mods=c(mods,fm25.2)}



#Tvar * Elevation##########################################################################
try((fm22=colext(psiformula=~Elevation,gammaformula=~Tvar*Elevation,epsilonformula=~Tvar*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm22")){mods=c(mods,fm22)}

try((fm22.1=colext(psiformula=~Elevation,gammaformula=~Tvar*Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm22.1")){mods=c(mods,fm22.1)}

try((fm22.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm22.2")){mods=c(mods,fm22.2)}


#Tvar + Elevation##########################################################################
try((fm26=colext(psiformula=~Elevation,gammaformula=~Tvar+Elevation,epsilonformula=~Tvar+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm26")){mods=c(mods,fm26)}

try((fm26.1=colext(psiformula=~Elevation,gammaformula=~Tvar+Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm26.1")){mods=c(mods,fm26.1)}

try((fm26.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm26.2")){mods=c(mods,fm26.2)}



#ForestLossCT * Elevation##########################################################################
try((fm23=colext(psiformula=~Elevation,gammaformula=~ForestLossCT*Elevation,epsilonformula=~ForestLossCT*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm23")){mods=c(mods,fm23)}

try((fm23.1=colext(psiformula=~Elevation,gammaformula=~ForestLossCT*Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm23.1")){mods=c(mods,fm23.1)}

try((fm23.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~ForestLossCT*Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm23.2")){mods=c(mods,fm23.2)}

#ForestLossCT + Elevation##########################################################################
try((fm27=colext(psiformula=~Elevation,gammaformula=~ForestLossCT+Elevation,epsilonformula=~ForestLossCT+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm27")){mods=c(mods,fm27)}

try((fm27.1=colext(psiformula=~Elevation,gammaformula=~ForestLossCT+Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm27.1")){mods=c(mods,fm27.1)}

try((fm27.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~ForestLossCT+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm27.2")){mods=c(mods,fm27.2)}




######################################
#Model Selection
######################################

models=fitList(fits=mods)

(ms <- modSel(models))

results.all[[k]]=ms
mods.all[[k]]=mods
rm(fm0,fm0.1,fm1,fm1.1,fm1.2,fm2,fm2.1,fm2.2,fm3,fm3.1,fm3.2,fm4,fm4.1,fm4.2,fm5,fm5.1,fm5.2,fm7,fm7.1,fm7.2,fm7.3,fm7.4,fm7.5,fm8,fm8.1,fm8.2,fm9,fm9.1,fm9.2,fm10,fm10.1,fm10.2,fm11,fm11.1,fm11.2,fm12,fm12.1,fm12.2,fm13,fm13.1,fm13.2,fm14,fm14.1,fm14.2,fm15,fm15.1,fm15.2,fm16,fm16.1,fm16.2,fm17,fm17.1,fm17.2,fm18,fm18.1,fm18.2,fm19,fm19.1,fm19.2,fm20,fm20.1,fm20.2,fm21,fm21.1,fm21.2,fm22,fm22.1,fm22.2,fm23,fm23.1,fm23.2,fm24,fm24.1,fm24.2,fm25,fm25.1,fm25.2,fm26,fm26.1,fm26.2,fm27,fm27.1,fm27.2,mods,ms)
}


save.image(file="sppAll_results_v2.RData")
names(results.all) <- nms

length(nms)
for(i in 1:length(nms)) {
  outputname <- paste(nms[i], "colextAIC", "csv", sep=".")
  output <- results.all[[i]]@Full
	write.csv(output, file=outputname)
} 