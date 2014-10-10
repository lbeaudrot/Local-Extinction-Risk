# Unmarked analysis of 32 TEAM populations with >8% detection rates at sites with >500 m elevation change
#load('/Volumes/SCIENCEWORK/Working_folder/UPR_Prof/Collaborations/TEAM/32spp/All_covs_scaled.RData')
#load('/Volumes/SCIENCEWORK/Working_folder/UPR_Prof/Collaborations/TEAM/32spp/All500m_covariate_species.RData')
library(unmarked)

# Matrices for each population are contained in the object "All500m_covariate_species"
nms=names(All500m_covariate_species)

results.all=list()
mods.all=list()

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
ForestLossPA <- unlist(as.matrix(sapply(covs, "[", 3)))

#yearlySiteCovs
Tmin <- as.data.frame(sapply(covs, "[", 4))
Tmax <- as.data.frame(sapply(covs, "[", 5))
Tvar <- as.data.frame(sapply(covs, "[", 6))

to=dim(Tmin)[2]
years=as.character(1:to)
years=matrix(years,nrow(species),to,byrow=TRUE)


site.covs<-data.frame(Elevation, ForestLossCT, ForestLossPA)

umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(year=years,Tmin=Tmin,Tmax=Tmax,Tvar=Tvar), siteCovs=site.covs, numPrimary=dim(Tmin)[2])

mods=list()

#Null
try((fm0=colext(psiformula=~1,gammaformula=~1,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm0")){mods=c(mods,fm0)}

#NULL + ELEVATION for initial occupancy
try((fm0.1=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm0.1")){mods=c(mods,fm0.1)}

#Elevation
try((fm1=colext(psiformula=~Elevation,gammaformula=~Elevation,epsilonformula=~Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm1")){mods=c(mods,fm1)}

try((fm1.1=colext(psiformula=~Elevation,gammaformula=~Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm1.1")){mods=c(mods,fm1.1)}



try((fm1.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm1.2")){mods=c(mods,fm1.2)}



#Tmin
try((fm2=colext(psiformula=~Elevation,gammaformula=~Tmin,epsilonformula=~Tmin,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2")){mods=c(mods,fm2)}


try((fm2.1=colext(psiformula=~Elevation,gammaformula=~Tmin,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.1")){mods=c(mods,fm2.1)}


try((fm2.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.2")){mods=c(mods,fm2.2)}



#Tmax
try((fm3=colext(psiformula=~Elevation,gammaformula=~Tmax,epsilonformula=~Tmax,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3")){mods=c(mods,fm3)}



try((fm3.1=colext(psiformula=~Elevation,gammaformula=~Tmax,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.1")){mods=c(mods,fm3.1)}



try((fm3.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.2")){mods=c(mods,fm3.2)}


#Tvar
   try((fm4=colext(psiformula=~Elevation,gammaformula=~Tvar,epsilonformula=~Tvar,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
   if(exists("fm4")){mods=c(mods,fm4)}
   
   
   
   try((fm4.1=colext(psiformula=~Elevation,gammaformula=~Tvar,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
   if(exists("fm4.1")){mods=c(mods,fm4.1)}
   
   
   
try((fm4.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm4.2")){mods=c(mods,fm4.2)}


#ForestLossCT
try((fm5=colext(psiformula=~Elevation,gammaformula=~ForestLossCT,epsilonformula=~ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5")){mods=c(mods,fm5)}



try((fm5.1=colext(psiformula=~Elevation,gammaformula=~ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5.1")){mods=c(mods,fm5.1)}



try((fm5.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm5.2")){mods=c(mods,fm5.2)}

#ForestLossPA
try((fm6=colext(psiformula=~Elevation,gammaformula=~ForestLossPA,epsilonformula=~ForestLossPA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6")){mods=c(mods,fm6)}



try((fm6.1=colext(psiformula=~Elevation,gammaformula=~ForestLossPA,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6.1")){mods=c(mods,fm6)}



try((fm6.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~ForestLossPA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6.2")){mods=c(mods,fm6.2)}


#Time-dependent
try((fm7=colext(psiformula=~1,gammaformula=~year-1,epsilonformula=~year-1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7")){mods=c(mods,fm7)}



try((fm7.1=colext(psiformula=~1,gammaformula=~year-1,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.1")){mods=c(mods,fm7.1)}



try((fm7.2=colext(psiformula=~1,gammaformula=~1,epsilonformula=~year-1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.2")){mods=c(mods,fm7.2)}




#Time+Tmin
try((fm8=colext(psiformula=~1,gammaformula=~(year-1)+Tmin,epsilonformula=~(year-1)+Tmin,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8")){mods=c(mods,fm8)}



try((fm8.1=colext(psiformula=~1,gammaformula=~(year-1)+Tmin,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8")){mods=c(mods,fm8)}



try((fm8.2=colext(psiformula=~1,gammaformula=~1,epsilonformula=~(year-1)+Tmin,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm8.2")){mods=c(mods,fm8.2)}


#Time+tmax
try((fm9=colext(psiformula=~1,gammaformula=~(year-1)+Tmax,epsilonformula=~(year-1)+Tmax,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9")){mods=c(mods,fm9)}



try((fm9.1=colext(psiformula=~1,gammaformula=~(year-1)+Tmax,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9.1")){mods=c(mods,fm9.1)}



try((fm9.2=colext(psiformula=~1,gammaformula=~1,epsilonformula=~(year-1)+Tmax,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm9.2")){mods=c(mods,fm9.2)}


#Time+tvar
try(
(fm10=colext(psiformula=~1,gammaformula=~(year-1)+Tvar,epsilonformula=~(year-1)+Tvar,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10")){mods=c(mods,fm10)}


try((fm10.1=colext(psiformula=~1,gammaformula=~(year-1)+Tvar,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10.1")){mods=c(mods,fm10.1)}



try((fm10.2=colext(psiformula=~1,gammaformula=~1,epsilonformula=~(year-1)+Tvar,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm10.2")){mods=c(mods,fm10.2)}


#Time+Elevation
try((fm11=colext(psiformula=~1,gammaformula=~(year-1)+Elevation,epsilonformula=~(year-1)+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11")){mods=c(mods,fm11)}



try((fm11.1=colext(psiformula=~1,gammaformula=~(year-1)+Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11.1")){mods=c(mods,fm11.1)}



try((fm11.2=colext(psiformula=~1,gammaformula=~1,epsilonformula=~(year-1)+Elevation,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm11.2")){mods=c(mods,fm11.2)}


#Time+ForestLossCT
try((fm12=colext(psiformula=~1,gammaformula=~(year-1)+ForestLossCT,epsilonformula=~(year-1)+ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm12")){mods=c(mods,fm12)}



try((fm12.1=colext(psiformula=~1,gammaformula=~(year-1)+ForestLossCT,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm12.1")){mods=c(mods,fm12.1)}



try((fm12.2=colext(psiformula=~1,gammaformula=~1,epsilonformula=~(year-1)+ForestLossCT,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm12.2")){mods=c(mods,fm12.2)}

#Time+ForestLossPA
try((fm13=colext(psiformula=~1,gammaformula=~(year-1)+ForestLossPA,epsilonformula=~(year-1)+ForestLossPA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm13")){mods=c(mods,fm13)}

try((fm13.1=colext(psiformula=~1,gammaformula=~(year-1)+ForestLossPA,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm13.1")){mods=c(mods,fm13.1)}

try((fm13.2=colext(psiformula=~1,gammaformula=~1,epsilonformula=~(year-1)+ForestLossPA,pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm13.2")){mods=c(mods,fm13.2)}

models=fitList(fits=mods)

(ms <- modSel(models))

results.all[[k]]=ms
mods.all[[k]]=mods
rm(fm0,fm0.1,fm1,fm1.1,fm1.2,fm2,fm2.1,fm2.2,fm3,fm3.1,fm3.2,fm4,fm4.1,fm4.2,fm5,fm5.1,fm5.2,fm6,fm6.1,fm6.2,fm7,fm7.1,fm7.2,fm8,fm8.1,fm8.2,fm9,fm9.1,fm9.2,fm10,fm10.1,fm10.2,fm11,fm11.1,fm11.2,fm12,fm12.1,fm12.2,fm13,fm13.1,fm13.2,mods,ms)
}

save.image(file="spp32_results.RData")
 