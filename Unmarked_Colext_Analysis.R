# Unmarked analysis of 32 TEAM populations with >8% detection rates at sites with >500 m elevation change

library(unmarked)

# Matrices for each population are contained in the object "All500m_covariate_species"
names(All500m_covariate_species)

# DEFINE SPECIES for analysis and site USING INDEX VALUE for list of all species (see previous call for list of species names)
index <- 1
sp.name <- names(All500m_covariate_species)[index]
sp.name
species <- All500m_covariate_species[[index]]
site <- substr(names(All500m_covariate_species)[[index]],1,3)
site

# Define covariate object based on site
site_covs <- paste(site, "covs", sep="_")
covs <- All_covs[names(All_covs)==site_covs]
Elevation <- data.frame(sapply(covs, "[", 1))
ForestLossCT <- data.frame(sapply(covs, "[", 2))
ForestLossPA <- data.frame(sapply(covs, "[", 3))
Tmin <- data.frame(sapply(covs, "[", 4))
Tmax <- data.frame(sapply(covs, "[", 5))
Tvar <- data.frame(sapply(covs, "[", 6))
d


###############################
#Scale and centering covariates
###############################
#Elevation=scale(Elevation)

#ForestLossCT=scale(ForestLossCT)

# NB ForestLossPA not scaled because each TEAM site has a single value

#Tmin=scale(Tmin)
#Tmin[is.na(Tmin)]<-0

#Tmax=scale(Tmax)
#Tmax[is.na(Tmax)]<-0

#Tvar=scale(Tvar)
#Tvar[is.na(Tvar)]<-0

site.covs<-data.frame(Elevation, ForestLossCT, ForestLossPA)

umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(Tmin=Tmin,Tmax=Tmax,Tvar=Tvar), siteCovs=site.covs, numPrimary=dim(Tmin)[2])

#Null
(fm0=colext(psiformula=~1,gammaformula=~1,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B"))

#NULL + ELEVATION for initial occupancy
(fm0.1=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B"))

#Elevation
(fm1=colext(psiformula=~Elevation,gammaformula=~Elevation,epsilonformula=~Elevation,pformula=~1,data=umf,method="L-BFGS-B"))

(fm1.1=colext(psiformula=~Elevation,gammaformula=~Elevation,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B"))

(fm1.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Elevation,pformula=~1,data=umf,method="L-BFGS-B"))


#Tmin
(fm2=colext(psiformula=~Elevation,gammaformula=~Tmin,epsilonformula=~Tmin,pformula=~1,data=umf,method="SANN"))

(fm2.1=colext(psiformula=~Elevation,gammaformula=~Tmin,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B"))

(fm2.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmin,pformula=~1,data=umf,method="L-BFGS-B"))



#Tmax
(fm3=colext(psiformula=~Elevation,gammaformula=~Tmax,epsilonformula=~Tmax,pformula=~1,data=umf,method="L-BFGS-B"))

(fm3.1=colext(psiformula=~Elevation,gammaformula=~Tmax,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B"))

(fm3.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tmax,pformula=~1,data=umf,method="L-BFGS-B"))



#Tvar
(fm4=colext(psiformula=~Elevation,gammaformula=~Tvar,epsilonformula=~Tvar,pformula=~1,data=umf,method="L-BFGS-B"))

(fm4.1=colext(psiformula=~Elevation,gammaformula=~Tvar,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B"))

(fm4.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~Tvar,pformula=~1,data=umf,method="L-BFGS-B"))


#ForestLoss
(fm5=colext(psiformula=~Elevation,gammaformula=~ForestLoss,epsilonformula=~ForestLoss,pformula=~1,data=umf,method="L-BFGS-B"))

(fm5.1=colext(psiformula=~Elevation,gammaformula=~ForestLoss,epsilonformula=~1,pformula=~1,data=umf,method="L-BFGS-B"))

(fm5.2=colext(psiformula=~Elevation,gammaformula=~1,epsilonformula=~ForestLoss,pformula=~1,data=umf,method="L-BFGS-B"))


#Time-dependent
#(fm6=colext(psiformula=~1,gammaformula=~year-1,epsilonformula=~year-1,pformula=~1,data=umf,method="BFGS"))

#Time+Tmin
#(fm7=colext(psiformula=~1,gammaformula=~(year-1)+Tmin,epsilonformula=~(year-1)+Tmin,pformula=~1,data=umf,method="L-BFGS-B"))

#Time+tmax
#(fm8=colext(psiformula=~1,gammaformula=~(year-1)+Tmax,epsilonformula=~(year-1)+Tmax,pformula=~1,data=umf,method="BFGS"))


#Time+tvar
#(fm9=colext(psiformula=~1,gammaformula=~(year-1)+Tvar,epsilonformula=~(year-1)+Tvar,pformula=~1,data=umf,method="BFGS"))

#Time+Elevation
#(fm10=colext(psiformula=~1,gammaformula=~(year-1)+Elevation,epsilonformula=~(year-1)+Elevation,pformula=~1,data=umf,method="SANN"))

#Time+ForestLoss
#(fm11=colext(psiformula=~1,gammaformula=~(year-1)+ForestLoss,epsilonformula=~(year-1)+ForestLoss,pformula=~1,data=umf,method="L-BFGS-B"))

models=fitList(
'fm0_psi(.),gam(.),eps(.),p(.)'=fm0,
'fm0_psi(Elev),gam(.),eps(.),p(.)'=fm0.1,
'fm1_psi(Elev),gam(Elev),eps(Elev),p(.)'=fm1,
'fm1.1_psi(Elev),gam(Elev),eps(.),p(.)'=fm1.1,
'fm1.2_psi(Elev),gam(.),eps(Elev),p(.)'=fm1.2,
'fm2_psi(Elev),gam(Tmin),eps(Tmin),p(.)'=fm2,
'fm2.1_psi(Elev),gam(Tmin),eps(.),p(.)'=fm2.1,
'fm2.2_psi(Elev),gam(.),eps(Tmin),p(.)'=fm2.2,
'fm3_psi(Elev),gam(Tmax),eps(Tmax),p(.)'=fm3,
'fm3.1_psi(Elev),gam(Tmax),eps(.),p(.)'=fm3.1,
'fm3.2_psi(Elev),gam(.),eps(Tmax),p(.)'=fm3.2,
'fm4_psi(Elev),gam(Tvar),eps(Tvar),p(.)'=fm4,
'fm4.1_psi(Elev),gam(Tvar),eps(.),p(.)'=fm4.1,
'fm4.2_psi(Elev),gam(.),eps(Tvar),p(.)'=fm4.2,
'fm5_psi(Elev),gam(ForestLoss),eps(ForestLoss),p(.)'=fm5,
'fm5.1_psi(Elev),gam(ForestLoss),eps(.),p(.)'=fm5.1,
'fm5.2_psi(Elev),gam(.),eps(ForestLoss),p(.)'=fm5.2,
#'fm6_psi(Elev),gam(Time),eps(Time),p(.)'=fm6,
#'fm6.1_psi(Elev),gam(Time),eps(.),p(.)'=fm6.1,
#'fm6.2_psi(Elev),gam(.),eps(Time),p(.)'=fm6.2,
#'fm7_psi(Elev),gam(Time+Timin),eps(Time+Tmin),p(.)'=fm7,
#'fm7.1_psi(Elev),gam(Time+Tmin),eps(.),p(.)'=fm7.1,
#'fm7.2_psi(Elev),gam(.),eps(Time+Tmin),p(.)'=fm7.2,
#'fm8_psi(Elev),gam(Time+Timax),eps(Time+Tmax),p(.)'=fm8,
#'fm8.1_psi(Elev),gam(Time+Tmax),eps(.),p(.)'=fm8.1,
#'fm8.2_psi(Elev),gam(.),eps(Time+Tmax),p(.)'=fm8.2,
#'fm9_psi(Elev),gam(Time+Tvar),eps(Time+Tvar),p(.)'=fm9
#'fm9.1_psi(Elev),gam(Time+Tvar),eps(.),p(.)'=fm9.1,
#'fm9.2_psi(Elev),gam(.),eps(Time+Tvar),p(.)'=fm9.2,
#'fm10_psi(Elev),gam(Time+Elev),eps(Time+Elev),p(.)'=fm10
#'fm10.1_psi(Elev),gam(Time+Elev),eps(.),p(.)'=fm10.1,
#'fm10.2_psi(Elev),gam(.),eps(Time+Elev),p(.)'=fm10.2,
#'fm11_psi(Elev),gam(Time+ForestLoss),eps(Time+ForestLoss),p(.)'=fm11,
#'fm11.1_psi(Elev),gam(Time+ForestLoss),eps(.),p(.)'=fm11.1,
#'fm11.2_psi(Elev),gam(.),eps(Time+ForestLoss),p(.)'=fm11.2
)

(ms <- modSel(models))


##############################
#FIGURE WITH BEST MODEL
##############################

ndm=data.frame(Tvar=seq(min(Tvar),max(Tvar),length=50))

E.eps=predict(fm4.2,type="ext",newdata=ndm,appendData=TRUE)

x=ndm$Tvar
y=E.eps$Predicted
y.Err=1.96*E.eps$SE
y.Up = y + y.Err; y.Dn =y-y.Err

plot(x,y,type="l",lty=1,ylim=c(0,0.6),xlim=c(-1.2,5),col="red",xlab="Scaled(Tvar)",ylab=expression(bold(hat(epsilon)),cex.lab=1.0,main="Tvar"))

polygon(c(x,rev(x)),c(y.Up,rev(y.Dn)),col=rgb(255,0,0,25,maxColorValue=255),border=rgb(255,0,0,25,maxColorValue=255))
# Use the rev commands so the border moves logically around the shaded area
lines(x,y,type="l",lty=1,col="red") # put the means back on top of the polygon 

 