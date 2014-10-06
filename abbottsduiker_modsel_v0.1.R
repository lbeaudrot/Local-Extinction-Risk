library(unmarked)

#laptop###################
#abbottsduiker=read.csv("~/Desktop/Working_folder/UF/DB/Dropbox/TEAM Extinction Risk/VB Species Data Files/VB.abbottsduiker_Panthera_onca.csv")

#laptop
#covs=read.csv("~/Desktop/Working_folder/UF/DB/Dropbox/TEAM Extinction Risk/VB Species Data Files/VB_Covariate_Data.csv")

# #iMac
 abbottsduiker=read.csv("~/Dropbox/TEAM Extinction Risk/UDZ Species Data Files/UDZ.Abbottsduiker_Cephalophus_spadix.csv")

 covs=read.csv("~/Dropbox/TEAM Extinction Risk/UDZ Species Data Files/UDZ.CovariateData_w_Interpolated.csv")


y.abbottsduiker=abbottsduiker[,2:616]
########################################################
#NOTE
#Number of samples in secondary periods is unbalanced
#Randomly remove excess of 100 samples
########################################################
set.seed(67)

# y.abbottsduiker.08=as.matrix(y.abbottsduiker[,1:111]) #total 111
# rem=round(runif(11,1,111))
# y.abbottsduiker.08=y.abbottsduiker.08[,-rem]

y.abbottsduiker.09=as.matrix(y.abbottsduiker[,1:124]) #total 124
rem=round(runif(4,1,124))
y.abbottsduiker.09=y.abbottsduiker.09[,-rem]

y.abbottsduiker.10=as.matrix(y.abbottsduiker[,125:247]) #total 123
rem=round(runif(3,1,123))
y.abbottsduiker.10=y.abbottsduiker.10[,-rem]

y.abbottsduiker.11=as.matrix(y.abbottsduiker[,248:370]) #total 123
rem=round(runif(3,1,123))
y.abbottsduiker.11=y.abbottsduiker.11[,-rem]

y.abbottsduiker.12=as.matrix(y.abbottsduiker[,371:494]) #total 124
rem=round(runif(4,1,124))
y.abbottsduiker.12=y.abbottsduiker.12[,-rem]

y.abbottsduiker.13=as.matrix(y.abbottsduiker[,495:615]) #total 121
rem=runif(1,1,121)
y.abbottsduiker.13=y.abbottsduiker.13[,-rem]

# y.abbottsduiker.14=as.matrix(y.abbottsduiker[,633:734]) #total 102
# rem=runif(2,1,102)
# y.abbottsduiker.14=y.abbottsduiker.14[,-rem]


y.abbottsduiker=cbind(y.abbottsduiker.09,y.abbottsduiker.10,y.abbottsduiker.11,y.abbottsduiker.12,y.abbottsduiker.13)

years=as.character(2010:2014)
years=matrix(years,nrow(abbottsduiker),5,byrow=TRUE)

Tmin=as.matrix((covs[,6:10]))
Tmax=as.matrix((covs[,11:15]))
Tvar=as.matrix((covs[,16:20]))

Elevation=as.vector((covs[,2]))

ForestLoss=as.vector((covs[,3]))

###############################
#Scale and centering covariates
###############################
Elevation=scale(Elevation)

ForestLoss=scale(ForestLoss)

Tmin=scale(Tmin)
Tmin[is.na(Tmin)]<-0

Tmax=scale(Tmax)
Tmax[is.na(Tmax)]<-0

Tvar=scale(Tvar)
Tvar[is.na(Tvar)]<-0

site.covs<-data.frame(Elevation,ForestLoss)



umf<-unmarkedMultFrame(y=y.abbottsduiker,yearlySiteCovs=list(year=years,Tmin=Tmin,Tmax=Tmax,Tvar=Tvar),siteCovs=site.covs,numPrimary=5)

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
(fm6=colext(psiformula=~1,gammaformula=~year-1,epsilonformula=~year-1,pformula=~1,data=umf,method="BFGS"))

#Time+Tmin
(fm7=colext(psiformula=~1,gammaformula=~(year-1)+Tmin,epsilonformula=~(year-1)+Tmin,pformula=~1,data=umf,method="L-BFGS-B"))

#Time+tmax
(fm8=colext(psiformula=~1,gammaformula=~(year-1)+Tmax,epsilonformula=~(year-1)+Tmax,pformula=~1,data=umf,method="BFGS"))


#Time+tvar
(fm9=colext(psiformula=~1,gammaformula=~(year-1)+Tvar,epsilonformula=~(year-1)+Tvar,pformula=~1,data=umf,method="BFGS"))

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
#'fm5.2_psi(Elev),gam(.),eps(ForestLoss),p(.)'=fm5.2,
'fm6_psi(Elev),gam(Time),eps(Time),p(.)'=fm6,
#'fm6.1_psi(Elev),gam(Time),eps(.),p(.)'=fm6.1,
#'fm6.2_psi(Elev),gam(.),eps(Time),p(.)'=fm6.2,
'fm7_psi(Elev),gam(Time+Timin),eps(Time+Tmin),p(.)'=fm7,
#'fm7.1_psi(Elev),gam(Time+Tmin),eps(.),p(.)'=fm7.1,
#'fm7.2_psi(Elev),gam(.),eps(Time+Tmin),p(.)'=fm7.2,
'fm8_psi(Elev),gam(Time+Timax),eps(Time+Tmax),p(.)'=fm8,
#'fm8.1_psi(Elev),gam(Time+Tmax),eps(.),p(.)'=fm8.1,
#'fm8.2_psi(Elev),gam(.),eps(Time+Tmax),p(.)'=fm8.2,
'fm9_psi(Elev),gam(Time+Tvar),eps(Time+Tvar),p(.)'=fm9
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

 