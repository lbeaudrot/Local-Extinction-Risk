# Plots for elevation shift species with interaction terms in top TB model

#### Colonization 
#Pan troglodytes BIF
#Eira barbara YAN
#Mazama americana YAN
#Dasypus novemcinctus VB

plot.colTmax <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmax, pch=19, las=1, ylab="", xlab="", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tmax, type="line", lwd=2)
}

extra.colTmax <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tmax, type="line", lwd=2)
}

plot.colTvar <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="", xlab="", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}

extra.colTvar <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}

# light green: col=rgb(0,200,100,75,maxColorValue=255)

############### BEGIN PLOTTING #######################
library(fields)
set.panel(4,1)
par(mar=c(2,1,1,1))

# Pan troglodytes (BIF) nms[27]
    # ~1 ~Tvar * Biotic ~1 ~1
    # BIF.PanTB <- fm40.1
nmsK <- 27
try((fm40.1=colext(psiformula=~1,
                   gammaformula=~Tvar * Biotic,
                   epsilonformula=~1,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm40.1")) {
  if(CondNum(fm40.1)<2000){
    if(CondNum(fm40.1)>0){mods=c(mods,fm40.1)}
} 
}}
inputmodel <- fm40.1
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTvar(inputmodel, inputdata)
#inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
#plot.colTvar(inputmodel, inputdata)
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTvar(inputmodel, inputdata)

#mtext("Pan troglodytes (BIF)", side=3, line=0, cex=0.7)


# ~1 ~Tmax*Biotic ~1 ~1

# Dasypus novemcinctus (VB) nms[3]
nmsK <- 3

# Mazama americana (YAN) nms[45]
nmsK <- 45
# Eira barbara (YAN)  nms[43]
nmsK <- 43


#mtext("Eira barbara (YAN)", side=3, line=1, cex=0.7)
#mtext("Mazama americana (YAN)", side=3, line=1, cex=0.7)
#mtext("Dasypus novemcinctus (VB)", side=3, line=1, cex=0.7)

for(k in nmsK:nmsK){
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
PCA1 <- data.frame(pca_covs)

to=dim(Tmin)[2]
years=as.character(1:to)
years=matrix(years,nrow(species),to,byrow=TRUE)

# ADD BIOTIC TO UMF COVARIATES
site.covs<-data.frame(Elevation, ForestLossCT, ForestGainCT, Biotic)

umf<-unmarkedMultFrame(y=species, yearlySiteCovs=list(Tmin=Tmin,Tmax=Tmax,Tvar=Tvar,Tsd=Tsd,Tmean=Tmean, PCA1=PCA1), siteCovs=site.covs, numPrimary=dim(Tmin)[2])

mods=list()



try((fm39.1=colext(psiformula=~1,
                   gammaformula=~Tmax * Biotic,
                   epsilonformula=~1,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm39.1")) {
  if(CondNum(fm39.1)<2000){
    if(CondNum(fm39.1)>0){mods=c(mods,fm39.1)}
} 
}}
inputmodel <- fm39.1
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
#inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(0, length=60), Int=rep(0, length=60))
#plot.colTmax(inputmodel, inputdata)
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTmax(inputmodel, inputdata)


