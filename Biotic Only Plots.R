# Plot partial relationships from colext models for populations with significant biotic top models (but not interactions)
# Show colonization and extinction probability as a function of covariates with confidence intervals



# Functions to plot biotic estimates + standard errors for colonization and extinction 
plot.col <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Biotic, pch=19, las=1, ylab="Colonization", xlab="", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Biotic, rev(outputdata$Biotic)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Biotic, type="line", lwd=2)
}

plot.ext <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Biotic, pch=19, las=1, ylab="Extinction", xlab="", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Biotic, rev(outputdata$Biotic)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Biotic, type="line", lwd=2)
    }


############### BEGIN PLOTTING #######################
library(fields)
set.panel(4,3)
par(mar=c(3,4,1,0))


#### COLONIZATION ####  

# Galidia elegans (RNF) nms[57]
# ~1 ~Biotic ~1 ~1 
nmsK <- 57
fm6.1=colext(psiformula=~1,
                    gammaformula=~Biotic,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm6.1   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60))
plot.col(inputmodel, inputdata)
mtext("Galidia elegans (RNF)", side=3, line=0, cex=0.7)


# Tapirus terrestris (YAN) nms[49]
# ~1 ~Biotic ~1 ~1 
nmsK <- 49
fm6.1=colext(psiformula=~1,
                    gammaformula=~Biotic,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm6.1   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60))
plot.col(inputmodel, inputdata)
mtext("Tapirus terrestris (YAN)", side=3, line=0, cex=0.7)


# Atherurus macrourus (PSH) nms[30]
# ~Elevation ~Biotic ~1 ~1 
nmsK <- 30
fm6.4=colext(psiformula=~Elevation,
                    gammaformula=~Biotic,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm6.4   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60))
plot.col(inputmodel, inputdata)
mtext("Atherurus macrourus (PSH)", side=3, line=0, cex=0.7)


# Mazama americana (YAN) nms[45]
# ~Elevation ~Biotic ~Biotic ~1 
nmsK <- 45
fm6.3=colext(psiformula=~Elevation,
                    gammaformula=~Biotic,
                    epsilonformula=~Biotic,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm6.3   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60))
plot.col(inputmodel, inputdata)
mtext("Mazama americana (YAN)", side=3, line=0, cex=0.7)


# Tragulus kanchil (NAK) nms[54]
# ~Elevation ~Tmax + Biotic ~1 ~1
nmsK <- 54
fm9.4=colext(psiformula=~Elevation,
                    gammaformula=~Tmax + Biotic,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm9.4   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60), Tmax=rep(0, length=60))
plot.col(inputmodel, inputdata)
mtext("Tragulus kanchil (NAK)", side=3, line=0, cex=0.7)

# Cercocebus sanjei (UDZ) nms[11]
# ~1 ~Biotic + Tmin ~1 ~1
nmsK <- 11
fm7.4=colext(psiformula=~1,
                    gammaformula=~Biotic + Tmin,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm7.4   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60), Tmin=rep(0, length=60))
plot.col(inputmodel, inputdata)
mtext("Cercocebus sanjei (UDZ)", side=3, line=0, cex=0.7)


# Cercocebus lhoesti (BIF) (UDZ) nms[24]
# ~Elevation ~Tvar + Biotic ~Tvar + Biotic ~1
nmsK <- 24
fm13.3=colext(psiformula=~Elevation,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm13.3   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60), Tvar=rep(0, length=60))
plot.col(inputmodel, inputdata)
mtext("Cercocebus lhoesti (BIF)", side=3, line=0, cex=0.7)

################################# 
########## EXTINCTION ###########
################################# 

# Bdeogale crassicauda (UDZ) nms[8]
# ~1 ~1 ~Biotic + Tmin ~1
nmsK <- 8
fm7.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm7.2   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60), Tmin=rep(0, length=60))
plot.ext(inputmodel, inputdata)
mtext("Bdeogale crassicauda (UDZ)", side=3, line=0, cex=0.7)


# Nasua nasua (YAN) nms[46]
# ~Elevation ~1 ~Biotic + Tmin ~1
nmsK <- 46
fm7.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm7.5   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60), Tmin=rep(0, length=60))
plot.ext(inputmodel, inputdata)
mtext("Nasua nasua (YAN)", side=3, line=0, cex=0.7)


# Potamochoerus larvatus (BIF) [29]
# ~Elevation ~1 ~Tmax + Biotic ~1
nmsK <- 29
fm9.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tmax + Biotic,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm9.5   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60), Tmax=rep(0, length=60))
plot.ext(inputmodel, inputdata)
mtext("Potamochoerus larvatus (BIF)", side=3, line=0, cex=0.7)


# Cercocebus lhoesti (BIF) (UDZ) nms[24]
# ~Elevation ~Tvar + Biotic ~Tvar + Biotic ~1
nmsK <- 24
fm13.3=colext(psiformula=~Elevation,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm13.3   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60), Tvar=rep(0, length=60))
plot.ext(inputmodel, inputdata)
mtext("Cercocebus lhoesti (BIF)", side=3, line=0, cex=0.7)



# Mazama americana (YAN) nms[45]
# ~Elevation ~Biotic ~Biotic ~1 
nmsK <- 45
fm6.3=colext(psiformula=~Elevation,
                    gammaformula=~Biotic,
                    epsilonformula=~Biotic,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm6.3   
inputdata <- data.frame(Biotic=seq(min(Biotic, na.rm=TRUE),max(Biotic, na.rm=TRUE),length=60))
plot.ext(inputmodel, inputdata)
mtext("Mazama americana (YAN)", side=3, line=0, cex=0.7)

 