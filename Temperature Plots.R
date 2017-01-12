# Revised for Nature Climate Change top models 
# Temperature only Plots (excludes species with interaction term in top model)

# Helper functions
plot.colTmax <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmax, pch=19, las=1, ylab="Colonization", xlab="Maximum", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tmax, type="line", lwd=2)
}

plot.colTmin <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, xlab="Minimum", ylab="Colonization", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

plot.colTvar <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="Colonization", xlab="Variance", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}


plot.extTmax <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmax, pch=19, las=1, ylab="Extinction", xlab="Maximum", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmax, rev(outputdata$Tmax)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tmax, type="line", lwd=2)
}

plot.extTmin <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, ylab="Extinction", xlab="Minimum", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

plot.extTvar <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="Extinction", xlab="Variance", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}



# Species to model

# Colonization

# Pecari tajacu (VB_) nms[6]
# ~Elevation ~Tvar ~1 ~1

# Potamochoerus larvatus (UDZ) nms[19]
# ~1 ~Tvar ~1 ~1

# Tragulus kanchil (NAK) nms[54]
# ~Elevation ~Tmax + Biotic ~1 ~1

# Cercocebus sanjei (UDZ) nms[11]
# ~1 ~Biotic + Tmin ~1 ~1

# Cercocebus lhoesti (BIF) (UDZ) nms[24]
# ~Elevation ~Tvar + Biotic ~Tvar + Biotic ~1

#### Extinction 

# Cricetomys gambianus  (UDZ) nms[13]
# ~1 ~1 ~Tmax ~1

# Fossa fossana (RNF) nms[56]
# ~Elevation ~1 ~Tmax ~1

# Pan troglodytes (BIF) nms[27]
# ~Elevation ~1 ~Tmin ~1

# Bdeogale crassicauda (UDZ) nms[8]
# ~1 ~1 ~Biotic + Tmin ~1

# Nasua nasua (YAN) nms[46]
# ~Elevation ~1 ~Biotic + Tmin ~1

# Potamochoerus larvatus (BIF) [29]
# ~Elevation ~1 ~Tmax + Biotic ~1


#### BEGIN PLOTTING ########

library(fields)
set.panel(4,3)
par(mar=c(4,4,3,0))

#### Colonization 

# Pecari tajacu (VB_) nms[6]
# ~Elevation ~Tvar ~1 ~1
nmsK <- 6
fm5.4=colext(psiformula=~Elevation,
                    gammaformula=~Tvar,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm5.4          
inputdata <- data.frame(Tvar=seq(min(Tvar, na.rm=TRUE),max(Tvar, na.rm=TRUE),length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Pecari tajacu (VB)", side=3, line=0, cex=0.7)  


# Potamochoerus larvatus (UDZ) nms[19]
# ~1 ~Tvar ~1 ~1
nmsK <- 19
fm5.1=colext(psiformula=~1,
                    gammaformula=~Tvar,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm5.1          
inputdata <- data.frame(Tvar=seq(min(Tvar, na.rm=TRUE),max(Tvar, na.rm=TRUE),length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Potamochoerus larvatus (UDZ)", side=3, line=0, cex=0.7)  


# Tragulus kanchil (NAK) nms[54]
# ~Elevation ~Tmax + Biotic ~1 ~1
nmsK <- 54

fm9.4=colext(psiformula=~Elevation,
                    gammaformula=~Tmax + Biotic,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm9.4   
inputdata <- data.frame(Tmax=seq(min(Tmax, na.rm=TRUE),max(Tmax, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
mtext("Tragulus kanchil (NAK)", side=3, line=0, cex=0.7)


# Cercocebus sanjei (UDZ) nms[11]
# ~1 ~Biotic + Tmin ~1 ~1
nmsK <- 11
fm7.4=colext(psiformula=~Elevation,
                    gammaformula=~Biotic + Tmin,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm7.4   
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Cercocebus sanjei (UDZ)", side=3, line=0, cex=0.7)


# Cercocebus lhoesti (BIF) nms[24]
# ~Elevation ~Tvar + Biotic ~Tvar + Biotic ~1
nmsK <- 24
fm13.3=colext(psiformula=~Elevation,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm13.3   
inputdata <- data.frame(Tvar=seq(min(Tvar, na.rm=TRUE),max(Tvar, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Cercocebus lhoesti (BIF)", side=3, line=0, cex=0.7)




#########################
#### Extinction #########

# Cricetomys gambianus  (UDZ) nms[13]
# ~1 ~1 ~Tmax ~1
nmsK <- 13
fm3.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Tmax,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm3.2   
inputdata <- data.frame(Tmax=seq(min(Tmax, na.rm=TRUE),max(Tmax, na.rm=TRUE),length=60))
plot.extTmax(inputmodel, inputdata)
mtext("Cricetomys gambianus (UDZ)", side=3, line=0, cex=0.7)


# Fossa fossana (RNF) nms[56]
# ~Elevation ~1 ~Tmax ~1
nmsK <- 56
fm3.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tmax,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm3.5   
inputdata <- data.frame(Tmax=seq(min(Tmax, na.rm=TRUE),max(Tmax, na.rm=TRUE),length=60))
plot.extTmax(inputmodel, inputdata)
mtext("Fossa fossana (RNF)", side=3, line=0, cex=0.7)

# Pan troglodytes (BIF) nms[27]
# ~Elevation ~1 ~Tmin ~1
nmsK <- 27
fm2.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm2.5   
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Pan troglodytes (BIF)", side=3, line=0, cex=0.7)


# Bdeogale crassicauda (UDZ) nms[8]
# ~1 ~1 ~Biotic + Tmin ~1
nmsK <- 8
fm7.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm7.2   
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Bdeogale crassicauda (UDZ)", side=3, line=0, cex=0.7)

# Nasua nasua (YAN) nms[46]
# ~Elevation ~1 ~Biotic + Tmin ~1
nmsK <- 46
fm7.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Biotic + Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm7.5   
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Nasua nasua (YAN)", side=3, line=0, cex=0.7)

# Potamochoerus larvatus (BIF) [29]
# ~Elevation ~1 ~Tmax + Biotic ~1
nmsK <- 29
fm9.5=colext(psiformula=~Elevation,
                    gammaformula=~1,
                    epsilonformula=~Tmax + Biotic,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm9.5   
inputdata <- data.frame(Tmax=seq(min(Tmax, na.rm=TRUE),max(Tmax, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.extTmax(inputmodel, inputdata)
mtext("Potamochoerus larvatus (BIF)", side=3, line=0, cex=0.7)

# Cercocebus lhoesti (BIF) nms[24]
# ~Elevation ~Tvar + Biotic ~Tvar + Biotic ~1
nmsK <- 24
fm13.3=colext(psiformula=~Elevation,
                     gammaformula=~Tvar + Biotic,
                     epsilonformula=~Tvar + Biotic,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm13.3   
inputdata <- data.frame(Tvar=seq(min(Tvar, na.rm=TRUE),max(Tvar, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.extTvar(inputmodel, inputdata)
mtext("Cercocebus lhoesti (BIF)", side=3, line=0, cex=0.7)


