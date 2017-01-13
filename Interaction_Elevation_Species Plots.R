# Plots for elevation shift species with interaction terms in top TB model

# Mazama temama (VB)
# ~1 ~Biotic*Tmin ~1 ~1

# Dasypus novemcinctus (VB)
# ~1 ~Tvar*Biotic ~Tvar*Biotic ~1

# Caracal aurata (BIF)
# ~1 ~1 ~Biotic*Tmin ~1

# Cephalophus spadix (UDZ)
# ~1 ~1 ~Biotic*Tmin ~1

# Genetta servalina (UDZ)
# ~Elevation ~Tmax*Biotic ~1 ~1


plot.colTmax <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmax, pch=19, las=1, ylab="Colonization", xlab="Tmax", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
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
    plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="Colonization", xlab="Tvar", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
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


plot.extTvar <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tvar, pch=19, las=1, ylab="Extinction", xlab="Tvar", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}

extra.extTvar <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    polygon(c(outputdata$Tvar, rev(outputdata$Tvar)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tvar, type="line", lwd=2)
}


plot.colTmin <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, ylab="Colonization", xlab="Tmin", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

extra.colTmin <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

plot.extTmin <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmin, pch=19, las=1, ylab="Extinction", xlab="Tmin", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(100,0,170,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

extra.extTmin <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    polygon(c(outputdata$Tmin, rev(outputdata$Tmin)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col=rgb(0,200,100,75,maxColorValue=255))
    lines(outputdata$Predicted ~ outputdata$Tmin, type="line", lwd=2)
}

# light green: col=rgb(0,200,100,75,maxColorValue=255)

############### BEGIN PLOTTING #######################
library(fields)
set.panel(3,2)
par(mar=c(4,4,2,2))


# Genetta servalina (UDZ)
# ~Elevation ~Tmax*Biotic ~1 ~1
nmsK <- 14
fm10.4=colext(psiformula=~Elevation,
                     gammaformula=~Tmax * Biotic,
                     epsilonformula=~1,
                     pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm10.4
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTmax(inputmodel, inputdata)
inputdata <- data.frame(Tmax=seq(min(Tmax),max(Tmax),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTmax(inputmodel, inputdata)
mtext("Genetta servalina (UDZ)", side=3, line=0, cex=0.7)


# Mazama temama (VB)
# ~1 ~Biotic*Tmin ~1 ~1
nmsK <- 5
fm8.1=colext(psiformula=~1,
                    gammaformula=~Biotic * Tmin,
                    epsilonformula=~1,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm8.1
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTmin(inputmodel, inputdata)
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTmin(inputmodel, inputdata)
mtext("Mazama temama (VB)", side=3, line=0, cex=0.7)


# Caracal aurata (BIF)
# ~1 ~1 ~Biotic*Tmin ~1
nmsK <- 21
fm8.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic * Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm8.2
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.extTmin(inputmodel, inputdata)
mtext("Caracal aurata (BIF)", side=3, line=0, cex=0.7)


# Cephalophus spadix (UDZ)
# ~1 ~1 ~Biotic*Tmin ~1
nmsK <- 10
fm8.2=colext(psiformula=~1,
                    gammaformula=~1,
                    epsilonformula=~Biotic * Tmin,
                    pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm8.2
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
inputdata <- data.frame(Tmin=seq(min(Tmin),max(Tmin),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.extTmin(inputmodel, inputdata)
mtext("Cephalophus spadix (UDZ)", side=3, line=0, cex=0.7)


# Dasypus novemcinctus (VB)
# ~1 ~Tvar*Biotic ~Tvar*Biotic ~1
nmsK <- 3
fm12=colext(psiformula=~1,
                   gammaformula=~Tvar * Biotic,
                   epsilonformula=~Tvar * Biotic,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))
inputmodel <- fm12
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.colTvar(inputmodel, inputdata)
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.colTvar(inputmodel, inputdata)
mtext("Dasypus novemcinctus (VB)", side=3, line=0, cex=0.7)

inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(1, length=60), Int=rep(0, length=60))
plot.extTvar(inputmodel, inputdata)
inputdata <- data.frame(Tvar=seq(min(Tvar),max(Tvar),length=60), Biotic=rep(-1, length=60), Int=rep(0, length=60))
extra.extTvar(inputmodel, inputdata)
mtext("Dasypus novemcinctus (VB)", side=3, line=0, cex=0.7)

