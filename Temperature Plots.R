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

plot.colMean <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="col", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmean, pch=19, las=1, ylab="Colonization", xlab="Mean", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmean, rev(outputdata$Tmean)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tmean, type="line", lwd=2)
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

plot.extMean <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tmean, pch=19, las=1, ylab="Extinction", xlab="Mean", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tmean, rev(outputdata$Tmean)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tmean, type="line", lwd=2)
}

plot.extTsd <- function(inputmodel, inputdata){
    outputdata <- predict(inputmodel, type="ext", newdata=inputdata, appendData=TRUE)
    plot(outputdata$Predicted ~ outputdata$Tsd, pch=19, las=1, ylab="Extinction", xlab="Standard deviation", bty="n", xlim=c(), ylim=c(0,1), cex=0.8, col="white")
    polygon(c(outputdata$Tsd, rev(outputdata$Tsd)),
        c(ifelse(outputdata$Predicted + 1.96*outputdata$SE>1, 1, outputdata$Predicted + 1.96*outputdata$SE), 
          ifelse(rev(outputdata$Predicted - 1.96*outputdata$SE)<0, 0, rev(outputdata$Predicted - 1.96*outputdata$SE))),
          col="gray75")
    lines(outputdata$Predicted ~ outputdata$Tsd, type="line", lwd=2)
    }


# Species to model
# Colonization
# Caracal aurata (BIF) nms[21]
    # ~1 ~ Tmin ~1 ~1
# Cercocebus sanjei (UDZ) nms[11]
    # ~1 ~ Biotic + Tmin ~1 ~1
# Cuniculus paca (YAN) nms[40]
     # ~1 ~Tmin ~Tmin ~1
# Dasyprocta punctata (VB) nms[2]
     # ~1 ~Tmin ~1 ~1
# Paraxerus vexillarius  (UDZ) nms[18]
     # ~1 ~Tmin ~1 ~1
# Pecari tajacu (VB_) nms[6]
     # ~1 ~Tvar+Biotic ~1 ~1


#### Extinction 

# Bdeogale crassicauda (UDZ) nms[8]
    # ~1 ~1 ~Biotic + Tmin ~1
# Cephalophus nigrifrons (BIF) nms[22]
    # ~1 ~1 ~Tmean + Tsd ~1
# Cephalophus spadix (UDZ) nms[10]
    # ~1 ~1 ~Tvar + Biotic ~1
# Cricetomys gambianus  (UDZ) nms[13]
     # ~1 ~1 ~Tmax ~1
# Cuniculus paca (YAN) nms[40]
     # ~1 ~Tmin ~Tmin ~1
# Potamochoerus larvatus (UDZ) nms[19]
     # ~1 ~1 ~Tvar ~1
# Tragulus kanchil (NAK) nms[54]
     # ~1 ~1 ~Tmean + Tsd ~1


#### BEGIN PLOTTING ########

library(fields)
set.panel(4,3)
par(mar=c(4,4,3,0))

#### Colonization 
# Caracal aurata (BIF) nms[21]
    # ~1 ~ Tmin ~1 ~1
nmsK <- 21

try((fm2.1=colext(psiformula=~1,
                  gammaformula=~Tmin,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.1")) {
  if(CondNum(fm2.1)<2000){
    if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
} 
}}
inputmodel <- fm2.1          
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Caracal aurata (BIF)", side=3, line=0, cex=0.7)


# Cercocebus sanjei (UDZ) nms[11]
    # ~1 ~ Biotic + Tmin ~1 ~1
nmsK <- 11

try((fm28.1=colext(psiformula=~1,
                   gammaformula=~Biotic + Tmin,
                   epsilonformula=~1,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm28.1")) {
  if(CondNum(fm28.1)<2000){
    if(CondNum(fm28.1)>0){mods=c(mods,fm28.1)}
} 
}}
inputmodel <- fm28.1
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Cercocebus sanjei (UDZ)", side=3, line=0, cex=0.7)


# Cuniculus paca (YAN) nms[40]
     # ~1 ~Tmin ~Tmin ~1
#nmsK <- 40

#try((fm2=colext(psiformula=~1,
#                gammaformula=~Tmin,
#                epsilonformula=~Tmin,
#                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm2")) {
#  if(CondNum(fm2)<2000){
#    if(CondNum(fm2)>0){mods=c(mods,fm2)}
#} 
#}}
#inputmodel <- fm2
#inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60))
#plot.colTmin(inputmodel, inputdata)
#mtext("Cuniculus paca (YAN)", side=3, line=0, cex=0.7)


# Dasyprocta punctata (VB) nms[2]
     # ~1 ~Tmin ~1 ~1
nmsK <- 2

try((fm2.1=colext(psiformula=~1,
                  gammaformula=~Tmin,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.1")) {
  if(CondNum(fm2.1)<2000){
    if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
} 
}}
inputmodel <- fm2.1
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Dasyprocta punctata (VB)", side=3, line=0, cex=0.7)


# Paraxerus vexillarius  (UDZ) nms[18]
     # ~1 ~Tmin ~1 ~1
nmsK <- 18

try((fm2.1=colext(psiformula=~1,
                  gammaformula=~Tmin,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm2.1")) {
  if(CondNum(fm2.1)<2000){
    if(CondNum(fm2.1)>0){mods=c(mods,fm2.1)}
} 
}}
inputmodel <- fm2.1
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60))
plot.colTmin(inputmodel, inputdata)
mtext("Paraxerus vexillarius  (UDZ)", side=3, line=0, cex=0.7)






#### Extinction 



# Bdeogale crassicauda (UDZ) nms[8]
    # ~1 ~1 ~Biotic + Tmin ~1
nmsK <- 8

try((fm28.2=colext(psiformula=~1,
                   gammaformula=~1,
                   epsilonformula=~Biotic + Tmin,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm28.2")) {
  if(CondNum(fm28.2)<2000){
    if(CondNum(fm28.2)>0){mods=c(mods,fm28.2)}
} 
}}
inputmodel <- fm28.2
inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.extTmin(inputmodel, inputdata)
mtext("Bdeogale crassicauda (UDZ)", side=3, line=0, cex=0.7)


# Cuniculus paca (YAN) nms[40]
     # ~1 ~Tmin ~Tmin ~1
#nmsK <- 40

#try((fm2=colext(psiformula=~1,
#                gammaformula=~Tmin,
#                epsilonformula=~Tmin,
#                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

#if(exists("fm2")) {
#  if(CondNum(fm2)<2000){
#    if(CondNum(fm2)>0){mods=c(mods,fm2)}
#} 
#}}
#inputmodel <- fm2
#inputdata <- data.frame(Tmin=seq(min(Tmin, na.rm=TRUE),max(Tmin, na.rm=TRUE),length=60))
#plot.extTmin(inputmodel, inputdata)
#mtext("Cuniculus paca (YAN)", side=3, line=0, cex=0.7)


# Cricetomys gambianus  (UDZ) nms[13]
     # ~1 ~1 ~Tmax ~1
nmsK <- 13

try((fm6.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmax,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm6.2")) {
  if(CondNum(fm6.2)<2000){
    if(CondNum(fm6.2)>0){mods=c(mods,fm6.2)}
} 
}}
inputmodel <- fm6.2
inputdata <- data.frame(Tmax=seq(min(Tmax, na.rm=TRUE),max(Tmax, na.rm=TRUE),length=60))
plot.extTmax(inputmodel, inputdata)
mtext("Cricetomys gambianus (UDZ)", side=3, line=0, cex=0.7)



# Pecari tajacu (VB_) nms[6]
     # ~1 ~Tvar+Biotic ~1 ~1
nmsK <- 6

try((fm42.1=colext(psiformula=~1,
                   gammaformula=~Tvar + Biotic,
                   epsilonformula=~1,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm42.1")) {
  if(CondNum(fm42.1)<2000){
    if(CondNum(fm42.1)>0){mods=c(mods,fm42.1)}
} 
}}
inputmodel <- fm42.1
inputdata <- data.frame(Tvar=seq(min(Tvar, na.rm=TRUE),max(Tvar, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.colTvar(inputmodel, inputdata)
mtext("Pecari tajacu (VB)", side=3, line=0, cex=0.7)


# Cephalophus spadix (UDZ) nms[10]
    # ~1 ~1 ~Tvar + Biotic ~1
nmsK <- 10

try((fm42.2=colext(psiformula=~1,
                   gammaformula=~1,
                   epsilonformula=~Tvar + Biotic,
                   pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm42.2")) {
  if(CondNum(fm42.2)<2000){
    if(CondNum(fm42.2)>0){mods=c(mods,fm42.2)}
} 
}}
inputmodel <- fm42.2
inputdata <- data.frame(Tvar=seq(min(Tvar, na.rm=TRUE),max(Tvar, na.rm=TRUE),length=60), Biotic=rep(0, length=60))
plot.extTvar(inputmodel, inputdata)
mtext("Cephalophus spadix (UDZ)", side=3, line=0, cex=0.7)


# Potamochoerus larvatus (UDZ) nms[19]
     # ~1 ~1 ~Tvar ~1
nmsK <- 19

try((fm7.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tvar,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm7.2")) {
  if(CondNum(fm7.2)<2000){
    if(CondNum(fm7.2)>0){mods=c(mods,fm7.2)}
} 
}}
inputmodel <- fm7.2
inputdata <- data.frame(Tvar=seq(min(Tvar, na.rm=TRUE),max(Tvar, na.rm=TRUE),length=60))
plot.extTvar(inputmodel, inputdata)
mtext("Potamochoerus larvatus (UDZ)", side=3, line=0, cex=0.7)


# Cephalophus nigrifrons (BIF) nms[22]
    # ~1 ~1 ~Tmean + Tsd ~1
nmsK <- 22

try((fm3.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmean + Tsd,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.2")) {
  if(CondNum(fm3.2)<2000){
    if(CondNum(fm3.2)>0){mods=c(mods,fm3.2)}
} 
}}
inputmodel <- fm3.2
#inputdata <- data.frame(Tmean=seq(min(Tmean, na.rm=TRUE),max(Tmean, na.rm=TRUE),length=60), Tsd=rep(0, length=60))
#plot.extMean(inputmodel, inputdata)
inputdata <- data.frame(Tsd=seq(min(Tsd, na.rm=TRUE),max(Tsd, na.rm=TRUE),length=60), Tmean=rep(0, length=60))
plot.extTsd(inputmodel, inputdata)
mtext("Cephalophus nigrifrons (BIF)", side=3, line=0, cex=0.7)





# Tragulus kanchil (NAK) nms[54]
     # ~1 ~1 ~Tmean + Tsd ~1
nmsK <- 54

try((fm3.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Tmean + Tsd,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm3.2")) {
  if(CondNum(fm3.2)<2000){
    if(CondNum(fm3.2)>0){mods=c(mods,fm3.2)}
} 
}}
inputmodel <- fm3.2
#inputdata <- data.frame(Tmean=seq(min(Tmean, na.rm=TRUE),max(Tmean, na.rm=TRUE),length=60), Tsd=rep(0, length=60))
#plot.extMean(inputmodel, inputdata)
inputdata <- data.frame(Tsd=seq(min(Tsd, na.rm=TRUE),max(Tsd, na.rm=TRUE),length=60), Tmean=rep(0, length=60))
plot.extTsd(inputmodel, inputdata)
mtext("Tragulus kanchil (NAK)", side=3, line=0, cex=0.7)


