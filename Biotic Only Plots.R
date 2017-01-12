# Plot partial relationships from colext models for populations with significant biotic top models (but not interactions)
# Show colonization and extinction probability as a function of covariates with confidence intervals

#### Colonization 
# Cercocebus sanjei (UDZ) nms[11]
    # ~1 ~ Biotic + Tmin ~1 ~1
    # UDZ.CerTB <- fm28.1
# Cercopithecus lhoesti (BIF) nms[24]
    # ~1 ~Biotic ~Biotic ~1
    # BIF.CerTB <- fm4
# Galidia elegans (RNF) nms[57]
    # ~1 ~Biotic ~1 ~1
    # RNF.GalTB <- fm4.1
# Pecari tajacu (VB_) nms[6]
    # ~1 ~Tvar+Biotic ~1 ~1
    # VB.PecTB <- fm42.1
# Potamochoerus larvatus (BIF) nms[19]
    # ~1 ~Biotic ~1 ~1
    # BIF.PotTB <- fm4.1


#### Extinction 

# Bdeogale crassicauda (UDZ) nms[8]
   # ~1 ~1 ~Biotic + Tmin ~1
   # UDZ.BdeTB <- fm28.2
# Cephalophus spadix (UDZ) nms[10]
     # ~1 ~1 ~Tvar + Biotic ~1
    # UDZ.CepTB <- fm42.2
# Cercopithecus lhoesti (BIF) nms[24]
   # ~1 ~Biotic ~Biotic ~1
   # BIF.CerTB <- fm4
# Mazama temama (VB_) nms[5]
    # ~1 ~1 ~Biotic ~1
    # VB.MazTB <- fm4.2



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
set.panel(3,3)
par(mar=c(3,4,1,0))





# Potamochoerus larvatus (BIF) nms[19]
    # ~1 ~Biotic ~1 ~1
    # BIF.PotTB <- fm4.1
nmsK <- 19
try((fm4.1=colext(psiformula=~1,
                  gammaformula=~Biotic,
                  epsilonformula=~1,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
if(exists("fm4.1")) {
  if(CondNum(fm4.1)<2000){
    if(CondNum(fm4.1)>0){mods=c(mods,fm4.1)}
} 
}}  
inputmodel <- fm4.1
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60))
plot.col(inputmodel, inputdata)
mtext("Potamochoerus larvatus (BIF)", side=3, line=0, cex=0.7)





# Galidia elegans (RNF) nms[57]
    # ~1 ~Biotic ~1 ~1
    # RNF.GalTB <- fm4.1
nmsK <- 57
try((fm4=colext(psiformula=~1,
                gammaformula=~Biotic,
                epsilonformula=~Biotic,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
if(exists("fm4")) {
  if(CondNum(fm4)<2000){
    if(CondNum(fm4)>0){mods=c(mods,fm4)}
} 
}}   
inputmodel <- fm4.1
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60))
plot.col(inputmodel, inputdata)
mtext("Galidia elegans (RNF)", side=3, line=0, cex=0.7)




# Cercopithecus lhoesti (BIF) nms[24]
    # ~1 ~Biotic ~Biotic ~1
    # BIF.CerTB <- fm4
nmsK <- 24
try((fm4=colext(psiformula=~1,
                gammaformula=~Biotic,
                epsilonformula=~Biotic,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
if(exists("fm4")) {
  if(CondNum(fm4)<2000){
    if(CondNum(fm4)>0){mods=c(mods,fm4)}
} 
}}   
inputmodel <- fm4
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60))
plot.col(inputmodel, inputdata)
mtext("Cercopithecus lhoesti (BIF)", side=3, line=0, cex=0.7)




# Pecari tajacu (VB_) nms[6]
    # ~1 ~Tvar+Biotic ~1 ~1
    # VB.PecTB <- fm42.1
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
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60), Tvar=rep(0, length=60))
plot.col(inputmodel, inputdata)
mtext("Pecari tajacu (VB)", side=3, line=0, cex=0.7)




# Cercocebus sanjei (UDZ) nms[11]
    # ~1 ~ Biotic + Tmin ~1 ~1
    # UDZ.CerTB <- fm28.1
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
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60), Tmin=rep(0, length=60))
plot.col(inputmodel, inputdata)
mtext("Cercocebus sanjei (UDZ)", side=3, line=0, cex=0.7)


#mtext("Ecological dissimilarity of local community", side=1, line=0, outer=TRUE)
#mtext("Estimated Probability of Colonization", side=2, line=1, outer=TRUE)
#######End Colonization #########################



#### Begin Extinction ##############


# Mazama temama (VB_) nms[5]
    # ~1 ~1 ~Biotic ~1
    # VB.MazTB <- fm4.2
nmsK <- 5
try((fm4.2=colext(psiformula=~1,
                  gammaformula=~1,
                  epsilonformula=~Biotic,
                  pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)

if(exists("fm4.2")) {
  if(CondNum(fm4.2)<2000){
    if(CondNum(fm4.2)>0){mods=c(mods,fm4.2)}
} 
}}
inputmodel <- fm4.2
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60))
plot.ext(inputmodel, inputdata)
mtext("Mazama temama (VB)", side=3, line=0, cex=0.7)




# Cephalophus spadix (UDZ) nms[10]
     # ~1 ~1 ~Tvar + Biotic ~1
    # UDZ.CepTB <- fm42.2
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
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60), Tvar=rep(0, length=60))
plot.ext(inputmodel, inputdata)
mtext("Cephalophus spadix (UDZ)", side=3, line=0, cex=0.7)




# Cercopithecus lhoesti (BIF) nms[24]
   # ~1 ~Biotic ~Biotic ~1
   # BIF.CerTB <- fm4
nmsK <- 24
try((fm4=colext(psiformula=~1,
                gammaformula=~Biotic,
                epsilonformula=~Biotic,
                pformula=~1,data=umf,method="L-BFGS-B",control=list(maxit=20000))),silent=TRUE)
   
if(exists("fm4")) {
  if(CondNum(fm4)<2000){
    if(CondNum(fm4)>0){mods=c(mods,fm4)}
} 
}   }
inputmodel <- fm4
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60))
plot.ext(inputmodel, inputdata)
mtext("Cercopithecus lhoesti (BIF)", side=3, line=0, cex=0.7)



# Bdeogale crassicauda (UDZ) nms[8]
   # ~1 ~1 ~Biotic + Tmin ~1
   # UDZ.BdeTB <- fm28.2
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
inputdata <- data.frame(Biotic=seq(min(BIOTIC_166[[nmsK]]),max(BIOTIC_166[[nmsK]]),length=60), Tmin=rep(0, length=60))
plot.ext(inputmodel, inputdata)
mtext("Bdeogale crassicauda (UDZ)", side=3, line=0, cex=0.7)

#mtext("Ecological dissimilarity of local community", side=1, line=0, outer=TRUE)
#mtext("Estimated Probability of Extinction", side=2, line=1, outer=TRUE)
