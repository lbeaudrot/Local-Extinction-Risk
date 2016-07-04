############## FIGURES ILLUSTRATING TEMPERATURE TRENDS ##############
# For Beaudrot, Acevedo, Lessard et al. (In prep)

library(denstrip)
library(reshape)
library(plyr)
library(fields)

############### SITE LEVEL DENSITY STRIPS #################
# Extract temperature data for all sites. Reduce PSH data to camera traps included in final analysis.

PSHdata <- eventsdata[eventsdata$Site.Code=="PSH",]
PSHdata <- data.frame(PSHdata$Sampling.Unit.Name, PSHdata$temp.degreesC)
PSHdata <- PSHdata[1:75936,] # Reduces dataset to only array 1
PSHdata <- PSHdata[-c(38052:43974),] # Removes camera trap 21 which has faulty temperature readings

place <- list(
           eventsdata$temp.degreesC[eventsdata$Site.Code=="RNF"],
           eventsdata$temp.degreesC[eventsdata$Site.Code=="NAK"],
           eventsdata$temp.degreesC[eventsdata$Site.Code=="UDZ"],
           eventsdata$temp.degreesC[eventsdata$Site.Code=="BIF"],
           eventsdata$temp.degreesC[eventsdata$Site.Code=="VB-"],
           eventsdata$temp.degreesC[eventsdata$Site.Code=="YAN"],
           PSHdata$PSHdata.temp.degreesC)

# Create graph with density strips for all 7 TEAM sites
y <- 5:40
par(mar=c(3,3,1,1))
plot(y, type="n", 
     xlim=c(0.75, 7.75), 
     bty="n", las=1, ylab="Temperature (C)", 
     cex.axis=1.5, bg="transparent", mgp=c(3,1,0), xaxt="n")
axis(1, at=c(1,2,3,4,5,6,7), 
     labels=c("RNF", "NAK", "UDZ", "BIF", "VB", "YAN", "PSH"), 
     cex.axis=1.5)
#mtext(text="Temperature (C)", side=2, line=2, cex=1.25)

for(i in 1:length(place)){
      denstrip(place[[i]], na.rm=TRUE, at=i, horiz=FALSE, width=0.5,
         ticks=c(min(place[[i]], na.rm=TRUE), max(place[[i]], na.rm=TRUE)), 
         mticks=median(place[[i]], na.rm=TRUE))
} 

# Plot only one site per graph (to add to map of TEAM sites)
y <- 5:40
par(mar=c(0,2.5,0,0))
plot(y, type="n", 
    xlim=c(1,1),
     bty="n", las=0, ylab="Temperature (C)", 
     cex.axis=1.9, bg="transparent", mgp=c(3,1,0), xaxt="n")

# Change k to number in place that corresponds to site
k=7
for(i in 1){
      denstrip(place[[k]], na.rm=TRUE, at=1, horiz=FALSE, width=0.5,
         ticks=c(min(place[[k]], na.rm=TRUE), max(place[[k]], na.rm=TRUE)), 
         mticks=median(place[[k]], na.rm=TRUE))
} 



################### DENSITY PLOTS ##########################
############ Create multi-panel figure with density plots for each site and variable ##########
# For each TEAM site and for each temperature variable (Tmin, Tmax, Tvar, etc.)
# Plots camera trap specific trends 

library(fields)
set.panel(7,4)
par(mar=c(2,3.5,2,1), oma=c(3,3,3,1))

###### BIF

ys <- as.numeric(colnames(BIF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

plot(density(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="Maximum", ylab="", xlab="", col="gray55", mgp=c(2,1,0), ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")
mtext("BIF", side=2, line=2.5, outer=FALSE, cex=1.2)

plot(density(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="Minimum", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="Mean", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="Variance", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

###### NAK

ys <- as.numeric(colnames(NAK.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

plot(density(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", mgp=c(2,1,0), ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")
mtext("NAK", side=2, line=2.5, outer=FALSE, cex=1.2)

plot(density(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


######## PSH

ys <- as.numeric(colnames(PSH.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

PSH.Tmax <- PSH.Tmax[1:30,]
PSH.Tmin <- PSH.Tmin[1:30,]
PSH.Mean <- PSH.Mean[1:30,]
PSH.Tvar <- PSH.Tvar[1:30,]

plot(density(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", mgp=c(2,1,0), ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")
mtext("PSH", side=2, line=2.5, outer=FALSE, cex=1.2)

plot(density(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")



########## RNF

ys <- as.numeric(colnames(RNF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

plot(density(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", mgp=c(2,1,0), ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")
mtext("RNF", side=2, line=2.5, outer=FALSE, cex=1.2)

plot(density(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


############## UDZ

ys <- as.numeric(colnames(UDZ.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

plot(density(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", mgp=c(2,1,0), ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")
mtext("UDZ", side=2, line=2.5, outer=FALSE, cex=1.2)

plot(density(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")



################# VB


ys <- as.numeric(colnames(VB.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

plot(density(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", mgp=c(2,1,0), ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")
mtext("VB", side=2, line=2.5, outer=FALSE, cex=1.2)

plot(density(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


###################### YAN

ys <- as.numeric(colnames(YAN.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

plot(density(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", mgp=c(2,1,0), ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")
mtext("YAN", side=2, line=2.5, outer=FALSE, cex=1.2)

plot(density(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

plot(density(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), lwd=3, las=1, 
     main="", ylab="", xlab="",  col="gray55", ylim=c(0,1.7), xlim=c(-8,10))
abline(v=median(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

###### add text labels to outer margins

mtext("Density of camera traps", side=2, line=1, outer=TRUE)
mtext("Temperature trend", side=1, line=1, outer=TRUE)




########################## EARLIER VERSIONS #################################
# Histograms (instead of density plots)
# Temperature trends as a function of elevation

############ Create multi-panel figure with histograms for each site and variable ##########
library(fields)
set.panel(7,4)
par(mar=c(2,3.5,2,1), oma=c(3,3,3,1))

###### BIF

ys <- as.numeric(colnames(BIF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="Maximum", ylab="BIF", xlab="", breaks=20, col="gray55", mgp=c(2,1,0), ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="Minimum", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="Mean", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="Variance", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

###### NAK

ys <- as.numeric(colnames(NAK.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="NAK", xlab="", breaks=20, col="gray55", mgp=c(2,1,0), ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


######## PSH

ys <- as.numeric(colnames(PSH.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

PSH.Tmax <- PSH.Tmax[1:30,]
PSH.Tmin <- PSH.Tmin[1:30,]
PSH.Mean <- PSH.Mean[1:30,]
PSH.Tvar <- PSH.Tvar[1:30,]

hist(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="PSH", xlab="", breaks=20, col="gray55", mgp=c(2,1,0), ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")



########## RNF

ys <- as.numeric(colnames(RNF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="RNF", xlab="", breaks=20, col="gray55", mgp=c(2,1,0), ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


############## UDZ

ys <- as.numeric(colnames(UDZ.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="UDZ", xlab="", breaks=20, col="gray55", mgp=c(2,1,0), ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")



################# VB


ys <- as.numeric(colnames(VB.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="VB", xlab="", breaks=20, col="gray55", mgp=c(2,1,0), ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


###################### YAN

ys <- as.numeric(colnames(YAN.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="YAN", xlab="", breaks=20, col="gray55", mgp=c(2,1,0), ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55", ylim=c(0,20), xlim=c(-6,6))
abline(v=median(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

###### add text labels to outer margins

mtext("Number of camera traps", side=2, line=1, outer=TRUE)
mtext("Camera trap temperature trend (slope of linear regression over time)", side=1, line=1, outer=TRUE)



############# TEMPERATURE ~ ELEVATION ###############

# Plot trends as a function of elevation
set.panel(7,4)
par(mar=c(2,4,2,1), oma=c(3,3,3,1))

ys <- as.numeric(colnames(BIF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="BIF", las=1, main="Maximum", mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)
median(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="", las=1, main="Minimum")
abline(coef(lm(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)
median(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="", las=1, main="Mean")
abline(coef(lm(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)
median(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="", las=1, main="Variance")
abline(coef(lm(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)
median(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

ys <- as.numeric(colnames(NAK.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="NAK", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)
median(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)
median(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)
median(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)
median(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE)


#plot(apply(as.matrix(PSH.Tvar[-21]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21])
PSH.Elev <- PSH.Elev[1:30]
ys <- as.numeric(colnames(PSH.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21], xlab="", ylab="PSH", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[1], coef(lm(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[2], lty=2)
median(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21], xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[1], coef(lm(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[2], lty=2)
median(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21], xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[1], coef(lm(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[2], lty=2)
median(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21], xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[1], coef(lm(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21]))[2], lty=2)
median(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

ys <- as.numeric(colnames(RNF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="RNF", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)
median(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)
median(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)
median(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)
median(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

ys <- as.numeric(colnames(UDZ.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="UDZ", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)
median(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)
median(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)
median(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)
median(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

ys <- as.numeric(colnames(VB.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="VB", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)
median(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)
median(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)
median(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)
median(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE)


ys <- as.numeric(colnames(YAN.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="YAN", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)
median(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)
median(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)
median(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

plot(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)
median(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE)

mtext("Camera trap temperature trend (slope of linear regression over time)", side=2, line=1, outer=TRUE)
mtext("Elevation (meters)", side=1, line=1, outer=TRUE)

