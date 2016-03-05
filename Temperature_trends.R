# For each TEAM site and for each temperature variable (Tmin, Tmax, Tvar, etc.)
# Plot histogram of camera trap specific trends 
library(fields)

ys <- as.numeric(colnames(VB.Tvar))
ys <- as.numeric(colnames(UDZ.Tvar))
ys <- as.numeric(colnames(BIF.Tvar))
ys <- as.numeric(colnames(PSH.Tvar))
ys <- as.numeric(colnames(NAK.Tvar))
ys <- as.numeric(colnames(RNF.Tvar))
ys <- as.numeric(colnames(YAN.Tvar))


mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
#apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn)



############ Create multi-panel figure with histograms for each site and variable ##########

set.panel(7,4)
par(mar=c(2,3.5,2,1), oma=c(3,3,3,1))

###### BIF

ys <- as.numeric(colnames(BIF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="Maximum", ylab="BIF", xlab="", breaks=20, col="gray55", mgp=c(2,1,0))
abline(v=median(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="Minimum", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="Mean", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="Variance", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

###### NAK

ys <- as.numeric(colnames(NAK.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="NAK", xlab="", breaks=20, col="gray55", mgp=c(2,1,0))
abline(v=median(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


######## PSH

ys <- as.numeric(colnames(PSH.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

PSH.Tmax <- PSH.Tmax[1:30,]
PSH.Tmin <- PSH.Tmin[1:30,]
PSH.Mean <- PSH.Mean[1:30,]
PSH.Tvar <- PSH.Tvar[1:30,]

hist(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="PSH", xlab="", breaks=20, col="gray55", mgp=c(2,1,0))
abline(v=median(apply(as.matrix(PSH.Tmax[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(PSH.Tmin[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(PSH.Mean[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(PSH.Tvar[-21,]), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")



########## RNF

ys <- as.numeric(colnames(RNF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="RNF", xlab="", breaks=20, col="gray55", mgp=c(2,1,0))
abline(v=median(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


############## UDZ

ys <- as.numeric(colnames(UDZ.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="UDZ", xlab="", breaks=20, col="gray55", mgp=c(2,1,0))
abline(v=median(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")



################# VB


ys <- as.numeric(colnames(VB.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="VB", xlab="", breaks=20, col="gray55", mgp=c(2,1,0))
abline(v=median(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")


###################### YAN

ys <- as.numeric(colnames(YAN.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]

hist(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="YAN", xlab="", breaks=20, col="gray55", mgp=c(2,1,0))
abline(v=median(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

hist(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn), las=1, 
     main="", ylab="", xlab="", breaks=20, col="gray55")
abline(v=median(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn), na.rm=TRUE), col="red")

###### add text labels to outer margins

mtext("Number of camera traps", side=2, line=1, outer=TRUE)
mtext("Camera trap temperature trend (slope of linear regression over time)", side=1, line=1, outer=TRUE)








# Plot trends as a function of elevation
set.panel(6,4)
par(mar=c(2,4,2,1), oma=c(3,3,3,1))

ys <- as.numeric(colnames(BIF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="BIF", las=1, main="Maximum", mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Tmax), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)

plot(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="", las=1, main="Minimum")
abline(coef(lm(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Tmin), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)

plot(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="", las=1, main="Mean")
abline(coef(lm(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Mean), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)

plot(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn) ~ BIF.Elev, xlab="", ylab="", las=1, main="Variance")
abline(coef(lm(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[1], coef(lm(apply(as.matrix(BIF.Tvar), MARGIN=1, FUN=mod.fn) ~ BIF.Elev))[2], lty=2)

ys <- as.numeric(colnames(NAK.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="NAK", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Tmax), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)

plot(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Tmin), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)

plot(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Mean), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)

plot(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn) ~ NAK.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[1], coef(lm(apply(as.matrix(NAK.Tvar), MARGIN=1, FUN=mod.fn) ~ NAK.Elev))[2], lty=2)


#plot(apply(as.matrix(PSH.Tvar[-21]), MARGIN=1, FUN=mod.fn) ~ PSH.Elev[-21])

ys <- as.numeric(colnames(RNF.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="RNF", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Tmax), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)

plot(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Tmin), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)

plot(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Mean), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)

plot(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn) ~ RNF.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[1], coef(lm(apply(as.matrix(RNF.Tvar), MARGIN=1, FUN=mod.fn) ~ RNF.Elev))[2], lty=2)

ys <- as.numeric(colnames(UDZ.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="UDZ", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Tmax), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)

plot(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Tmin), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)

plot(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Mean), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)

plot(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[1], coef(lm(apply(as.matrix(UDZ.Tvar), MARGIN=1, FUN=mod.fn) ~ UDZ.Elev))[2], lty=2)

ys <- as.numeric(colnames(VB.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="VB", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Tmax), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)

plot(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Tmin), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)

plot(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Mean), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)

plot(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn) ~ VB.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[1], coef(lm(apply(as.matrix(VB.Tvar), MARGIN=1, FUN=mod.fn) ~ VB.Elev))[2], lty=2)


ys <- as.numeric(colnames(YAN.Tvar))
mod.fn <- function(x) coef(lm(x ~ ys), na.omit=TRUE)[2]
plot(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="YAN", las=1, mgp=c(2,1,0))
abline(coef(lm(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Tmax), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)

plot(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Tmin), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)

plot(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Mean), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)

plot(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn) ~ YAN.Elev, xlab="", ylab="", las=1)
abline(coef(lm(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[1], coef(lm(apply(as.matrix(YAN.Tvar), MARGIN=1, FUN=mod.fn) ~ YAN.Elev))[2], lty=2)

mtext("Camera trap temperature trend (slope of linear regression over time)", side=2, line=1, outer=TRUE)
mtext("Elevation (meters)", side=1, line=1, outer=TRUE)


