# Plot partial relationships instead of the point estimate
# Show colonization probability as a function of elevation with confidence intervals

# Interested in estimates for 9 populations with significant range shifts

#### Colonization (fm2.1)
# Pan troglodytes (BIF) nms[27]
    # BIF.Pan <- fm2.1
# Eira barbara (YAN)  nms[43]
    # YAN.Eir <- fm2.1
# Tragulus kanchil (NAK) nms[54]
    # NAK.Tra <- fm2.1
# Mazama americana (YAN) nms[45]
    # YAN.Maz <- fm2.1

#### Extinction (fm2.2)
# Pecari tajacu (VB_) nms[6]
    # VB.Pec <- fm2.2
# Dasypus novemcinctus (VB_) nms[3]
    # VB.Das <- fm2.2
# Potamochoerus larvatus (UDZ) nms[19]
    # UDZ.Pot <- fm2.2
# Cricetomys gambianus (UDZ) nms[13]
    # UDZ.Cri <- fm2.2
# Cephalophus nigrifrons (BIF) nms[22]
    # BIF.Cep <- fm2.2

#E.ext <- predict(fm2.2, type="ext", newdata=ndm, appendData=TRUE)

set.panel(3,3)
par(mar=c(3,4,1,0))

# Pan troglodytes (BIF) nms[27] Colonization
    # BIF.Pan <- fm2.1
nwd <- data.frame(BIF.Elev)
pan <- predict(BIF.Pan, type="col", newdata=nwd, appendData=TRUE)
plot(BIF.Elev, pan$Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(1400,2400), ylim=c(0, 1))
points(BIF.Elev, pan$Predicted + 1.96*pan$SE, pch=19, col="gray55", cex=1.1)
points(BIF.Elev, ifelse(pan$Predicted - 1.96*pan$SE<0, 0, pan$Predicted - 1.96*pan$SE), pch=19, col="gray55", cex=1.1)
points(BIF.Elev, pan$Predicted, pch=19, cex=0.7)
mtext("Pan troglodytes (BIF)", side=3, line=0, cex=0.7)
#text(1700,0.4, "Pan troglodytes (BIF)", adj=0)


# Tragulus kanchil (NAK) nms[54]
    # NAK.Tra <- fm2.1
nwd <- data.frame(NAK.Elev)
tra <- predict(NAK.Tra, type="col", newdata=nwd, appendData=TRUE)
plot(NAK.Elev, tra$Predicted, pch=19, las=1, ylab="Colonization",
     xlab="", bty="n", xlim=c(300,1200), ylim=c(0, 1))
points(NAK.Elev, tra$Predicted + 1.96*tra$SE, pch=19, col="gray55", cex=1.1)
points(NAK.Elev, ifelse(tra$Predicted - 1.96*tra$SE<0, 0, tra$Predicted - 1.96*tra$SE), pch=19, col="gray55", cex=1.1)
points(NAK.Elev, tra$Predicted, pch=19, cex=0.7)
mtext("Tragulus kanchil (NAK)", side=3, line=0, cex=0.7)

# Cephalophus nigrifrons (BIF) nms[22] Extinction
     # BIF.Cep <- fm2.2
nwd <- data.frame(BIF.Elev)
cep <- predict(BIF.Cep, type="ext", newdata=nwd, appendData=TRUE)
plot(BIF.Elev, cep$Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(1400, 2400), ylim=c(0, 1))
points(BIF.Elev, ifelse(cep$Predicted + 1.96*cep$SE>1, 1, cep$Predicted + 1.96*cep$SE), pch=19, col="gray55", cex=1.1)
points(BIF.Elev, ifelse(cep$Predicted - 1.96*cep$SE<0, 0, cep$Predicted - 1.96*cep$SE), pch=19, col="gray55", cex=1.1)
points(BIF.Elev, cep$Predicted, pch=19, cex=0.7)
mtext("Cephalophus nigrifrons (BIF)", side=3, line=0, cex=0.7)
#text(1700,0.4, "Cephalophus nigrifrons (BIF)", adj=0)
#mtext(expression(hat(epsilon)), side=2, srt=180)

# Eira barbara (YAN)  nms[43]
    # YAN.Eir <- fm2.1
nwd <- data.frame(YAN.Elev)
eir <- predict(YAN.Eir, type="col", newdata=nwd, appendData=TRUE)
plot(YAN.Elev, eir$Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(300,1200), ylim=c(0, 1))
points(YAN.Elev, ifelse(eir$Predicted + 1.96*eir$SE>1, 1, eir$Predicted + 1.96*eir$SE), pch=19, col="gray55", cex=1.1)
points(YAN.Elev, ifelse(eir$Predicted - 1.96*eir$SE<0, 0, eir$Predicted - 1.96*eir$SE), pch=19, col="gray55", cex=1.1)
points(YAN.Elev, eir$Predicted, pch=19, cex=0.7)
mtext("Eira barbara (YAN)", side=3, line=0, cex=0.7)

# Mazama americana (YAN) nms[45]
    # YAN.Maz <- fm2.1
nwd <- data.frame(YAN.Elev)
maz <- predict(YAN.Maz, type="col", newdata=nwd, appendData=TRUE)
plot(YAN.Elev, maz$Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(300,1200), ylim=c(0, 1))
points(YAN.Elev, ifelse(maz$Predicted + 1.96*maz$SE>1, 1, maz$Predicted + 1.96*maz$SE), pch=19, col="gray55", cex=1.1)
points(YAN.Elev, ifelse(maz$Predicted - 1.96*maz$SE<0, 0, maz$Predicted - 1.96*maz$SE), pch=19, col="gray55", cex=1.1)
points(YAN.Elev, maz$Predicted, pch=19, cex=0.7)
mtext("Mazama americana (YAN)", side=3, line=0, cex=0.7)
#mtext(expression(hat(gamma)), side=2, srt=180)


# Cricetomys gambianus (UDZ) nms[13]
    # UDZ.Cri <- fm2.2
nwd <- data.frame(UDZ.Elev)
cri <- predict(UDZ.Cri, type="ext", newdata=nwd, appendData=TRUE)
plot(UDZ.Elev, cri$Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(400, 1800), ylim=c(0, 1))
points(UDZ.Elev, ifelse(cri$Predicted + 1.96*cri$SE>1, 1, cri$Predicted + 1.96*cri$SE), pch=19, col="gray55", cex=1.1)
points(UDZ.Elev, ifelse(cri$Predicted - 1.96*cri$SE<0, 0, cri$Predicted - 1.96*cri$SE), pch=19, col="gray55", cex=1.1)
points(UDZ.Elev, cri$Predicted, pch=19, cex=0.7)
mtext("Cricetomys gambianus (UDZ)", side=3, line=0, cex=0.7)



# Pecari tajacu (VB_) nms[6]
    # VB.Pec <- fm2.2
nwd <- data.frame(VB.Elev)
pec <- predict(VB.Pec, type="ext", newdata=nwd, appendData=TRUE)
plot(VB.Elev, pec$Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(0, 2600), ylim=c(0, 1))
points(VB.Elev, ifelse(pec$Predicted + 1.96*pec$SE>1, 1, pec$Predicted + 1.96*pec$SE), pch=19, col="gray55", cex=1.1)
points(VB.Elev, ifelse(pec$Predicted - 1.96*pec$SE<0, 0, pec$Predicted - 1.96*pec$SE), pch=19, col="gray55", cex=1.1)
points(VB.Elev, pec$Predicted, pch=19, cex=0.7)
mtext("Pecari tajacu (VB)", side=3, line=0, cex=0.7)


# Dasypus novemcinctus (VB_) nms[3]
    # VB.Das <- fm2.2
nwd <- data.frame(VB.Elev)
das <- predict(VB.Das, type="ext", newdata=nwd, appendData=TRUE)
plot(VB.Elev, das$Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(0, 2600), ylim=c(0, 1))
points(VB.Elev, ifelse(das$Predicted + 1.96*das$SE>1, 1, das$Predicted + 1.96*das$SE), pch=19, col="gray55", cex=1.1)
points(VB.Elev, ifelse(das$Predicted - 1.96*das$SE<0, 0, das$Predicted - 1.96*das$SE), pch=19, col="gray55", cex=1.1)
points(VB.Elev, das$Predicted, pch=19, cex=0.7)
mtext("Dasypus novemcinctus (VB)", side=3, line=0, cex=0.7)


# Potamochoerus larvatus (UDZ) nms[19]
    # UDZ.Pot <- fm2.2
nwd <- data.frame(UDZ.Elev)
pot <- predict(UDZ.Pot, type="ext", newdata=nwd, appendData=TRUE)
plot(UDZ.Elev, pot$Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(400, 1800), ylim=c(0, 1))
points(UDZ.Elev, ifelse(pot$Predicted + 1.96*pot$SE>1, 1, pot$Predicted + 1.96*pot$SE), pch=19, col="gray55", cex=1.1)
points(UDZ.Elev, ifelse(pot$Predicted - 1.96*pot$SE<0, 0, pot$Predicted - 1.96*pot$SE), pch=19, col="gray55", cex=1.1)
points(UDZ.Elev, pot$Predicted, pch=19, cex=0.7)
mtext("Potamochoerus larvatus (UDZ)", side=3, line=0, cex=0.7)



mtext("Elevation (m)", side=1, line=0, outer=TRUE)
mtext("Estimated Probability", side=2, line=1, outer=TRUE)






############### Code from Miguel Acevedo ##############
#This is some of the code I used to make the figure I showed you. It might help you.

ndm=data.frame(A=rep(a,length=50),Phorophyte=rep(0,length=50),Sfm=seq(min(Sfm),max(Sfm),length=50),Sm=seq(min(Sm),max(Sm),length=50))

#Colonization Rock

E.gam_asym=predict(modasym,type="col",newdata=ndm,appendData=TRUE)
E.gam_sym=predict(modsym,type="col",newdata=ndm,appendData=TRUE)

x=ndm$Sfm
y=E.gam_asym$Predicted
y.Err=1.96*E.gam_asym$SE
y.Up = y + y.Err; y.Dn =y-y.Err

x2=ndm$Sm
y2=E.gam_sym$Predicted 
y.Err2=1.96*E.gam_sym$SE 
y.Up2 = y2 + y.Err2; y.Dn2 =y2-y.Err2 

plot(x,y,type="l",lty=1,ylim=c(0,0.5),xlim=c(-1,2.5),col="black",ylab=expression(hat(gamma)),xlab="",cex.lab=1.0,main="",axes=FALSE)
axis(2,las=1)
axis(side = 1, at = x, labels = FALSE, tck = 0)

polygon(c(x,rev(x)),c(y.Up,rev(y.Dn)),col=rgb(64,64,64,25,maxColorValue=255),border=rgb(192,192,192,25,maxColorValue=255))
# Use the rev commands so the border moves logically around the shaded area
lines(x,y,type="l",lty=1,col="black") # put the means back on top of the polygon 


polygon(c(x2,rev(x2)),c(y.Up2,rev(y.Dn2)),col=rgb(192,192,192,25,maxColorValue=255),border=rgb(64,64,64,25,maxColorValue=255))
# Use the rev commands so the border moves logically around the shaded area
lines(x2,y2,type="l",lty=2,col="black") # put t
#text(2.15,0.099,"Rock phorophyte")
text(-0.05,0.48,"a) Rock phorophyte")
 legend("topright",c(expression(S[asym]),expression(S[sym])),lty=c(1,2),col=c("black","black"),bty='n')
