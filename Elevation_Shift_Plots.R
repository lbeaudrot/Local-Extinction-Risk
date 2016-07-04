# Plot partial relationships from colext models for 9 populations with significant range shifts
# Show colonization probability as a function of elevation with confidence intervals

#### DIRECTIONS #######
# To create plot for each species within the overall plot matrix: 
#   1) Run "Unmarked Colext Analysis Elevation Only.R" for nms[k] in the loop specified below for species i
#   2) Run subset of plotting code below for species i
#   3) Repeat until all species have been modeled and plotted

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


############### BEGIN PLOTTING #######################
library(fields)
set.panel(3,3)
par(mar=c(3,4,1,0))

# Pan troglodytes (BIF) nms[27] Colonization
    # BIF.Pan <- fm2.1
pan.nwd <- data.frame(seq(min(BIF.Elev),max(BIF.Elev),length=60))
pan <- predict(BIF.Pan, type="col", newdata=pan.nwd, appendData=TRUE)
pan.sort <- data.frame(BIF.Elev, pan$Predicted)
pan.sort <- pan.sort[do.call(order, pan.sort),]
plot(pan.sort$BIF.Elev, pan.sort$pan.Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(1400,2400), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(sort(BIF.Elev),rev(sort(BIF.Elev))),c(sort(pan$Predicted + 1.96*pan$SE),
                          rev(sort(ifelse(pan$Predicted - 1.96*pan$SE<0, 0, pan$Predicted - 1.96*pan$SE)))), 
        col=rgb(0,100,255,75,maxColorValue=255))
mtext("Pan troglodytes (BIF)", side=3, line=0, cex=0.7)


# Tragulus kanchil (NAK) nms[54]
    # NAK.Tra <- fm2.1
tra.nwd <- data.frame(seq(min(NAK.Elev),max(NAK.Elev),length=60))
tra <- predict(NAK.Tra, type="col", newdata=data.frame(tra.nwd), appendData=TRUE)
tra.sort <- data.frame(NAK.Elev, tra$Predicted)
tra.sort <- tra.sort[do.call(order, tra.sort),]
plot(tra.sort$NAK.Elev, tra.sort$tra.Predicted, pch=19, las=1, ylab="Colonization",
     xlab="", bty="n", xlim=c(300,1200), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(NAK.Elev)),sort(NAK.Elev)),c(sort(tra$Predicted + 1.96*tra$SE)
                          ,rev(sort(ifelse(tra$Predicted - 1.96*tra$SE<0, 0, tra$Predicted - 1.96*tra$SE)))), 
        col=rgb(0,100,255,75,maxColorValue=255))
mtext("Tragulus kanchil (NAK)", side=3, line=0, cex=0.7)

# Cephalophus nigrifrons (BIF) nms[22] Extinction
     # BIF.Cep <- fm2.2
cep.nwd <- data.frame(seq(min(BIF.Elev),max(BIF.Elev),length=60))
cep <- predict(BIF.Cep, type="ext", newdata=cep.nwd, appendData=TRUE)
cep.sort <- data.frame(BIF.Elev, cep$Predicted)
cep.sort <- cep.sort[do.call(order, cep.sort),]
plot(cep.sort$BIF.Elev, cep.sort$cep.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(1400, 2400), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(BIF.Elev)),sort(BIF.Elev)),c(sort(ifelse(cep$Predicted + 1.96*cep$SE>1, 1, cep$Predicted + 1.96*cep$SE))
                          ,rev(sort(ifelse(cep$Predicted - 1.96*cep$SE<0, 0, cep$Predicted - 1.96*cep$SE)))), 
        col=rgb(225,95,0,125,maxColorValue=255))
mtext("Cephalophus nigrifrons (BIF)", side=3, line=0, cex=0.7)
#mtext(expression(hat(epsilon)), side=2, srt=180)

# Eira barbara (YAN)  nms[43]
    # YAN.Eir <- fm2.1
eir.nwd <- data.frame(seq(min(YAN.Elev),max(YAN.Elev),length=60))
eir <- predict(YAN.Eir, type="col", newdata=eir.nwd, appendData=TRUE)
eir.sort <- data.frame(YAN.Elev, eir$Predicted)
eir.sort <- eir.sort[do.call(order, eir.sort),]
plot(eir.sort$YAN.Elev, eir.sort$eir.Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(300,1200), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(sort(YAN.Elev),rev(sort(YAN.Elev))),c(sort(ifelse(eir$Predicted + 1.96*eir$SE>1, 1, eir$Predicted + 1.96*eir$SE))
                          ,rev(sort(ifelse(eir$Predicted - 1.96*eir$SE<0, 0, eir$Predicted - 1.96*eir$SE)))), 
        col=rgb(0,100,255,75,maxColorValue=255))
mtext("Eira barbara (YAN)", side=3, line=0, cex=0.7)

# Mazama americana (YAN) nms[45]
    # YAN.Maz <- fm2.1
maz.nwd <- data.frame(YAN.Elev)
maz <- predict(YAN.Maz, type="col", newdata=maz.nwd, appendData=TRUE)
maz.sort <- data.frame(YAN.Elev, maz$Predicted)
maz.sort <- maz.sort[do.call(order, maz.sort),]
plot(maz.sort$YAN.Elev, maz.sort$maz.Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(300,1200), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(YAN.Elev)),(sort(YAN.Elev))),c(sort(ifelse(maz$Predicted + 1.96*maz$SE>1, 1, maz$Predicted + 1.96*maz$SE))
                          ,rev(sort(ifelse(maz$Predicted - 1.96*maz$SE<0, 0, maz$Predicted - 1.96*maz$SE)))), 
        col=rgb(0,100,255,75,maxColorValue=255))
mtext("Mazama americana (YAN)", side=3, line=0, cex=0.7)
#mtext(expression(hat(gamma)), side=2, srt=180)


# Cricetomys gambianus (UDZ) nms[13]
    # UDZ.Cri <- fm2.2
cri.nwd <- data.frame(UDZ.Elev)
cri <- predict(UDZ.Cri, type="ext", newdata=cri.nwd, appendData=TRUE)
cri.sort <- data.frame(UDZ.Elev, cri$Predicted)
cri.sort <- cri.sort[do.call(order, cri.sort),]
plot(cri.sort$UDZ.Elev, cri.sort$cri.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(400, 1800), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(UDZ.Elev)),(sort(UDZ.Elev))),c(sort(ifelse(cri$Predicted + 1.96*cri$SE>1, 1, cri$Predicted + 1.96*cri$SE))
                          ,rev(sort(ifelse(cri$Predicted - 1.96*cri$SE<0, 0, cri$Predicted - 1.96*cri$SE)))), 
        col=rgb(225,95,0,125,maxColorValue=255))
mtext("Cricetomys gambianus (UDZ)", side=3, line=0, cex=0.7)


# Potamochoerus larvatus (UDZ) nms[19]
    # UDZ.Pot <- fm2.2
pot.nwd <- data.frame(UDZ.Elev)
pot <- predict(UDZ.Pot, type="ext", newdata=pot.nwd, appendData=TRUE)
pot.sort <- data.frame(UDZ.Elev, pot$Predicted)
pot.sort <- pot.sort[do.call(order, pot.sort),]
plot(pot.sort$UDZ.Elev, pot.sort$pot.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(400, 1800), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c((sort(UDZ.Elev)),rev(sort(UDZ.Elev))),c(sort(ifelse(pot$Predicted + 1.96*pot$SE>1, 1, pot$Predicted + 1.96*pot$SE))
                          ,rev(sort(ifelse(pot$Predicted - 1.96*pot$SE<0, 0, pot$Predicted - 1.96*pot$SE)))), 
        col=rgb(225,95,0,125,maxColorValue=255))
mtext("Potamochoerus larvatus (UDZ)", side=3, line=0, cex=0.7)

# Dasypus novemcinctus (VB_) nms[3]
    # VB.Das <- fm2.2
das.nwd <- data.frame(VB.Elev)
das <- predict(VB.Das, type="ext", newdata=das.nwd, appendData=TRUE)
das.sort <- data.frame(VB.Elev, das$Predicted)
das.sort <- das.sort[do.call(order, das.sort),]
plot(das.sort$VB.Elev, das.sort$das.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(0, 2600), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c((sort(VB.Elev)),rev(sort(VB.Elev))),c(sort(ifelse(das$Predicted + 1.96*das$SE>1, 1, das$Predicted + 1.96*das$SE))
                          ,rev(sort(ifelse(das$Predicted - 1.96*das$SE<0, 0, das$Predicted - 1.96*das$SE)))), 
        col=rgb(225,95,0,125,maxColorValue=255))
mtext("Dasypus novemcinctus (VB)", side=3, line=0, cex=0.7)


# Pecari tajacu (VB_) nms[6]
    # VB.Pec <- fm2.2
pec.nwd <- data.frame(VB.Elev)
pec <- predict(VB.Pec, type="ext", newdata=pec.nwd, appendData=TRUE)
pec.sort <- data.frame(VB.Elev, pec$Predicted)
pec.sort <- pec.sort[do.call(order, pec.sort),]
plot(pec.sort$VB.Elev, pec.sort$pec.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(0, 2600), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c((sort(VB.Elev)),rev(sort(VB.Elev))),c(sort(ifelse(pec$Predicted + 1.96*pec$SE>1, 1, pec$Predicted + 1.96*pec$SE))
                          ,rev(sort(ifelse(pec$Predicted - 1.96*pec$SE<0, 0, pec$Predicted - 1.96*pec$SE)))), 
        col=rgb(225,95,0,125,maxColorValue=255))
mtext("Pecari tajacu (VB)", side=3, line=0, cex=0.7)

# Outer margin text
mtext("Elevation (m)", side=1, line=0, outer=TRUE)
mtext("Estimated Probability", side=2, line=1, outer=TRUE)

############## END PLOTTING #########################




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
