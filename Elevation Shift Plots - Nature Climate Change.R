### COLONIZATION (i.e. Elevation Expansion) ########
# Top Model: ~1 ~Elevation ~1 ~1 (fm1.1)
# Tragulus kanchil (NAK) nms[54]
# Sus scrofa (NAK) nms[53]

# Top Model: ~Elevation ~Elevation ~1 ~1 (fm2.1)
# Nasua nasua (YAN) nms[46]

### EXTINCTION (i.e. Elevation Retraction) ##########
# Top Model: ~1 ~1 ~Elevation ~1 (fm1.2)
# Dasypus novemcinctus (VB_) nms[3]
# Pecari tajacu (VB_) nms[6]
# Cricetomys gambianus (UDZ) nms[13]
# Potamochoerus larvatus (UDZ) nms[19]
# Cephalophus nigrifrons (BIF) nms[22]

# Top Model: ~Elevation ~1 ~Elevation ~1 (fm2.2)
# Pan troglodytes (BIF) nms[27]
# Fossa fossana (RNF) nms[56]

par(mar=c(2, 3, 1, 1))


# Tragulus kanchil (NAK) nms[54]
NAK.tra <- fm1.1
tra.nwd <- data.frame(seq(min(NAK.Elev),max(NAK.Elev),length=60))
tra <- predict(NAK.tra, type="col", newdata=tra.nwd, appendData=TRUE)
tra.sort <- data.frame(NAK.Elev, tra$Predicted)
tra.sort <- tra.sort[do.call(order, tra.sort),]
plot(tra.sort$NAK.Elev, tra.sort$tra.Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(350,1200), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(NAK.Elev)), sort(NAK.Elev)),c(sort(ifelse(tra$Predicted + 1.96*tra$SE > 1, 1, tra$Predicted + 1.96*tra$SE)),
                                                rev(sort(ifelse(tra$Predicted - 1.96*tra$SE<0, 0, tra$Predicted - 1.96*tra$SE)))), 
                                                col=rgb(0,100,255,75,maxColorValue=255))


# Sus scrofa (NAK) nms[53]
NAK.sus <- fm1.1
sus.nwd <- data.frame(seq(min(NAK.Elev),max(NAK.Elev),length=60))
sus <- predict(NAK.sus, type="col", newdata=sus.nwd, appendData=TRUE)
sus.sort <- data.frame(NAK.Elev, sus$Predicted)
sus.sort <- sus.sort[do.call(order, sus.sort),]
plot(sus.sort$NAK.Elev, sus.sort$sus.Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(350,1200), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(NAK.Elev)), sort(NAK.Elev)),c(sort(ifelse(sus$Predicted + 1.96*sus$SE > 1, 1, sus$Predicted + 1.96*sus$SE)),
                                                rev(sort(ifelse(sus$Predicted - 1.96*sus$SE<0, 0, sus$Predicted - 1.96*sus$SE)))), 
                                                col=rgb(0,100,255,75,maxColorValue=255))

# Nasua nasua (YAN) nms[46]
YAN.nas <- fm2.1
nas.nwd <- data.frame(seq(min(YAN.Elev),max(YAN.Elev),length=61))
nas <- predict(YAN.nas, type="col", newdata=nas.nwd, appendData=TRUE)
nas.sort <- data.frame(YAN.Elev, nas$Predicted)
nas.sort <- nas.sort[do.call(order, nas.sort),]
plot(nas.sort$YAN.Elev, nas.sort$nas.Predicted, pch=19, las=1, ylab="Colonization", 
     xlab="", bty="n", xlim=c(400,1200), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(sort(YAN.Elev), rev(sort(YAN.Elev))),c(sort(ifelse(nas$Predicted + 1.96*nas$SE > 1, 1, nas$Predicted + 1.96*nas$SE)),
                                                rev(sort(ifelse(nas$Predicted - 1.96*nas$SE<0, 0, nas$Predicted - 1.96*nas$SE)))), 
                                                col=rgb(0,100,255,75,maxColorValue=255))

# Dasypus novemcinctus (VB_) nms[3]
VB.das <- fm1.2
das.nwd <- data.frame(seq(min(VB.Elev),max(VB.Elev),length=60))
das <- predict(VB.das, type="ext", newdata=das.nwd, appendData=TRUE)
das.sort <- data.frame(VB.Elev, das$Predicted)
das.sort <- das.sort[do.call(order, das.sort),]
plot(das.sort$VB.Elev, das.sort$das.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(0,2600), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(sort(VB.Elev), rev(sort(VB.Elev))),c(sort(ifelse(das$Predicted + 1.96*das$SE > 1, 1, das$Predicted + 1.96*das$SE)),
                                                rev(sort(ifelse(das$Predicted - 1.96*das$SE<0, 0, das$Predicted - 1.96*das$SE)))), 
                                                col=rgb(225,95,0,125,maxColorValue=255))


# Pecari tajacu (VB_) nms[6]
VB.pec <- fm1.2
pec.nwd <- data.frame(seq(min(VB.Elev),max(VB.Elev),length=60))
pec <- predict(VB.pec, type="ext", newdata=pec.nwd, appendData=TRUE)
pec.sort <- data.frame(VB.Elev, pec$Predicted)
pec.sort <- pec.sort[do.call(order, pec.sort),]
plot(pec.sort$VB.Elev, pec.sort$pec.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(0,2600), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(sort(VB.Elev), rev(sort(VB.Elev))),c(sort(ifelse(pec$Predicted + 1.96*pec$SE > 1, 1, pec$Predicted + 1.96*pec$SE)),
                                                rev(sort(ifelse(pec$Predicted - 1.96*pec$SE<0, 0, pec$Predicted - 1.96*pec$SE)))), 
                                                col=rgb(225,95,0,125,maxColorValue=255))


# Cricetomys gambianus (UDZ) nms[13]
UDZ.cri <- fm1.2
cri.nwd <- data.frame(seq(min(UDZ.Elev),max(UDZ.Elev),length=61))
cri <- predict(UDZ.cri, type="ext", newdata=cri.nwd, appendData=TRUE)
cri.sort <- data.frame(UDZ.Elev, cri$Predicted)
cri.sort <- cri.sort[do.call(order, cri.sort),]
plot(cri.sort$UDZ.Elev, cri.sort$cri.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(400,1800), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(UDZ.Elev)), sort(UDZ.Elev)),c(sort(ifelse(cri$Predicted + 1.96*cri$SE > 1, 1, cri$Predicted + 1.96*cri$SE)),
                                                rev(sort(ifelse(cri$Predicted - 1.96*cri$SE<0, 0, cri$Predicted - 1.96*cri$SE)))), 
                                                col=rgb(225,95,0,125,maxColorValue=255))

# Potamochoerus larvatus (UDZ) nms[19]
UDZ.pot <- fm1.2
pot.nwd <- data.frame(seq(min(UDZ.Elev),max(UDZ.Elev),length=61))
pot <- predict(UDZ.pot, type="ext", newdata=pot.nwd, appendData=TRUE)
pot.sort <- data.frame(UDZ.Elev, pot$Predicted)
pot.sort <- pot.sort[do.call(order, pot.sort),]
plot(pot.sort$UDZ.Elev, pot.sort$pot.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(400,1800), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(sort(UDZ.Elev), rev(sort(UDZ.Elev))),c(sort(ifelse(pot$Predicted + 1.96*pot$SE > 1, 1, pot$Predicted + 1.96*pot$SE)),
                                                rev(sort(ifelse(pot$Predicted - 1.96*pot$SE<0, 0, pot$Predicted - 1.96*pot$SE)))), 
                                                col=rgb(225,95,0,125,maxColorValue=255))

# Cephalophus nigrifrons (BIF) nms[22]
BIF.cep <- fm1.2
cep.nwd <- data.frame(seq(min(BIF.Elev),max(BIF.Elev),length=60))
cep <- predict(BIF.cep, type="ext", newdata=cep.nwd, appendData=TRUE)
cep.sort <- data.frame(BIF.Elev, cep$Predicted)
cep.sort <- cep.sort[do.call(order, cep.sort),]
plot(cep.sort$BIF.Elev, cep.sort$cep.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(1400,2400), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(BIF.Elev)),sort(BIF.Elev)),c(sort(ifelse(cep$Predicted + 1.96*cep$SE > 1, 1, cep$Predicted + 1.96*cep$SE)),
                                                rev(sort(ifelse(cep$Predicted - 1.96*cep$SE<0, 0, cep$Predicted - 1.96*cep$SE)))), 
                                                col=rgb(225,95,0,125,maxColorValue=255))

# Pan troglodytes (BIF) nms[27]
BIF.Pan <- fm2.2
pan.nwd <- data.frame(seq(min(BIF.Elev),max(BIF.Elev),length=60))
pan <- predict(BIF.Pan, type="ext", newdata=pan.nwd, appendData=TRUE)
pan.sort <- data.frame(BIF.Elev, pan$Predicted)
pan.sort <- pan.sort[do.call(order, pan.sort),]
plot(pan.sort$BIF.Elev, pan.sort$pan.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(1400,2400), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(BIF.Elev)),sort(BIF.Elev)),c(sort(ifelse(pan$Predicted + 1.96*pan$SE > 1, 1, pan$Predicted + 1.96*pan$SE)),
                                                rev(sort(ifelse(pan$Predicted - 1.96*pan$SE<0, 0, pan$Predicted - 1.96*pan$SE)))), 
                                                col=rgb(225,95,0,125,maxColorValue=255))
#pdf(file="BIF_Pan_troglodytes_fm2.2.pdf", width=4, height=4)
    
# Fossa fossana (RNF) nms[56]
RNF.Fos <- fm2.2
fos.nwd <- data.frame(seq(min(RNF.Elev),max(RNF.Elev),length=60))
fos <- predict(RNF.Fos, type="ext", newdata=fos.nwd, appendData=TRUE)
fos.sort <- data.frame(RNF.Elev, fos$Predicted)
fos.sort <- fos.sort[do.call(order, fos.sort),]
plot(fos.sort$RNF.Elev, fos.sort$fos.Predicted, pch=19, las=1, ylab="Extinction", 
     xlab="", bty="n", xlim=c(700,1300), ylim=c(0, 1), cex=0.8, type="line", lwd=2)
polygon(c(rev(sort(RNF.Elev)),sort(RNF.Elev)),c(sort(ifelse(fos$Predicted + 1.96*fos$SE > 1, 1, fos$Predicted + 1.96*fos$SE)),
                                                rev(sort(ifelse(fos$Predicted - 1.96*fos$SE<0, 0, fos$Predicted - 1.96*fos$SE)))), 
                                                col=rgb(225,95,0,125,maxColorValue=255))
