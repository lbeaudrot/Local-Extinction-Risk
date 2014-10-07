library(lubridate)
library(reshape)

#ctdata <- f.teamdb.query("camera trap")
#load(ctdata)
#alldata <- ctdata

#alldata<-f.fix.data2(alldata)
#Site.Code <- substr(alldata$Sampling.Unit.Name,4,6)
#alldata <- cbind(alldata, Site.Code)

# Set threshold to one day (24 hours * 60 minutes = 1440 minutes) and the number of replicates for each camera trap as 30 days
# Note that f.separate.events takes a good bit of time to run on complete CT dataset. Avoid replicating where possible.
#alldata <- f.order.data(alldata)
#eventsdata <- f.separate.events(alldata, thresh=(1440))
#allevents <- f.separate.events(alldata, thresh=(1))


CT.dates <- unique(eventsdata$Photo.Date[eventsdata$Site.Code=="VB-"])
CT.names <- unique(eventsdata$Sampling.Unit.Name[eventsdata$Site.Code=="VB-"])

empty <- matrix(999, length(CT.names), length(CT.dates))
rownames(empty) <- CT.names
colnames(empty) <- as.character(as.Date(CT.dates))

# Code from EstimateSpeciesRichness.R to determine # sampling days per CT
# Control for the number of sampling day per camera trap
#effortdays <- as.matrix(cbind(data.use$Sampling.Unit.Name, f.start.minus.end(data.use)))
#effortdays <- as.matrix(unique(effortdays))
#effortdays <- effortdays[match(colnames(X), effortdays[,1]),]
#edays <- as.numeric(effortdays[,2])
#X <- round(t(t(X)*(edays/30)))




f.matrix.creatorLB<-function(data,year){
  #results object
  res<-list()
  
  #get the dimensions of the matrix
  
  #list if sanpling units
  cams<-unique(data$Sampling.Unit.Name)
  cams<-sort(cams)
  rows<-length(cams)
  species<-unique(data$bin)
  #start and end dates of sampling periods
  data<-data[data$Sampling.Period==year,]
  min<-min(data$Start.Date)
  max<-max(data$End.Date)
  cols<-max-min+1
  
  #sampling period
  date.header<-seq(from=min,to=max, by="days")
  mat<-matrix(NA,rows,cols,dimnames=list(cams,as.character(date.header)))
  
  #for all cameras, determine the open and close date and mark in the matrix
  start.dates<-tapply(as.character(data$Start.Date),data$Sampling.Unit.Name,unique)
  nms<-names(start.dates)
  start.dates<-ymd(start.dates)
  names(start.dates)<-nms
  end.dates<-tapply(as.character(data$End.Date),data$Sampling.Unit.Name,unique)
  end.dates<-ymd(end.dates)
  names(end.dates)<-nms
  
  #outline the sampling periods for each camera j
  for(j in 1:length(start.dates)){
    #for each camera beginning and end of sampling
    low<-which(date.header==start.dates[j])
    hi<-which(date.header==end.dates[j])
    if(length(low)+length(hi)>0){
      indx<-seq(from=low,to=hi)
      mat[names(start.dates)[j],indx]<-0
    } else next
  }
  mat.template<-mat
  #get the species
  #species<-unique(data$bin)
  #construct the matrix for each species i
  for(i in 1:length(species)){
    indx<-which(data$bin==species[i])
    #dates and cameras when/where the species was photographed
    dates<-data$Photo.Date[indx]
    cameras<-data$Sampling.Unit.Name[indx]
    dates.cameras<-data.frame(dates,cameras)
    #unique combination of dates and cameras 
    dates.cameras<-unique(dates.cameras)
    #fill in the matrix
    for(j in 1:length(dates.cameras[,1])){
      col<-which(date.header==dates.cameras[j,1])
      row<-which(cams==dates.cameras[j,2])
      mat[row,col]<-1
    }
    mat.nas<-is.na(mat)
    sum.nas<-apply(mat.nas,2,sum)
    indx.nas<-which(sum.nas==rows)
    if(length(indx.nas)>0){
      mat<-mat[,-indx.nas]
    }
    
    res<-c(res,list(mat))
    #return the matrix to its original form
    mat<-mat.template
  }
  
  names(res)<-species
  #res<-lapply(res,f.dum)
  res
  
}

# Change site code to site of interest; check table to see sampling periods per site
table(eventsdata$Site.Code, eventsdata$Sampling.Period)
# Based on sampling periods (which aren't necessarily the dates sampled):
# VB runs 2007-2013 (Photo dates are one year ahead of Sampling.Period; Photos 2008-2014)
# UDZ runs 2009-2013
# BIF runs 2009-2012 (Photo dates are one year ahead of Sampling.Period; Photos 2010-2013)
# PSH runs 2011-2013
# YAN runs 2011-2013
# NAK runs 2009-2012 (Photo dates are one year ahead of Sampling.Period; Photos 2010-2013)
# RNF runs 2010-2013
# NB VB Sampling.Period does not align with actual sampling dates (Photo dates are one year ahead of Sampling.Period)
VBMatrix2008 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2007.01")
VBMatrix2009 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2008.01")
VBMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2009.01")
VBMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2010.01")
VBMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2011.01")
VBMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2012.01")
VBMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2013.01")

UDZMatrix2009 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2009.01")
UDZMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2010.01")
UDZMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2011.01")
UDZMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2012.01")
UDZMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2013.01")

# NB BIF Sampling.Period does not align with actual sampling dates (Photo dates are one year ahead of Sampling.Period)
BIFMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2009.01")
BIFMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2010.01")
BIFMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2011.01")
BIFMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2012.01")

PSHMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2011.01")
PSHMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2012.01")
PSHMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2013.01")

YANMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2011.01")
YANMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2012.01")
YANMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2013.01")

# NB NAK Sampling.Period does not align with actual sampling dates (Photo dates are one year ahead of Sampling.Period)
NAKMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2009.01")
NAKMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2010.01")
NAKMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2011.01")
NAKMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2012.01")

RNFMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2010.01")
RNFMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2011.01")
RNFMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2012.01")
RNFMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2013.01")

#test <- list()
#test[[i]] <- f.shrink.matrix.to15(Matrix2010[[i]])

VBMatrix2008.15 <- list()
VBMatrix2009.15 <- list()
VBMatrix2010.15 <- list()
VBMatrix2011.15 <- list()
VBMatrix2012.15 <- list()
VBMatrix2013.15 <- list()
VBMatrix2014.15 <- list()

UDZMatrix2009.15 <- list()
UDZMatrix2010.15 <- list()
UDZMatrix2011.15 <- list()
UDZMatrix2012.15 <- list()
UDZMatrix2013.15 <- list()

BIFMatrix2010.15 <- list()
BIFMatrix2011.15 <- list()
BIFMatrix2012.15 <- list()
BIFMatrix2013.15 <- list()

PSHMatrix2011.15 <- list()
PSHMatrix2012.15 <- list()
PSHMatrix2013.15 <- list()

YANMatrix2011.15 <- list()
YANMatrix2012.15 <- list()
YANMatrix2013.15 <- list()

NAKMatrix2010.15 <- list()
NAKMatrix2011.15 <- list()
NAKMatrix2012.15 <- list()
NAKMatrix2013.15 <- list()

RNFMatrix2010.15 <- list()
RNFMatrix2011.15 <- list()
RNFMatrix2012.15 <- list()
RNFMatrix2013.15 <- list()

for(i in 1:length(VBMatrix2008)){
  VBMatrix2008.15[[i]] <- f.shrink.matrix.to15(VBMatrix2008[[i]])
  VBMatrix2009.15[[i]] <- f.shrink.matrix.to15(VBMatrix2009[[i]])
  VBMatrix2010.15[[i]] <- f.shrink.matrix.to15(VBMatrix2010[[i]])
  VBMatrix2011.15[[i]] <- f.shrink.matrix.to15(VBMatrix2011[[i]])
  VBMatrix2012.15[[i]] <- f.shrink.matrix.to15(VBMatrix2012[[i]])
  VBMatrix2013.15[[i]] <- f.shrink.matrix.to15(VBMatrix2013[[i]])
  VBMatrix2014.15[[i]] <- f.shrink.matrix.to15(VBMatrix2014[[i]])
}

for(i in 1:length(UDZMatrix2010)){
  UDZMatrix2009.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2009[[i]])
  UDZMatrix2010.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2010[[i]])
  UDZMatrix2011.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2011[[i]])
  UDZMatrix2012.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2012[[i]])
  UDZMatrix2013.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2013[[i]])
}

for(i in 1:length(BIFMatrix2010)){
  BIFMatrix2010.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2010[[i]])
  BIFMatrix2011.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2011[[i]])
  BIFMatrix2012.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2012[[i]])
  BIFMatrix2013.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2013[[i]])
}  

for(i in 1:length(PSHMatrix2012)){
  PSHMatrix2011.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2011[[i]])
  PSHMatrix2012.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2012[[i]])
  PSHMatrix2013.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2013[[i]])
}

for(i in 1:length(YANMatrix2012)){
  YANMatrix2011.15[[i]] <- f.shrink.matrix.to15(YANMatrix2011[[i]])
  YANMatrix2012.15[[i]] <- f.shrink.matrix.to15(YANMatrix2012[[i]])
  YANMatrix2013.15[[i]] <- f.shrink.matrix.to15(YANMatrix2013[[i]])
}

for(i in 1:length(NAKMatrix2010)){
  NAKMatrix2010.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2010[[i]])
  NAKMatrix2011.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2011[[i]])
  NAKMatrix2012.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2012[[i]])
  NAKMatrix2013.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2013[[i]])
}

for(i in 1:length(RNFMatrix2011)){
  RNFMatrix2010.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2010[[i]])
  RNFMatrix2011.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2011[[i]])
  RNFMatrix2012.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2012[[i]])
  RNFMatrix2013.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2013[[i]])
}

# Code used to extract sample data for analysis of 32 species modeled with WPI covariates
# 5 VB species
names(VBMatrix2008)
VB.Pecari_tajacu <- data.frame("2008"=VBMatrix2008.15[[1]], "2009"=VBMatrix2009.15[[1]], "2010"=VBMatrix2010.15[[1]], "2011"=VBMatrix2011.15[[1]], "2012"=VBMatrix2012.15[[1]], "2013"=VBMatrix2013.15[[1]], "2014"=VBMatrix2014.15[[1]])
VB.Dasyprocta_punctata <- data.frame("2008"=VBMatrix2008.15[[3]], "2009"=VBMatrix2009.15[[3]], "2010"=VBMatrix2010.15[[3]], "2011"=VBMatrix2011.15[[3]], "2012"=VBMatrix2012.15[[3]], "2013"=VBMatrix2013.15[[3]], "2014"=VBMatrix2014.15[[3]])
VB.Tapirus_bairdii <- data.frame("2008"=VBMatrix2008.15[[18]], "2009"=VBMatrix2009.15[[18]], "2010"=VBMatrix2010.15[[18]], "2011"=VBMatrix2011.15[[18]], "2012"=VBMatrix2012.15[[18]], "2013"=VBMatrix2013.15[[18]], "2014"=VBMatrix2014.15[[18]])
VB.Cuniculus_paca <- data.frame("2008"=VBMatrix2008.15[[6]], "2009"=VBMatrix2009.15[[6]], "2010"=VBMatrix2010.15[[6]], "2011"=VBMatrix2011.15[[6]], "2012"=VBMatrix2012.15[[6]], "2013"=VBMatrix2013.15[[6]], "2014"=VBMatrix2014.15[[6]])
VB.Dasypus_novemcinctus <- data.frame("2008"=VBMatrix2008.15[[4]], "2009"=VBMatrix2009.15[[4]], "2010"=VBMatrix2010.15[[4]], "2011"=VBMatrix2011.15[[4]], "2012"=VBMatrix2012.15[[4]], "2013"=VBMatrix2013.15[[4]], "2014"=VBMatrix2014.15[[4]])

VB_covariate_species <- list(VB.Pecari_tajacu=VB.Pecari_tajacu,
                             VB.Dasyprocta_punctata=VB.Dasyprocta_punctata,
                             VB.Tapirus_bairdii=VB.Tapirus_bairdii,
                             VB.Cuniculus_paca=VB.Cuniculus_paca,
                             VB.Dasypus_novemcinctus,VB.Dasypus_novemcinctus)

# 10 UDZ species
names(UDZMatrix2009)
UDZ.Guttera_pucherani <- data.frame("2009"=UDZMatrix2009.15[[19]], "2010"=UDZMatrix2010.15[[19]], "2011"=UDZMatrix2011.15[[19]], "2012"=UDZMatrix2012.15[[19]], "2013"=UDZMatrix2013.15[[19]])
UDZ.Bdeogale_crassicauda <- data.frame("2009"=UDZMatrix2009.15[[3]], "2010"=UDZMatrix2010.15[[3]], "2011"=UDZMatrix2011.15[[3]], "2012"=UDZMatrix2012.15[[3]], "2013"=UDZMatrix2013.15[[3]])
UDZ.Cricetomys_gambianus <- data.frame("2009"=UDZMatrix2009.15[[1]], "2010"=UDZMatrix2010.15[[1]], "2011"=UDZMatrix2011.15[[1]], "2012"=UDZMatrix2012.15[[1]], "2013"=UDZMatrix2013.15[[1]])
UDZ.Cephalophus_spadix <- data.frame("2009"=UDZMatrix2009.15[[6]], "2010"=UDZMatrix2010.15[[6]], "2011"=UDZMatrix2011.15[[6]], "2012"=UDZMatrix2012.15[[6]], "2013"=UDZMatrix2013.15[[6]])
UDZ.Genetta_servalina <- data.frame("2009"=UDZMatrix2009.15[[20]], "2010"=UDZMatrix2010.15[[20]], "2011"=UDZMatrix2011.15[[20]], "2012"=UDZMatrix2012.15[[20]], "2013"=UDZMatrix2013.15[[20]])
UDZ.Cephalophus_harveyi <- data.frame("2009"=UDZMatrix2009.15[[5]], "2010"=UDZMatrix2010.15[[5]], "2011"=UDZMatrix2011.15[[5]], "2012"=UDZMatrix2012.15[[5]], "2013"=UDZMatrix2013.15[[5]])
UDZ.Paraxerus_vexillarius <- data.frame("2009"=UDZMatrix2009.15[[18]], "2010"=UDZMatrix2010.15[[18]], "2011"=UDZMatrix2011.15[[18]], "2012"=UDZMatrix2012.15[[18]], "2013"=UDZMatrix2013.15[[18]])
UDZ.Cercocebus_sanjei <- data.frame("2009"=UDZMatrix2009.15[[4]], "2010"=UDZMatrix2010.15[[4]], "2011"=UDZMatrix2011.15[[4]], "2012"=UDZMatrix2012.15[[4]], "2013"=UDZMatrix2013.15[[4]])
UDZ.Rhynchocyon_udzungwensis <- data.frame("2009"=UDZMatrix2009.15[[30]], "2010"=UDZMatrix2010.15[[30]], "2011"=UDZMatrix2011.15[[30]], "2012"=UDZMatrix2012.15[[30]], "2013"=UDZMatrix2013.15[[30]])
UDZ.Nesotragus_moschatus <- data.frame("2009"=UDZMatrix2009.15[[11]], "2010"=UDZMatrix2010.15[[11]], "2011"=UDZMatrix2011.15[[11]], "2012"=UDZMatrix2012.15[[11]], "2013"=UDZMatrix2013.15[[11]])

UDZ_covariate_species <- list(UDZ.Guttera_pucherani=UDZ.Guttera_pucherani, 
                              UDZ.Bdeogale_crassicauda=UDZ.Bdeogale_crassicauda, 
                              UDZ.Cricetomys_gambianus=UDZ.Cricetomys_gambianus, 
                              UDZ.Cephalophus_spadix=UDZ.Cephalophus_spadix, 
                              UDZ.Genetta_servalina=UDZ.Genetta_servalina, 
                              UDZ.Cephalophus_harveyi=UDZ.Cephalophus_harveyi, 
                              UDZ.Paraxerus_vexillarius=UDZ.Paraxerus_vexillarius, 
                              UDZ.Cercocebus_sanjei,UDZ.Cercocebus_sanjei,
                              UDZ.Rhynchocyon_udzungwensis=UDZ.Rhynchocyon_udzungwensis, 
                              UDZ.Nesotragus_moschatus=UDZ.Nesotragus_moschatus)

# 3 BIF species
names(BIFMatrix2011)
BIF.Cercopithecus_lhoesti <- data.frame("2010"=BIFMatrix2010.15[[8]], "2011"=BIFMatrix2011.15[[8]], "2012"=BIFMatrix2012.15[[8]], "2013"=BIFMatrix2013.15[[8]])
BIF.Cephalophus_silvicultor <- data.frame("2010"=BIFMatrix2010.15[[6]], "2011"=BIFMatrix2011.15[[6]], "2012"=BIFMatrix2012.15[[6]], "2013"=BIFMatrix2013.15[[6]])
BIF.Cephalophus_nigrifrons <- data.frame("2010"=BIFMatrix2010.15[[1]], "2011"=BIFMatrix2011.15[[1]], "2012"=BIFMatrix2012.15[[1]], "2013"=BIFMatrix2013.15[[1]])

BIF_covariate_species <- list(BIF.Cercopithecus_lhoesti=BIF.Cercopithecus_lhoesti,
                              BIF.Cephalophus_silvicultor=BIF.Cephalophus_silvicultor,
                              BIF.Cephalophus_nigrifrons=BIF.Cephalophus_nigrifrons)
                              

# 5 PSH species
names(PSHMatrix2011)

PSH.Leopoldamys_sabanus <- data.frame("2011"=PSHMatrix2011.15[[18]], "2012"=PSHMatrix2012.15[[18]], "2013"=PSHMatrix2013.15[[18]])
PSH.Muntiacus_muntjak <- data.frame("2011"=PSHMatrix2011.15[[15]], "2012"=PSHMatrix2012.15[[15]], "2013"=PSHMatrix2013.15[[15]])
PSH.Macaca_nemestrina <- data.frame("2011"=PSHMatrix2011.15[[2]], "2012"=PSHMatrix2012.15[[2]], "2013"=PSHMatrix2013.15[[2]])
PSH.Sus_scrofa <- data.frame("2011"=PSHMatrix2011.15[[6]], "2012"=PSHMatrix2012.15[[6]], "2013"=PSHMatrix2013.15[[6]])
PSH.Tragulus_kanchil <- data.frame("2011"=PSHMatrix2011.15[[9]], "2012"=PSHMatrix2012.15[[9]], "2013"=PSHMatrix2013.15[[9]])

PSH_covariate_species <- list(PSH.Leopoldamys_sabanus=PSH.Leopoldamys_sabanus,
                              PSH.Muntiacus_muntjak=PSH.Muntiacus_muntjak,
                              PSH.Macaca_nemestrina=PSH.Macaca_nemestrina,
                              PSH.Sus_scrofa=PSH.Sus_scrofa,
                              PSH.Tragulus_kanchil=PSH.Tragulus_kanchil)

# 5 YAN species
names(YANMatrix2011)
YAN.Cuniculus_paca <- data.frame("2011"=YANMatrix2011.15[[3]], "2012"=YANMatrix2012.15[[3]], "2013"=YANMatrix2013.15[[3]])
YAN.Dasyprocta_fuliginosa <- data.frame("2011"=YANMatrix2011.15[[7]], "2012"=YANMatrix2012.15[[7]], "2013"=YANMatrix2013.15[[7]])
YAN.Dasypus_novemcinctus <- data.frame("2011"=YANMatrix2011.15[[5]], "2012"=YANMatrix2012.15[[5]], "2013"=YANMatrix2013.15[[5]])
YAN.Mitu_tuberosum <- data.frame("2011"=YANMatrix2011.15[[2]], "2012"=YANMatrix2012.15[[2]], "2013"=YANMatrix2013.15[[2]])
YAN.Tapirus_terrestris <- data.frame("2011"=YANMatrix2011.15[[11]], "2012"=YANMatrix2012.15[[11]], "2013"=YANMatrix2013.15[[11]])

YAN_covariate_species <- list(YAN.Cuniculus_paca=YAN.Cuniculus_paca,
                              YAN.Dasyprocta_fuliginosa=YAN.Dasyprocta_fuliginosa,
                              YAN.Dasypus_novemcinctus=YAN.Dasypus_novemcinctus,
                              YAN.Mitu_tuberosum=YAN.Mitu_tuberosum,
                              YAN.Tapirus_terrestris=YAN.Tapirus_terrestris)

# 2 NAK species
names(NAKMatrix2011)
NAK.Muntiacus_muntjak <- data.frame("2010"=NAKMatrix2010.15[[1]], "2011"=NAKMatrix2011.15[[1]], "2012"=NAKMatrix2012.15[[1]], "2013"=NAKMatrix2013.15[[1]])
NAK.Atherurus_macrourus <- data.frame("2010"=NAKMatrix2010.15[[10]], "2011"=NAKMatrix2011.15[[10]], "2012"=NAKMatrix2012.15[[10]], "2013"=NAKMatrix2013.15[[10]])

NAK_covariate_species <- list(NAK.Muntiacus_muntjak=NAK.Muntiacus_muntjak,
                              NAK.Atherurus_macrourus=NAK.Atherurus_macrourus)

# 2 RNF species
names(RNFMatrix2011)
RNF.Nesomys_rufus <- data.frame("2010"=RNFMatrix2010.15[[4]], "2011"=RNFMatrix2011.15[[4]], "2012"=RNFMatrix2012.15[[4]], "2013"=RNFMatrix2013.15[[4]])
RNF.Fossa_fossana <- data.frame("2010"=RNFMatrix2010.15[[3]], "2011"=RNFMatrix2011.15[[3]], "2012"=RNFMatrix2012.15[[3]], "2013"=RNFMatrix2013.15[[3]])

RNF_covariate_species <- list(RNF.Nesomys_rufus=RNF.Nesomys_rufus,
                              RNF.Fossa_fossana=RNF.Fossa_fossana)







# Code used to extract sample data for Miguel Acevedo for Skype meeting on September 5, 2014; modified for analysis of 32 species modeled with WPI covariates
#names(Matrix2010)
#UDZ.crestedguineafowl <- cbind(Matrix2010[[19]], Matrix2011[[19]], Matrix2012[[19]], Matrix2013[[19]], Matrix2014[[19]])
#UDZ.mongoose <- cbind(Matrix2010[[3]], Matrix2011[[3]], Matrix2012[[3]], Matrix2013[[3]], Matrix2014[[3]])
#UDZ.giantrat <- cbind(Matrix2010[[1]], Matrix2011[[1]], Matrix2012[[1]], Matrix2013[[1]], Matrix2014[[1]])
#UDZ.Abbottsduiker <- cbind(Matrix2010[[6]], Matrix2011[[6]], Matrix2012[[6]], Matrix2013[[6]], Matrix2014[[6]])
#UDZ.genet <- cbind(Matrix2010[[20]], Matrix2011[[20]], Matrix2012[[20]], Matrix2013[[20]], Matrix2014[[20]])
#UDZ.Harveysduiker <- cbind(Matrix2010[[5]], Matrix2011[[5]], Matrix2012[[5]], Matrix2013[[5]], Matrix2014[[5]])
#UDZ.bushsquirrel <- cbind(Matrix2010[[18]], Matrix2011[[18]], Matrix2012[[18]], Matrix2013[[18]], Matrix2014[[18]])
#UDZ.mangabey <- cbind(Matrix2010[[4]], Matrix2011[[4]], Matrix2012[[4]], Matrix2013[[4]], Matrix2014[[4]])
#UDZ.elephantshrew <- cbind(Matrix2010[[30]], Matrix2011[[30]], Matrix2012[[30]], Matrix2013[[30]], Matrix2014[[30]])
#UDZ.suni <- cbind(Matrix2010[[11]], Matrix2011[[11]], Matrix2012[[11]], Matrix2013[[11]], Matrix2014[[11]])
# Not modeled with covariates in WPI analysis
##UDZ.bluemonkey <- cbind(Matrix2010[[22]], Matrix2011[[22]], Matrix2012[[22]], Matrix2013[[22]], Matrix2014[[22]])
##UDZ.elephant <- cbind(Matrix2010[[14]], Matrix2011[[14]], Matrix2012[[14]], Matrix2013[[14]], Matrix2014[[14]])
##UDZ.honeybadger <- cbind(Matrix2010[[23]], Matrix2011[[23]], Matrix2012[[23]], Matrix2013[[23]], Matrix2014[[23]])
##UDZ.palmcivet <- cbind(Matrix2010[[16]], Matrix2011[[16]], Matrix2012[[16]], Matrix2013[[16]], Matrix2014[[16]])
##UDZ.bushpig <- cbind(Matrix2010[[10]], Matrix2011[[10]], Matrix2012[[10]], Matrix2013[[10]], Matrix2014[[10]])

write.csv(UDZ.crestedguineafowl, file="UDZ.crestedguineafowl_Guttera_pucherani.csv")
write.csv(UDZ.mongoose, file="UDZ.mongoose_Bdeogale_crassicauda.csv")
write.csv(UDZ.giantrat, file="UDZ.giantrat_Cricetomys_gambianus.csv")
write.csv(UDZ.Abbottsduiker, file="UDZ.Abbottsduiker_Cephalophus_spadix.csv")
write.csv(UDZ.genet, file="UDZ.genet_Genetta_servalina.csv")
write.csv(UDZ.Harveysduiker, file="UDZ.Harveysduiker_Cephalophus_harveyi.csv")
write.csv(UDZ.bushsquirrel, file="UDZ.bushsquirrel_Paraxerus_vexillarius.csv")
write.csv(UDZ.mangabey, file="UDZ.mangabey_Cercocebus_sanjei.csv")
write.csv(UDZ.elephantshrew, file="UDZ.elephantshrew_Rhynchocyon_udzungwensis.csv")
write.csv(UDZ.suni, file="UDZ.suni_Nesotragus_moschatus.csv")
# UDZ species not modeled with covariates in WPI analysis
##write.csv(UDZ.bluemonkey, file="UDZ.bluemonkey_Cercopithecus_mitis.csv")
##write.csv(UDZ.elephant, file="UDZ.elephant_Loxodonta_africana.csv")
##write.csv(UDZ.honeybadger, file="UDZ.honeybadger_Mellivora_capensis.csv")
##write.csv(UDZ.palmcivet, file="UDZ.palmcivet_Nandinia_binotata.csv")
##write.csv(UDZ.bushpig, file="UDZ.bushpig_Potamochoerus_larvatus.csv")




# Code used to extract sample data for Miguel Acevedo for Skype meeting on August 22, 2014
#VB.pecari <- cbind(Matrix2008[[1]], Matrix2009[[1]], Matrix2010[[1]], Matrix2011[[1]], Matrix2012[[1]], Matrix2013[[1]], Matrix2014[[1]])
#VB.agouti <- cbind(Matrix2008[[3]], Matrix2009[[3]], Matrix2010[[3]], Matrix2011[[3]], Matrix2012[[3]], Matrix2013[[3]], Matrix2014[[3]])
#VB.Bairdstapir <- cbind(Matrix2008[[18]], Matrix2009[[18]], Matrix2010[[18]], Matrix2011[[18]], Matrix2012[[18]], Matrix2013[[18]], Matrix2014[[18]])
#VB.paca <- cbind(Matrix2008[[6]], Matrix2009[[6]], Matrix2010[[6]], Matrix2011[[6]], Matrix2012[[6]], Matrix2013[[6]], Matrix2014[[6]])
#VB.armadillo <- cbind(Matrix2008[[4]], Matrix2009[[4]], Matrix2010[[4]], Matrix2011[[4]], Matrix2012[[4]], Matrix2013[[4]], Matrix2014[[4]])
# VB species not modeled with covariates in WPI analysis
#VB.puma <- cbind(Matrix2008[[14]], Matrix2009[[14]], Matrix2010[[14]], Matrix2011[[14]], Matrix2012[[14]], Matrix2013[[14]], Matrix2014[[14]])
#VB.coati <- cbind(Matrix2008[[15]], Matrix2009[[15]], Matrix2010[[15]], Matrix2011[[15]], Matrix2012[[15]], Matrix2013[[15]], Matrix2014[[15]])
#VB.opossum <- cbind(Matrix2008[[5]], Matrix2009[[5]], Matrix2010[[5]], Matrix2011[[5]], Matrix2012[[5]], Matrix2013[[5]], Matrix2014[[5]])
#VB.ocelot <- cbind(Matrix2008[[12]], Matrix2009[[12]], Matrix2010[[12]], Matrix2011[[12]], Matrix2012[[12]], Matrix2013[[12]], Matrix2014[[12]])
#VB.jaguar <- cbind(Matrix2008[[31]], Matrix2009[[31]], Matrix2010[[31]], Matrix2011[[31]], Matrix2012[[31]], Matrix2013[[31]], Matrix2014[[31]])
#VB.tinamou <- cbind(Matrix2008[[2]], Matrix2009[[2]], Matrix2010[[2]], Matrix2011[[2]], Matrix2012[[2]], Matrix2013[[2]], Matrix2014[[2]])
#VB.redbrocket <- cbind(Matrix2008[[19]], Matrix2009[[19]], Matrix2010[[19]], Matrix2011[[19]], Matrix2012[[19]], Matrix2013[[19]], Matrix2014[[19]])
#VB.curassow <- cbind(Matrix2008[[8]], Matrix2009[[8]], Matrix2010[[8]], Matrix2011[[8]], Matrix2012[[8]], Matrix2013[[8]], Matrix2014[[8]])

write.csv(VB.pecari, file="VB.pecari_Pecari_tajacu.csv")
write.csv(VB.agouti, file="VB.agouti_Dasyprocta_punctata.csv")
write.csv(VB.Bairdstapir, file="VB.Bairdstapir_Tapirus_bairdii.csv")
write.csv(VB.paca, file="VB.paca_Cuniculus_paca.csv")
write.csv(VB.armadillo, file="VB.armadillo_Dasypus_novemcinctus.csv")
# VB species not modeled with covariates in WPI analysis
##write.csv(VB.puma, file="VB.puma_Puma_concolor.csv")
##write.csv(VB.coati, file="VB.coati_Nasua_narica.csv")
##write.csv(VB.opossum, file="VB.opossum_Didelphis_marsupialis.csv")
##write.csv(VB.ocelot, file="VB.ocelot_Leopardus_pardalis.csv")
##write.csv(VB.jaguar, file="VB.jaguar_Panthera_onca.csv")
##write.csv(VB.tinamou, file="VB.tinamou_Tinamus_major.csv")
##write.csv(VB.redbrocket, file="VB.redbrocket_Mazama_temama.csv")
##write.csv(VB.curassow, file="VB.curassow_Crax_rubra.csv")



# Code used to extract sample data for Miguel Acevedo for initial Skype meeting on August 8, 2014

Pecari.2008 <- Matrix2008[[1]]
Pecari.2009 <- Matrix2009[[1]]
Pecari.2010 <- Matrix2010[[1]]
Pecari.2011 <- Matrix2011[[1]]
Pecari.2012 <- Matrix2012[[1]]
Pecari.2013 <- Matrix2013[[1]]
Pecari.2014 <- Matrix2014[[1]]


Pecari.2008c <- rbind(Pecari.2008[1:20,2:31], Pecari.2008[21:40,37:66], Pecari.2008[41:60,80:109])
colnames(Pecari.2008c) <- paste0(rep("2008.Day",30),1:30)
Pecari.2009c <- rbind(Pecari.2009[1:20,2:31], Pecari.2009[21:40,37:66], Pecari.2009[41:60,80:109])
colnames(Pecari.2009c) <- paste0(rep("2009.Day",30),1:30)
Pecari.2010c <- rbind(Pecari.2010[1:20,2:31], Pecari.2010[21:40,37:66], Pecari.2010[41:60,72:101])
colnames(Pecari.2010c) <- paste0(rep("2010.Day",30),1:30)
Pecari.2011c <- rbind(Pecari.2011[1:20,2:31], Pecari.2011[21:40,37:66], Pecari.2011[41:60,72:101])
colnames(Pecari.2011c) <- paste0(rep("2011.Day",30),1:30)
Pecari.2012c <- rbind(Pecari.2012[1:20,2:31], Pecari.2012[21:40,34:63], Pecari.2012[41:60,72:101])
colnames(Pecari.2012c) <- paste0(rep("2012.Day",30),1:30)
Pecari.2013c <- rbind(Pecari.2013[1:20,2:31], Pecari.2013[21:40,33:62], Pecari.2013[41:60,71:100])
colnames(Pecari.2013c) <- paste0(rep("2013.Day",30),1:30)
Pecari.2014c <- rbind(Pecari.2014[1:20,2:31], Pecari.2014[21:40,40:69], Pecari.2014[41:60,71:100])
colnames(Pecari.2014c) <- paste0(rep("2014.Day",30),1:30)

Pecari_allyears <- cbind(Pecari.2008c, Pecari.2009c, Pecari.2010c, Pecari.2011c, Pecari.2012c, Pecari.2013c, Pecari.2014c)



Pecari <- subset(eventsdata, Site.Code=="VB-" & bin=="Pecari tajacu")
Dasyprocta <- subset(eventsdata, Site.Code=="VB-" & bin=="Dasyprocta punctata")

species <- Pecari

X <- table(species$Sampling.Unit.Name, species$Photo.Date)
X <- as.matrix(ifelse(X>0,1,0))


ELEV <- read.csv("CT_edgedist_elevation_final.csv")
ELEVsub <- ELEV[match(rownames(X), ELEV$Sampling.Unit.Name),]

X <- cbind(Elevation=ELEVsub$Elevation, EdgeDist=ELEVsub$EdgeDist, X)

#write.csv(X, file="VolcanBarva_Pecari.tajacu.csv")
#write.csv(X, file="VolcanBarva_Dasyprocta.punctata.csv")

#Next step is to include the camera traps and dates with no observations (above table only includes obs)
#Also, decide whether to collapse the number of columns down into Day1...Day30 rather than keeping dates
#for now, use same approach to control for the variation in # of days as used in carbon policy paper

test <- f.matrix.creator2(alldata,"2008.01")




################## Camera Trap Temperature Data Formatting ####################
#Examine distribution of NA temperature values across camera traps
table(eventsdata$Sampling.Unit.Name, is.na(eventsdata$Temperature)==TRUE)

# Clean temperature data and convert F measurements to C
library(stringr)
eventsdata$Temperature <- str_trim(eventsdata$Temperature, side="both")
temp.degrees <- sub(" .*","", eventsdata$Temperature)
temp.unit <- str_sub(eventsdata$Temperature, nchar(eventsdata$Temperature), nchar(eventsdata$Temperature))
temp.unit <- str_trim(temp.unit, side="both")
#Convert F temperature values to Celcius
#Change C Temperatures to numeric format from character
temp.degreesC <- as.numeric(ifelse(temp.unit=="F", f.FtoC(as.numeric(temp.degrees)), temp.degrees))
eventsdata <- cbind(eventsdata, temp.degreesC)


# Determine the annual min, max and variance of the non-calibrated temperature data for each CT

Temp.Min <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=min)
names(Temp.Min) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Min")
Temp.Max <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=max)
names(Temp.Max) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Max")
Temp.Var <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=var)
names(Temp.Var) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Var")

CT.Temp <- cbind(Temp.Min, Temp.Max$Temp.Max, Temp.Var$Temp.Var)
names(CT.Temp) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Min", "Temp.Max", "Temp.Var")

CT.Temp.VB <- melt(CT.Temp[CT.Temp$Site.Code=="VB-",])
Tmin <- as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Sampling.Period ~ variable)[,,1])
names(Tmin) <- paste0("Tmin.",2008:2014)
Tmax <- as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Sampling.Period ~ variable)[,,2])
names(Tmax) <- paste0("Tmax.",2008:2014)
Tvar <- round(as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Sampling.Period ~ variable)[,,3]),2)
names(Tvar) <- paste0("Tvar.",2008:2014)

ELEV <- read.csv("CT_edgedist_elevation_final.txt")
ELEVsub <- ELEV[match(rownames(Tmin), ELEV$Sampling.Unit.Name),]

Pecari.Data <- cbind(Elevation=ELEVsub$Elevation, Tmin, Tmax, Tvar, Pecari_allyears)
#write.csv(Pecari.Data, file="VB_allyears.temp_Pecari.tajacu.csv")


CT.Temp.UDZ <- melt(CT.Temp[CT.Temp$Site.Code=="UDZ",])
Tmin <- as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Sampling.Period ~ variable)[,,1])
names(Tmin) <- paste0("Tmin.",2010:2014)
Tmax <- as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Sampling.Period ~ variable)[,,2])
names(Tmax) <- paste0("Tmax.",2010:2014)
Tvar <- round(as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Sampling.Period ~ variable)[,,3]),2)
names(Tvar) <- paste0("Tvar.",2010:2014)
UDZ.CovariateData <- cbind(Elevation=ELEVsub$Elevation, Tmin, Tmax, Tvar)
#write.csv(UDZ.CovariateData, file="UDZ.CovariateData.csv")
