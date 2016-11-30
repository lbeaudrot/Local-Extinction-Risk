rm(list=ls())

library(lubridate)
library(reshape)
library(stringr)
library(plyr)
library(fields)


# Load Jorge's file "camera trap analysis functions.R"
source("camera trap analysis functions.R")
# Load matrix creator function at bottom of this script "f.matrix.creatorLB"
source("matrix creator LB.R")
# Load TEAM data, fix data, and add site code; no need to separate events because matrix creator functions turns all detections into trinary matrix
#ctdata <- f.teamdb.query("camera trap")
#load("ct_data2014-10-31.gzip")
#load("ct_data2015-07-30.gzip")
load("ct_data2016-11-26.gzip")

alldata <- cam_trap_data
alldata<-f.fix.data2(alldata)
Site.Code <- substr(alldata$Sampling.Unit.Name,4,6)
alldata <- cbind(alldata, Site.Code)


splist<-read.csv("master_species_list_updated_7April2014.csv",h=T) #master list
sitelist<-unique(alldata$bin) 
newsplist<-subset(splist,splist$Unique_Name %in% sitelist & splist$Include==1)
subdata<-subset(alldata, alldata$bin %in% newsplist$Unique_Name) #this is the original camera trap data subsetted to these species
subdata<-f.correct.DF(subdata)
eventsdata <- subdata

# Reduce dataset to mammals only
eventsdata <- eventsdata[eventsdata$Class=="MAMMALIA",]

# Try ordering data by new site-species column to see if downstream data will all follow alphabetical order
SITE.SP <- paste(eventsdata$Site.Code, eventsdata$bin, sep="-")
eventsdata <- data.frame(eventsdata, SITE.SP=SITE.SP)
eventsdata <- eventsdata[order(eventsdata$SITE.SP),]

# Remove Atelocynus microtis from dataset so that new 2014 YAN observations are not included and Speothos venaticus too.
eventsdata <- eventsdata[eventsdata$bin!="Atelocynus microtis",]
eventsdata <- eventsdata[eventsdata$bin!="Speothos venaticus",]

# Need to add Muntiacus montanus so that it appears in the UID

# Remove CT-YAN-1-04 (inactive camera trap)
eventsdata <- eventsdata[eventsdata$Sampling.Unit.Name!="CT-YAN-1-04",]

##################### CREATE INPUT DATA FOR UNMARKED ANALYSIS WITH 15 SECONDARY SAMPLING PERIODS #################
# Create matrices with 15 secondary sampling periods for each year of data collection for each species for sites with >500 m elevation gradients
# Check table to see sampling periods per site. Based on sampling periods:
table(eventsdata$Site.Code, eventsdata$Sampling.Period)
# VB runs 2007-2013 (Photo dates are one year ahead of Sampling.Period; Photos 2008-2014)
# UDZ runs 2009-2013
# BIF runs 2009-2012 (Photo dates are one year ahead of Sampling.Period; Photos 2010-2013)
# PSH runs 2011-2013
# YAN runs 2011-2013
# NAK runs 2009-2012 (Photo dates are one year ahead of Sampling.Period; Photos 2010-2013)
# RNF runs 2010-2013

# Create lists of full sized matrices for each site and year
  # NB VB Sampling.Period does not align with actual sampling dates (Photo dates are one year ahead of Sampling.Period)
  VBMatrix2008 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2007.01")
  VBMatrix2009 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2008.01")
  VBMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2009.01")
  VBMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2010.01")
  VBMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2011.01")
  VBMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2012.01")
  VBMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2013.01")
  VBMatrix2015 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2014.01") #Added 7/30/2015
  VBMatrix2016 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="VB-",], "2015.01") #Added 11/26/2016


  UDZMatrix2009 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2009.01")
  UDZMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2010.01")
  UDZMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2011.01")
  UDZMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2012.01")
  UDZMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2013.01")
  UDZMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2014.01")
  #UDZMatrix2015 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="UDZ",], "2015.01") #Added 11/26/2016; Removed bc temp data missing



  # NB BIF Sampling.Period does not align with actual sampling dates (Photo dates are one year ahead of Sampling.Period)
  BIFMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2009.01")
  BIFMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2010.01")
  BIFMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2011.01")
  BIFMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2012.01")
  BIFMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2013.01") #Added 11/28/2016
  BIFMatrix2015 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2014.01") #Added 7/30/2015
  #BIFMatrix2016 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="BIF",], "2015.01") #Added 11/26/2016; Removed bc temp data missing


  PSHMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2011.01")
  PSHMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2012.01")
  PSHMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2013.01")
  PSHMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2014.01")
  PSHMatrix2015 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2015.01")
  PSHMatrix2016 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="PSH",], "2016.01") #Added 11/26/2016



  YANMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2011.01")
  YANMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2012.01")
  YANMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2013.01")
  YANMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2014.01") #Added 7/30/2015
  YANMatrix2015 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="YAN",], "2015.01") #Added 11/26/2016
  
  
  # NB NAK Sampling.Period does not align with actual sampling dates (Photo dates are one year ahead of Sampling.Period)
  NAKMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2009.01")
  NAKMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2010.01")
  NAKMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2011.01")
  NAKMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2012.01")
  NAKMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2013.01") #Changed from 2014.01 to 2013.01 on 7/30/2015
  #NAKMatrix2015 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="NAK",], "2014.01") #Added 11/26/2016; Removed bc temp data missing


  RNFMatrix2010 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2010.01")
  RNFMatrix2011 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2011.01")
  RNFMatrix2012 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2012.01")
  RNFMatrix2013 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2013.01")
  RNFMatrix2014 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2014.01") #Added 11/26/2016 - few obs? check with Eileen
  RNFMatrix2015 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2015.01") #Added 11/26/2016 - few obs? check with Eileen
  #RNFMatrix2016 <- f.matrix.creatorLB(eventsdata[eventsdata$Site.Code=="RNF",], "2016.01") #Added 11/26/2016 - few obs? check with Eileen; Removed bc temp data only available for 6 CTs


# Create empty objects to fill with matrices collapsed to 15 secondary sampling periods
  VBMatrix2008.15 <- list()
  VBMatrix2009.15 <- list()
  VBMatrix2010.15 <- list()
  VBMatrix2011.15 <- list()
  VBMatrix2012.15 <- list()
  VBMatrix2013.15 <- list()
  VBMatrix2014.15 <- list()
  VBMatrix2015.15 <- list()
  VBMatrix2016.15 <- list()



  UDZMatrix2009.15 <- list()
  UDZMatrix2010.15 <- list()
  UDZMatrix2011.15 <- list()
  UDZMatrix2012.15 <- list()
  UDZMatrix2013.15 <- list()
  UDZMatrix2014.15 <- list()
  #UDZMatrix2015.15 <- list()



  BIFMatrix2010.15 <- list()
  BIFMatrix2011.15 <- list()
  BIFMatrix2012.15 <- list()
  BIFMatrix2013.15 <- list()
  BIFMatrix2014.15 <- list()
  BIFMatrix2015.15 <- list()
  #BIFMatrix2016.15 <- list()




  PSHMatrix2011.15 <- list()
  PSHMatrix2012.15 <- list()
  PSHMatrix2013.15 <- list()
  PSHMatrix2014.15 <- list()
  PSHMatrix2015.15 <- list()
  PSHMatrix2016.15 <- list()



  YANMatrix2011.15 <- list()
  YANMatrix2012.15 <- list()
  YANMatrix2013.15 <- list()
  YANMatrix2014.15 <- list()
  YANMatrix2015.15 <- list()



  NAKMatrix2010.15 <- list()
  NAKMatrix2011.15 <- list()
  NAKMatrix2012.15 <- list()
  NAKMatrix2013.15 <- list()
  NAKMatrix2014.15 <- list()
  #NAKMatrix2015.15 <- list()



  RNFMatrix2010.15 <- list()
  RNFMatrix2011.15 <- list()
  RNFMatrix2012.15 <- list()
  RNFMatrix2013.15 <- list()
  RNFMatrix2014.15 <- list()
  RNFMatrix2015.15 <- list()
  #RNFMatrix2016.15 <- list()

# Shrink data matrices to 15 secondary sampling periods
for(i in 1:length(VBMatrix2008)){
  VBMatrix2008.15[[i]] <- f.shrink.matrix.to15(VBMatrix2008[[i]])
  VBMatrix2009.15[[i]] <- f.shrink.matrix.to15(VBMatrix2009[[i]])
  VBMatrix2010.15[[i]] <- f.shrink.matrix.to15(VBMatrix2010[[i]])
  VBMatrix2011.15[[i]] <- f.shrink.matrix.to15(VBMatrix2011[[i]])
  VBMatrix2012.15[[i]] <- f.shrink.matrix.to15(VBMatrix2012[[i]])
  VBMatrix2013.15[[i]] <- f.shrink.matrix.to15(VBMatrix2013[[i]])
  VBMatrix2014.15[[i]] <- f.shrink.matrix.to15(VBMatrix2014[[i]])
  VBMatrix2015.15[[i]] <- f.shrink.matrix.to15(VBMatrix2015[[i]])
  VBMatrix2016.15[[i]] <- f.shrink.matrix.to15(VBMatrix2016[[i]])
}

for(i in 1:length(UDZMatrix2010)){
  UDZMatrix2009.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2009[[i]])
  UDZMatrix2010.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2010[[i]])
  UDZMatrix2011.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2011[[i]])
  UDZMatrix2012.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2012[[i]])
  UDZMatrix2013.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2013[[i]])
  UDZMatrix2014.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2014[[i]])
  #UDZMatrix2015.15[[i]] <- f.shrink.matrix.to15(UDZMatrix2015[[i]])
}

for(i in 1:length(BIFMatrix2010)){
  BIFMatrix2010.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2010[[i]])
  BIFMatrix2011.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2011[[i]])
  BIFMatrix2012.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2012[[i]])
  BIFMatrix2013.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2013[[i]])
  BIFMatrix2014.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2014[[i]])
  BIFMatrix2015.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2015[[i]])
  #BIFMatrix2016.15[[i]] <- f.shrink.matrix.to15(BIFMatrix2016[[i]])

}  

for(i in 1:length(PSHMatrix2012)){
  PSHMatrix2011.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2011[[i]])
  PSHMatrix2012.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2012[[i]])
  PSHMatrix2013.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2013[[i]])
  PSHMatrix2014.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2014[[i]])
  PSHMatrix2015.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2015[[i]])
  PSHMatrix2016.15[[i]] <- f.shrink.matrix.to15(PSHMatrix2016[[i]])
}

for(i in 1:length(YANMatrix2012)){
  YANMatrix2011.15[[i]] <- f.shrink.matrix.to15(YANMatrix2011[[i]])
  YANMatrix2012.15[[i]] <- f.shrink.matrix.to15(YANMatrix2012[[i]])
  YANMatrix2013.15[[i]] <- f.shrink.matrix.to15(YANMatrix2013[[i]])
  YANMatrix2014.15[[i]] <- f.shrink.matrix.to15(YANMatrix2014[[i]])
  YANMatrix2015.15[[i]] <- f.shrink.matrix.to15(YANMatrix2015[[i]])
}

for(i in 1:length(NAKMatrix2010)){
  NAKMatrix2010.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2010[[i]])
  NAKMatrix2011.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2011[[i]])
  NAKMatrix2012.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2012[[i]])
  NAKMatrix2013.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2013[[i]])
  NAKMatrix2014.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2014[[i]])
  #NAKMatrix2015.15[[i]] <- f.shrink.matrix.to15(NAKMatrix2015[[i]])
}

for(i in 1:length(RNFMatrix2011)){
  RNFMatrix2010.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2010[[i]])
  RNFMatrix2011.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2011[[i]])
  RNFMatrix2012.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2012[[i]])
  RNFMatrix2013.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2013[[i]])
  RNFMatrix2014.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2014[[i]])
  RNFMatrix2015.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2015[[i]])
  #RNFMatrix2016.15[[i]] <- f.shrink.matrix.to15(RNFMatrix2016[[i]])
}

# Extract all species for each site

    # Extract all VB species
    VB_.species <- list()
    for(i in 1:length(VBMatrix2008)){
      VB_.species[[i]] <- data.frame("2008"=VBMatrix2008.15[[i]], "2009"=VBMatrix2009.15[[i]], "2010"=VBMatrix2010.15[[i]], 
                                     "2011"=VBMatrix2011.15[[i]], "2012"=VBMatrix2012.15[[i]], "2013"=VBMatrix2013.15[[i]], 
                                     "2014"=VBMatrix2014.15[[i]], "2015"=VBMatrix2015.15[[i]], "2016"=VBMatrix2016.15[[i]])
    }
    names(VB_.species) <- paste("VB_", names(VBMatrix2008), sep=".")

    # Extract all UDZ species
    UDZ.species <- list()
    for(i in 1:length(UDZMatrix2009)){
      UDZ.species[[i]] <- data.frame("2009"=UDZMatrix2009.15[[i]], "2010"=UDZMatrix2010.15[[i]], "2011"=UDZMatrix2011.15[[i]], 
                                     "2012"=UDZMatrix2012.15[[i]], "2013"=UDZMatrix2013.15[[i]], "2014"=UDZMatrix2014.15[[i]])
                                     #"2015"=UDZMatrix2015.15[[i]])
    }
    names(UDZ.species) <- paste("UDZ", names(UDZMatrix2009), sep=".")

    # Extract all BIF species
    BIF.species <- list()
    for(i in 1:length(BIFMatrix2011)){
      BIF.species[[i]] <- data.frame("2010"=BIFMatrix2010.15[[i]], "2011"=BIFMatrix2011.15[[i]], "2012"=BIFMatrix2012.15[[i]], 
                                     "2013"=BIFMatrix2013.15[[i]], "2014"=BIFMatrix2014.15[[i]], "2015"=BIFMatrix2015.15[[i]])
                                     #"2016"=BIFMatrix2016.15[[i]])
    }
    names(BIF.species) <- paste("BIF", names(BIFMatrix2011), sep=".")
    
## All Pasoh data (camera trap and covariates) need to be limited to Array 1 only
    # Extract all PSH species
    PSH.species <- list()
    for(i in 1:length(PSHMatrix2011)){
      PSH.species[[i]] <- data.frame("2011"=PSHMatrix2011.15[[i]], "2012"=PSHMatrix2012.15[[i]], "2013"=PSHMatrix2013.15[[i]], 
                                     "2014"=PSHMatrix2014.15[[i]], "2015"=PSHMatrix2015.15[[i]], "2016"=PSHMatrix2016.15[[i]])[1:30,]
    }
    names(PSH.species) <- paste("PSH", names(PSHMatrix2011), sep=".")

    # Extract all YAN species
    YAN.species <- list()
    for(i in 1:length(YANMatrix2011)){
      YAN.species[[i]] <- data.frame("2011"=YANMatrix2011.15[[i]], "2012"=YANMatrix2012.15[[i]], "2013"=YANMatrix2013.15[[i]], 
                                     "2014"=YANMatrix2014.15[[i]], "2015"=YANMatrix2015.15[[i]])
    }
    names(YAN.species) <- paste("YAN", names(YANMatrix2011), sep=".")

    # Extract all NAK species
    NAK.species <- list()
    for(i in 1:length(NAKMatrix2011)){
      NAK.species[[i]] <- data.frame("2010"=NAKMatrix2010.15[[i]], "2011"=NAKMatrix2011.15[[i]], "2012"=NAKMatrix2012.15[[i]], 
                                     "2013"=NAKMatrix2013.15[[i]], "2014"=NAKMatrix2014.15[[i]])
                                     #"2015"=NAKMatrix2015.15[[i]])
    }
    names(NAK.species) <- paste("NAK", names(NAKMatrix2011), sep=".")

    # Extract all RNF species
    RNF.species <- list()
    for(i in 1:length(RNFMatrix2011)){
      RNF.species[[i]] <- data.frame("2010"=RNFMatrix2010.15[[i]], "2011"=RNFMatrix2011.15[[i]], "2012"=RNFMatrix2012.15[[i]], 
                                     "2013"=RNFMatrix2013.15[[i]], "2014"=RNFMatrix2014.15[[i]], "2015"=RNFMatrix2015.15[[i]])
                                     #"2016"=RNFMatrix2016.15[[i]])
    }
    names(RNF.species) <- paste("RNF", names(RNFMatrix2011), sep=".")


# Remove all populations of "Dendrohyrax arboreus", "Tragulus javanicus", "Tragulus napu" and "Muntiacus muntjak"

#UDZ.species <- UDZ.species[-10] # Remove Dendrohyrax arboreus from UDZ
PSH.species <- PSH.species[-33]  # Remove Tragulus napu from PSH
PSH.species <- PSH.species[-17]  # Remove Muntiacus muntjak from PSH
NAK.species <- NAK.species[-22]  # Remove Tragulus javanicus from NAK
NAK.species <- NAK.species[-16] # Remove Muntiacus muntjak from NAK

All_species7sites <- c(VB_.species=VB_.species,
                          UDZ.species=UDZ.species,
                          BIF.species=BIF.species,
                          PSH.species=PSH.species,
                          YAN.species=YAN.species,
                          NAK.species=NAK.species,
                          RNF.species=RNF.species)

save(All_species7sites, file="All_species7sites.RData")

###################### EXTRACT 32 POPULATIONS MODELED WITH COVARIATES IN WPI ANALYSIS ########################
# Covariate species list available in file "WPI_Covariate_Populations.xlsx"; Can be recreated using "WPI_Analysis.R" file

# 5 VB species
#  names(VBMatrix2008)
#  VB_.Pecari_tajacu <- VB_.species$"VB_.Pecari tajacu"
#  VB_.Dasyprocta_punctata <- VB_.species$"VB_.Dasyprocta punctata"
#  VB_.Tapirus_bairdii <- VB_.species$"VB_.Tapirus bairdii" 
#  VB_.Cuniculus_paca <-  VB_.species$"VB_.Cuniculus paca" 
#  VB_.Dasypus_novemcinctus <- VB_.species$"VB_.Dasypus novemcinctus"

#  VB_covariate_species <- list(VB_.Pecari_tajacu=VB_.Pecari_tajacu,
#                             VB_.Dasyprocta_punctata=VB_.Dasyprocta_punctata,
#                             VB_.Tapirus_bairdii=VB_.Tapirus_bairdii,
#                             VB_.Cuniculus_paca=VB_.Cuniculus_paca,
#                             VB_.Dasypus_novemcinctus=VB_.Dasypus_novemcinctus)

# 10 UDZ species
#  names(UDZMatrix2009)
#  UDZ.Guttera_pucherani <- UDZ.species$"UDZ.Guttera pucherani" 
#  UDZ.Bdeogale_crassicauda <- UDZ.species$"UDZ.Bdeogale crassicauda"
#  UDZ.Cricetomys_gambianus <- UDZ.species$"UDZ.Cricetomys gambianus"
#  UDZ.Cephalophus_spadix <- UDZ.species$"UDZ.Cephalophus spadix"
#  UDZ.Genetta_servalina <- UDZ.species$"UDZ.Genetta servalina"
#  UDZ.Cephalophus_harveyi <- UDZ.species$"UDZ.Cephalophus harveyi"
#  UDZ.Paraxerus_vexillarius <- UDZ.species$"UDZ.Paraxerus vexillarius"
#  UDZ.Cercocebus_sanjei <- UDZ.species$"UDZ.Cercocebus sanjei"
#  UDZ.Rhynchocyon_udzungwensis <- UDZ.species$"UDZ.Rhynchocyon udzungwensis"
  #UDZ.Nesotragus_moschatus <- #NB mismatch of names causing species to be absent from matched list

#  UDZ_covariate_species <- list(UDZ.Guttera_pucherani=UDZ.Guttera_pucherani, 
#                              UDZ.Bdeogale_crassicauda=UDZ.Bdeogale_crassicauda, 
#                             UDZ.Cricetomys_gambianus=UDZ.Cricetomys_gambianus, 
#                              UDZ.Cephalophus_spadix=UDZ.Cephalophus_spadix, 
#                              UDZ.Genetta_servalina=UDZ.Genetta_servalina, 
#                              UDZ.Cephalophus_harveyi=UDZ.Cephalophus_harveyi, 
#                              UDZ.Paraxerus_vexillarius=UDZ.Paraxerus_vexillarius, 
#                              UDZ.Cercocebus_sanjei=UDZ.Cercocebus_sanjei,
#                              #UDZ.Nesotragus_moschatus=UDZ.Nesotragus_moschatus,
#                              UDZ.Rhynchocyon_udzungwensis=UDZ.Rhynchocyon_udzungwensis)


# 3 BIF species
#  names(BIFMatrix2011)
#  BIF.Cercopithecus_lhoesti <- BIF.species$"BIF.Cercopithecus lhoesti"
#  BIF.Cephalophus_silvicultor <- BIF.species$"BIF.Cephalophus silvicultor"
#  BIF.Cephalophus_nigrifrons <- BIF.species$"BIF.Cephalophus nigrifrons"

#  BIF_covariate_species <- list(BIF.Cercopithecus_lhoesti=BIF.Cercopithecus_lhoesti,
#                              BIF.Cephalophus_silvicultor=BIF.Cephalophus_silvicultor,
#                              BIF.Cephalophus_nigrifrons=BIF.Cephalophus_nigrifrons)

# 5 PSH species
#  names(PSHMatrix2011)
#  PSH.Leopoldamys_sabanus <- PSH.species$"PSH.Leopoldamys sabanus" 
#  PSH.Muntiacus_muntjak <- PSH.species$"PSH.Muntiacus muntjak"
#  PSH.Macaca_nemestrina <- PSH.species$"PSH.Macaca nemestrina"
#  PSH.Sus_scrofa <- PSH.species$"PSH.Sus scrofa" 
#  PSH.Tragulus_kanchil <- PSH.species$"PSH.Tragulus kanchil"

#  PSH_covariate_species <- list(PSH.Leopoldamys_sabanus=PSH.Leopoldamys_sabanus,
#                              PSH.Muntiacus_muntjak=PSH.Muntiacus_muntjak,
#                              PSH.Macaca_nemestrina=PSH.Macaca_nemestrina,
#                              PSH.Sus_scrofa=PSH.Sus_scrofa,
#                              PSH.Tragulus_kanchil=PSH.Tragulus_kanchil)

# 5 YAN species
#  names(YANMatrix2011)
#  YAN.Cuniculus_paca <- YAN.species$"YAN.Cuniculus paca"
#  YAN.Dasyprocta_fuliginosa <- YAN.species$"YAN.Dasyprocta fuliginosa"
#  YAN.Dasypus_novemcinctus <- YAN.species$"YAN.Dasypus novemcinctus"
#  YAN.Mitu_tuberosum <- YAN.species$"YAN.Mitu tuberosum"
#  YAN.Tapirus_terrestris <- YAN.species$"YAN.Tapirus terrestris"

#  YAN_covariate_species <- list(YAN.Cuniculus_paca=YAN.Cuniculus_paca,
#                              YAN.Dasyprocta_fuliginosa=YAN.Dasyprocta_fuliginosa,
#                              YAN.Dasypus_novemcinctus=YAN.Dasypus_novemcinctus,
#                              YAN.Mitu_tuberosum=YAN.Mitu_tuberosum,
#                              YAN.Tapirus_terrestris=YAN.Tapirus_terrestris)

# 2 NAK species
#  names(NAKMatrix2011)
#  NAK.Muntiacus_muntjak <- NAK.species$"NAK.Muntiacus muntjak"
#  NAK.Atherurus_macrourus <- NAK.species$"NAK.Atherurus macrourus"

#  NAK_covariate_species <- list(NAK.Muntiacus_muntjak=NAK.Muntiacus_muntjak,
#                              NAK.Atherurus_macrourus=NAK.Atherurus_macrourus)


# 2 RNF species
#  names(RNFMatrix2011)
#  RNF.Nesomys_rufus <- RNF.species$"RNF.Nesomys rufus"
#  RNF.Fossa_fossana <- RNF.species$"RNF.Fossa fossana"

#  RNF_covariate_species <- list(RNF.Nesomys_rufus=RNF.Nesomys_rufus,
#                              RNF.Fossa_fossana=RNF.Fossa_fossana)


################## COMBINE ALL 32 SPECIES MATRICES INTO A SINGLE LIST ########################
#All500m_covariate_species <- c(VB_covariate_species=VB_covariate_species,
#                                  UDZ_covariate_species=UDZ_covariate_species, 
#                                  BIF_covariate_species=BIF_covariate_species,
#                                  PSH_covariate_species=PSH_covariate_species,
#                                  YAN_covariate_species=YAN_covariate_species,
#                                  NAK_covariate_species=NAK_covariate_species,
#                                  RNF_covariate_species=RNF_covariate_species)

#All500m_covariate_species <- list(VB_.Pecari_tajacu=VB_.Pecari_tajacu,
#                              VB_.Dasyprocta_punctata=VB_.Dasyprocta_punctata,
#                              VB_.Tapirus_bairdii=VB_.Tapirus_bairdii,
#                              VB_.Cuniculus_paca=VB_.Cuniculus_paca,
#                              VB_.Dasypus_novemcinctus=VB_.Dasypus_novemcinctus,
#                              UDZ.Guttera_pucherani=UDZ.Guttera_pucherani, 
#                              UDZ.Bdeogale_crassicauda=UDZ.Bdeogale_crassicauda, 
#                              UDZ.Cricetomys_gambianus=UDZ.Cricetomys_gambianus, 
#                              UDZ.Cephalophus_spadix=UDZ.Cephalophus_spadix, 
#                              UDZ.Genetta_servalina=UDZ.Genetta_servalina, 
#                              UDZ.Cephalophus_harveyi=UDZ.Cephalophus_harveyi, 
#                              UDZ.Paraxerus_vexillarius=UDZ.Paraxerus_vexillarius, 
#                              UDZ.Cercocebus_sanjei=UDZ.Cercocebus_sanjei,
#                              UDZ.Rhynchocyon_udzungwensis=UDZ.Rhynchocyon_udzungwensis, 
#                              UDZ.Nesotragus_moschatus=UDZ.Nesotragus_moschatus,
#                              BIF.Cercopithecus_lhoesti=BIF.Cercopithecus_lhoesti,
#                              BIF.Cephalophus_silvicultor=BIF.Cephalophus_silvicultor,
#                              BIF.Cephalophus_nigrifrons=BIF.Cephalophus_nigrifrons,
#                              PSH.Leopoldamys_sabanus=PSH.Leopoldamys_sabanus,
#                              PSH.Muntiacus_muntjak=PSH.Muntiacus_muntjak,
#                              PSH.Macaca_nemestrina=PSH.Macaca_nemestrina,
#                              PSH.Sus_scrofa=PSH.Sus_scrofa,
#                              PSH.Tragulus_kanchil=PSH.Tragulus_kanchil,
#                              YAN.Cuniculus_paca=YAN.Cuniculus_paca,
#                              YAN.Dasyprocta_fuliginosa=YAN.Dasyprocta_fuliginosa,
#                              YAN.Dasypus_novemcinctus=YAN.Dasypus_novemcinctus,
#                              YAN.Mitu_tuberosum=YAN.Mitu_tuberosum,
#                              YAN.Tapirus_terrestris=YAN.Tapirus_terrestris,
#                              NAK.Muntiacus_muntjak=NAK.Muntiacus_muntjak,
#                              NAK.Atherurus_macrourus=NAK.Atherurus_macrourus,
#                              RNF.Nesomys_rufus=RNF.Nesomys_rufus,
#                              RNF.Fossa_fossana=RNF.Fossa_fossana)


#save(All500m_covariate_species, file="All500m_covariate_species.RData")


########### Format presence absence matrices for sites
# Reduce overall data to the 7 sites only
Sites7dataDF <- eventsdata[eventsdata$Site.Code=="VB-"|eventsdata$Site.Code=="UDZ"|eventsdata$Site.Code=="BIF"|eventsdata$Site.Code=="PSH"|eventsdata$Site.Code=="YAN"|eventsdata$Site.Code=="NAK"|eventsdata$Site.Code=="RNF",]
spnames <- c(names(VBMatrix2008), names(UDZMatrix2009), names(BIFMatrix2011), names(PSHMatrix2011), names(YANMatrix2011), names(NAKMatrix2011), names(RNFMatrix2011))
#spnames[order(spnames)]
#Sites7dataDF[order(spnames)]


SitesdataVB <- eventsdata[eventsdata$Site.Code=="VB-",]
SitesdataUDZ <- eventsdata[eventsdata$Site.Code=="UDZ",]
SitesdataBIF <- eventsdata[eventsdata$Site.Code=="BIF",]
SitesdataPSH <- eventsdata[eventsdata$Site.Code=="PSH",]
SitesdataYAN <- eventsdata[eventsdata$Site.Code=="YAN",]
SitesdataNAK <- eventsdata[eventsdata$Site.Code=="NAK",]
SitesdataRNF <- eventsdata[eventsdata$Site.Code=="RNF",]


Sites7data <- list(SitesdataVB, SitesdataUDZ, SitesdataBIF, SitesdataPSH, SitesdataYAN, SitesdataNAK, SitesdataRNF)
#Sitenames <- list("SitesdataVB", "SitesdataUDZ", "SitesdataBIF", "SitesdataPSH", "SitesdataYAN", "SitesdataNAK", "SitesdataRNF")
Sitenames <- list("VB", "UDZ", "BIF", "PSH", "YAN", "NAK", "RNF")


spnames <- list(names(VBMatrix2008), names(UDZMatrix2009), names(BIFMatrix2011), names(PSHMatrix2011), names(YANMatrix2011), names(NAKMatrix2011), names(RNFMatrix2011))


# CREATE A UNIQUE IDENTIFICATION TO INDEX ALL 126 SPECIES
# Make list of unique species names for the 7 sites
# Add Muntiacus montanus here to correspond with UID list from Montreal
UID <- unique(c(names(VBMatrix2008), names(UDZMatrix2009), names(BIFMatrix2011), 
                names(PSHMatrix2011), names(YANMatrix2011), names(NAKMatrix2011), 
                names(RNFMatrix2011), "Muntiacus montanus"))# Remove 

UID <- UID[order(UID)]

#### NEW FOR 7/30/2015 Data - need to remove Atelocynus microtis from data used in analysis and from UID; New species detection in YAN

# Remove 4 extra species "Dendrohyrax arboreus", "Tragulus javanicus", "Tragulus napu" and "Muntiacus muntjak"
UID <- UID[-132]
UID <- UID[-130]
UID <- UID[-75]
UID <- data.frame(UID, 1:length(UID))
colnames(UID) <- c("bin", "UID")

#UID$UID[match(as.factor(rownames(SitesBinary[[1]])), UID$bin)]

# Sort by rownames then add table to list and save as an RData object
# 
SitesBinary <- list()
SpID <- vector()

for(i in 1:length(Sites7data)){
  sptable <- table(Sites7data[[i]]$bin, Sites7data[[i]]$Sampling.Unit.Name)
  sptable <- sptable[match(spnames[[i]], rownames(sptable)),]
  sptable <- ifelse(sptable>0,1,0)
  sptable <- sptable[order(rownames(sptable)),]
  SpID    <- UID$UID[match(as.factor(rownames(sptable)), UID$bin)]
  sptable <- cbind(SpID, sptable)
  SitesBinary[[i]] <- sptable
  #outputname <- paste(Sitenames[i], "Binary", "csv", sep=".")
  #write.csv(sptable, file=outputname)
}

# Remove PSH camera traps from arrays 2 and 3
SitesBinary[[4]] <- SitesBinary[[4]][,1:31]

# Remove all populations of "Dendrohyrax arboreus", "Tragulus javanicus", "Tragulus napu" and "Muntiacus muntjak"
#SitesBinary[[2]] <- SitesBinary[[2]][-10,] # Remove Dendrohyrax arboreus from UDZ
SitesBinary[[4]] <- SitesBinary[[4]][-33,] # Remove Tragulus napu from PSH
SitesBinary[[4]] <- SitesBinary[[4]][-17,] # Remove Muntiacus muntjak from PSH
SitesBinary[[6]] <- SitesBinary[[6]][-22,] # Remove Tragulus javanicus from NAK
SitesBinary[[6]] <- SitesBinary[[6]][-16,] # Remove Muntiacus muntjak from NAK

names(SitesBinary) <- Sitenames
save(SitesBinary, file="SitesBinary.RData")

# Create Presence-Absence matrices for each CT for each year

temp <- list()
SitesBinaryAnnual <- list()
SpID <- vector()

SiteYears <- length(table(Sites7data[[1]]$Sampling.Period)) + length(table(Sites7data[[2]]$Sampling.Period)) +
             length(table(Sites7data[[3]]$Sampling.Period)) + length(table(Sites7data[[4]]$Sampling.Period)) + 
             length(table(Sites7data[[5]]$Sampling.Period)) + length(table(Sites7data[[6]]$Sampling.Period)) + 
             length(table(Sites7data[[7]]$Sampling.Period))  

# First get inner loop to work (i.e. get it to make a separate matrix for each year at one site)
# Then add outer loop to loop over all sites

for(i in 1:length(Sites7data)){
  for(j in 1:length(table(Sites7data[[i]]$Sampling.Period))){
  sptable <- table(Sites7data[[i]]$bin, Sites7data[[i]]$Sampling.Unit.Name, Sites7data[[i]]$Sampling.Period)[,,j]
  sptable <- sptable[match(spnames[[i]], rownames(sptable)),]
  sptable <- ifelse(sptable>0,1,0)
  sptable <- sptable[order(rownames(sptable)),]
  SpID    <- UID$UID[match(as.factor(rownames(sptable)), UID$bin)]
  sptable <- cbind(SpID, sptable)
  temp[[j]] <- sptable
  }
  names(temp) <- names(table(Sites7data[[i]]$Sampling.Period))
  SitesBinaryAnnual[[i]] <- temp
  rm(temp)
  temp <- list()
  #outputname <- paste(Sitenames[i], "Binary", "csv", sep=".")
  #write.csv(sptable, file=outputname)
}
names(SitesBinaryAnnual) <- Sitenames

SitesBinaryAnnual[[4]] <- lapply(SitesBinaryAnnual[[4]], "[", ,1:31) # Remove arrays 2 and 3 at PSH
#SitesBinaryAnnual[[2]] <- lapply(SitesBinaryAnnual[[2]], "[", -10,)  # Remove Dendrohyrax arboreus from UDZ
SitesBinaryAnnual[[4]] <- lapply(SitesBinaryAnnual[[4]], "[", -33,)  # Remove Tragulus napu from PSH
SitesBinaryAnnual[[4]] <- lapply(SitesBinaryAnnual[[4]], "[", -17,)  # Remove Muntiacus muntjak from PSH
SitesBinaryAnnual[[6]] <- lapply(SitesBinaryAnnual[[6]], "[", -22,)  # Remove Tragulus javanicus from NAK
SitesBinaryAnnual[[6]] <- lapply(SitesBinaryAnnual[[6]], "[", -16,)  # Remove Muntiacus muntjak from NAK

save(SitesBinaryAnnual, file="SitesBinaryAnnual.RData")



# Figure out why there are 4 mismatched species between what I sent JP for phylogeny and overall matrix

# WORKS BELOW for overall data - no longer works b/c missing spnames to match on; try using UID instead
# Table species by camera trap
#sptable <- table(Sites7data$bin, Sites7data$Sampling.Unit.Name)

# Reduce the table to only the species for which we are interested
#sptable <- sptable[match(spnames, rownames(sptable)),]
#sptable <- ifelse(sptable>0,1,0)

#write.csv(sptable, file="Species_PresenceAbsenceTable.csv")
#rowSums(sptable)

########### EXAMINE THE OCCUPANCY, # COLONIZATION and # EXTINCTION EVENTS FOR EACH SPECIES AT EACH SITE ##############
# We want the number of times that a species at a camera trap changes from 0 to 1; then also from 1 to 0
# These objects give the transition between year t and t+1



lapply(SitesBinaryAnnual[[4]], "rowSums") - lapply(SitesBinaryAnnual[[4]], "[", ,1)

# Gives the number of camera traps at which a species is present in a year (not generalized)
data.frame(lapply(SitesBinaryAnnual[[4]], "rowSums")[1]) - data.frame(lapply(SitesBinaryAnnual[[4]], "[", ,1)[1])


##################### This gives the number of col and ext events for site 4 for the first transition
T1 <- data.frame(SitesBinaryAnnual[[4]][2]) - data.frame(SitesBinaryAnnual[[4]][1])
T1col <- T1
T1col[T1col==-1] <- 0
T1col <- rowSums(T1col)

T1ext <- T1
T1ext[T1ext==1] <- 0
T1ext <- rowSums(T1ext)

T1.ID <- data.frame(SitesBinaryAnnual[[4]][2])[,1] # reassign species ID values 
data.frame(UID=T1.ID, T1col, T1ext)

####################### End code for the col and ext events for site 4 for the first transition

##################### Generalize number of col and ext events for all sites and for all annual transitions; include occupancy (# CTs with species present in a primary period)  
Ti <- list()
Col <- list()
Ext <- list()
Occ <- list()
Tsum <- list()
#Species <- data.frame(SitesBinaryAnnual[[4]][2])[,1]
Species <- data.frame()
Percent <- data.frame()
OccColExt.Raw <- list()
OccColExt.Per <- list()

for(j in 1:length(SitesBinaryAnnual)){
    Species <- data.frame(SitesBinaryAnnual[[j]][2])[,1]
for(i in 1:(length(SitesBinaryAnnual[[j]])-1)){
    Ti[[i]] <- data.frame(SitesBinaryAnnual[[j]][i+1]) - data.frame(SitesBinaryAnnual[[j]][i])
    
    Occ[[i]] <- data.frame(lapply(SitesBinaryAnnual[[j]], "rowSums")[i+1]) - data.frame(lapply(SitesBinaryAnnual[[j]], "[", ,1)[1]) # Note initial year of occupancy data is omitted
    
    Col[[i]] <- Ti[[i]]
    Col[[i]][Col[[i]]==-1] <- 0
    Col[[i]] <- rowSums(Col[[i]])

    Ext[[i]] <- Ti[[i]]
    Ext[[i]][Ext[[i]]==1] <- 0
    Ext[[i]] <- rowSums(Ext[[i]])

    Tsum[[i]] <- data.frame(Col[[i]], Ext[[i]])
    
    Species <- data.frame(Species, Occ[[i]], Col[[i]], Ext[[i]])
    Percent <- data.frame(Species=Species[,1], round((Species[,2:dim(Species)[2]]/(dim(Ti[[1]])[2]-1))*100,1))
  }

  OccColExt.Raw[[j]] <- Species
  OccColExt.Per[[j]] <- Percent
  #rm(temp)
  #temp <- list()
}

names(OccColExt.Raw) <- Sitenames
names(OccColExt.Per) <- Sitenames

# above loop works and creates a list of occupancy, col, ext events for each time period transition for all sites

save(OccColExt.Raw, file="OccColExt.Raw.RData")
save(OccColExt.Per, file="OccColExt.Per.RData")





####################### End code for the col and ext events 






################## CAMERA TRAP TEMPERATURE DATA FORMATTING ####################

# Clean temperature data and convert F measurements to C
  library(stringr)
  eventsdata$Temperature <- str_trim(eventsdata$Temperature, side="both")
  temp.degrees <- sub(" .*","", eventsdata$Temperature)
  temp.unit <- str_sub(eventsdata$Temperature, nchar(eventsdata$Temperature), nchar(eventsdata$Temperature))
  temp.unit <- str_trim(temp.unit, side="both")
  #Convert F temperature values to Celcius; Change C Temperatures to numeric format from character
  temp.degreesC <- as.numeric(ifelse(temp.unit=="F", f.FtoC(as.numeric(temp.degrees)), temp.degrees))
  eventsdata <- cbind(eventsdata, temp.degreesC)

# Note that CT-PSH-1-21 has problematic temperature values (Max=67 degrees) 
# Change all temp values to NA for CT-PSH-1-21
#  eventsdata$temp.degreesC <- ifelse(eventsdata$Sampling.Unit.Name=="CT-PSH-1-21", NA, eventsdata$temp.degreesC)
# Note that by changing "CT-PSH-1-21" values to NA, the camera trap is omitted from aggregated values, which is problematic downstream.
# Try changing value downstream instead of ahead of aggregation.

# Determine the annual min, max and variance of the non-calibrated temperature data for each CT without using interpolated data
  Temp.Min <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=min)
  names(Temp.Min) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Min")
  Temp.Max <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=max)
  names(Temp.Max) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Max")
  Temp.Var <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=var)
  names(Temp.Var) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Var")
  Temp.Mean <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=mean)
  names(Temp.Mean) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Mean")
  Temp.SD <- aggregate(eventsdata$temp.degreesC ~ eventsdata$Site.Code + eventsdata$Sampling.Unit.Name + eventsdata$Sampling.Period, FUN=sd)
  names(Temp.SD) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.SD")
  CT.Year <- as.integer(substr(Temp.SD$Sampling.Period,1,4))
  CT.Temp <- cbind(Temp.Min, Temp.Max$Temp.Max, Temp.Var$Temp.Var, Temp.Mean$Temp.Mean, Temp.SD$Temp.SD, Year=CT.Year)
  names(CT.Temp) <- c("Site.Code", "Sampling.Unit.Name", "Sampling.Period", "Temp.Min", "Temp.Max", "Temp.Var", "Temp.Mean", "Temp.SD", "Year")

# Note that CT-PSH-1-21 has problematic temperature values (Max=67 degrees) 
# Change all temp values to NA for CT-PSH-1-21
CT.Temp$Temp.Min <- ifelse(CT.Temp$Sampling.Unit.Name=="CT-PSH-1-21", NA, CT.Temp$Temp.Min)
CT.Temp$Temp.Max <- ifelse(CT.Temp$Sampling.Unit.Name=="CT-PSH-1-21", NA, CT.Temp$Temp.Max)
CT.Temp$Temp.Var <- ifelse(CT.Temp$Sampling.Unit.Name=="CT-PSH-1-21", NA, CT.Temp$Temp.Var)
CT.Temp$Temp.Mean <- ifelse(CT.Temp$Sampling.Unit.Name=="CT-PSH-1-21", NA, CT.Temp$Temp.Mean)
CT.Temp$Temp.SD <- ifelse(CT.Temp$Sampling.Unit.Name=="CT-PSH-1-21", NA, CT.Temp$Temp.SD)








########### NB: THE FOLLOW TWO FILES FOR ELEV AND FOREST LOSS WILL NEED TO BE UPDATED WITH NEW FILES FROM ALEX ################
# Bring in new elevation data from Alex Zvoleff as of 10/21/2014
  load("ct_pts_elev.RData")

# Create object with temperature covariate data to use for all sites with 500 m elevation gradients 

  Alltemp500 <- data.frame(Site.Code=CT.Temp$Site.Code, 
                         Sampling.Unit.Name=CT.Temp$Sampling.Unit.Name,
                         Year=as.factor(CT.Temp$Year),
                         Temp.Min=CT.Temp$Temp.Min, 
                         Temp.Max=CT.Temp$Temp.Max,
                         Temp.Var=CT.Temp$Temp.Var, 
                         Temp.SD=CT.Temp$Temp.SD,
                         Temp.Mean=CT.Temp$Temp.Mean)

######################## FORMAT OTHER COVARIATE DATA SOURCES (i.e. Elevation, forest loss) ###################
# Read in elevation data
  #ELEV <- read.csv("CT_edgedist_elevation_final.txt") # Elevation data from Melissa Rosa to be replaced with Data from Alex Zvoleff
  ELEV <- ct_pts_elev
  

# Read in CT specific forest loss data from Alex (spans 2000-2013)
  traps_fc <- read.csv("traps_fc.csv")
  #load("ct_fc_loss.RData")
  #traps_fc <- ct_loss  

  ftraps120 <- traps_fc[traps_fc$buffer_m=="120m buffer",]
  ftraps120 <- cbind(ftraps120, ELEV[match(ftraps120$trap_ID, ELEV$ct_ID),])
  ftraps120_500m <- ftraps120[ftraps120$sitecode=="VB"|
                                ftraps120$sitecode=="UDZ"|
                                ftraps120$sitecode=="BIF"|
                                ftraps120$sitecode=="PSH"|
                                ftraps120$sitecode=="YAN"|
                                ftraps120$sitecode=="NAK"|
                                ftraps120$sitecode=="RNF",]

  ELEV_FL <- data.frame(Site.Code=ftraps120_500m$sitecode, 
                             Sampling.Unit.Name=ftraps120_500m$ct_ID, 
                             Elevation=ftraps120_500m$elev_m,
                             FL120=ftraps120_500m$fc_frac_loss, 
                             FG120=ftraps120_500m$fc_frac_gain)
# Remove inactive camera traps from YAN to preseve integrity of covariate structure below  
  ELEV_FL <- ELEV_FL[-437,] # Remove CT-YAN-1-04
  ELEV_FL <- ELEV_FL[-395,] # Remove CT-YAN-2-16   


# Read in site level forest loss calculations from Alex (spans 5 years prior to CT sampling start at each site)
# NB manually update VB sitecode to "VB-" in csv file to enable merging
  #site_fc <- read.csv("20141004_forest_loss.csv")
# Extract site level forest loss for protected areas
  #FL_PA <- site_fc[site_fc$aoi=="PA",]


######################## COMBINE TEMPERATURE COVARIATE DATA WITH OTHER COVARIATE DATA SOURCES (i.e. Elevation, forest loss) ###################


# Create site specific covariate lists with scaled covariates
  CT.Temp.VB <- melt(Alltemp500[Alltemp500$Site.Code=="VB-",])
  ELEV_FL.VB <- melt(ELEV_FL[ELEV_FL$Site.Code=="VB",])
    VB.Elev <- cast(ELEV_FL.VB, Sampling.Unit.Name ~ variable)[,2]
    VB.FL120 <- cast(ELEV_FL.VB, Sampling.Unit.Name ~ variable)[,3]
    VB.FG120 <- cast(ELEV_FL.VB, Sampling.Unit.Name ~ variable)[,4] 
    VB.Tmin <- as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Year ~ variable)[,,1])
    VB.Tmax <- as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Year ~ variable)[,,2])
    VB.Tvar <- round(as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Year ~ variable)[,,3]),2)
    VB.SD <- round(as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Year ~ variable)[,,4]),2)
    VB.Mean <- as.data.frame(cast(CT.Temp.VB, Sampling.Unit.Name ~ Year ~ variable)[,,5])


      VB.Elev <- scale(VB.Elev)
      VB.Elev[is.na(VB.Elev)] <- 0

      VB.FL120 <- scale(VB.FL120)
      VB.FL120[is.na(VB.FL120)] <- 0

      VB.FG120 <- scale(VB.FG120)
      VB.FG120[is.na(VB.FG120)] <- 0

      VB.Tmin <- scale(VB.Tmin)
      VB.Tmin[is.na(VB.Tmin)] <- 0

      VB.Tmax <- scale(VB.Tmax)
      VB.Tmax[is.na(VB.Tmax)] <- 0

      VB.Tvar <- scale(VB.Tvar)
      VB.Tvar[is.na(VB.Tvar)] <- 0

      VB.SD <- scale(VB.SD)
      VB.SD[is.na(VB.SD)] <- 0

      VB.Mean <- scale(VB.Mean)
      VB.Mean[is.na(VB.Mean)] <- 0


          VB__covs <- list(Elevation=VB.Elev,
                ForestLossCT=VB.FL120,
                ForestGainCT=VB.FG120,
                Tmin=VB.Tmin,
                Tmax=VB.Tmax,
                Tvar=VB.Tvar,
                Tsd=VB.SD,
                Tmean=VB.Mean)

    CT.Temp.UDZ <- melt(Alltemp500[Alltemp500$Site.Code=="UDZ",])
  ELEV_FL.UDZ <- melt(ELEV_FL[ELEV_FL$Site.Code=="UDZ",])
    UDZ.Elev <- cast(ELEV_FL.UDZ, Sampling.Unit.Name ~ variable)[,2]
    UDZ.FL120 <- cast(ELEV_FL.UDZ, Sampling.Unit.Name ~ variable)[,3]
    UDZ.FG120 <- cast(ELEV_FL.UDZ, Sampling.Unit.Name ~ variable)[,4] 
    UDZ.Tmin <- as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Year ~ variable)[,,1])
    UDZ.Tmax <- as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Year ~ variable)[,,2])
    UDZ.Tvar <- round(as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Year ~ variable)[,,3]),2)
    UDZ.SD <- round(as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Year ~ variable)[,,4]),2)
    UDZ.Mean <- as.data.frame(cast(CT.Temp.UDZ, Sampling.Unit.Name ~ Year ~ variable)[,,5])


      UDZ.Elev <- scale(UDZ.Elev)
      UDZ.Elev[is.na(UDZ.Elev)] <- 0

      UDZ.FL120 <- scale(UDZ.FL120)
      UDZ.FL120[is.na(UDZ.FL120)] <- 0

      UDZ.FG120 <- scale(UDZ.FG120)
      UDZ.FG120[is.na(UDZ.FG120)] <- 0

      UDZ.Tmin <- scale(UDZ.Tmin)
      UDZ.Tmin[is.na(UDZ.Tmin)] <- 0

      UDZ.Tmax <- scale(UDZ.Tmax)
      UDZ.Tmax[is.na(UDZ.Tmax)] <- 0

      UDZ.Tvar <- scale(UDZ.Tvar)
      UDZ.Tvar[is.na(UDZ.Tvar)] <- 0

      UDZ.SD <- scale(UDZ.SD)
      UDZ.SD[is.na(UDZ.SD)] <- 0

      UDZ.Mean <- scale(UDZ.Mean)
      UDZ.Mean[is.na(UDZ.Mean)] <- 0


          UDZ_covs <- list(Elevation=UDZ.Elev,
                ForestLossCT=UDZ.FL120,
                ForestGainCT=UDZ.FG120,
                Tmin=UDZ.Tmin,
                Tmax=UDZ.Tmax,
                Tvar=UDZ.Tvar,
                Tsd=UDZ.SD,
                Tmean=UDZ.Mean)




  CT.Temp.BIF <- melt(Alltemp500[Alltemp500$Site.Code=="BIF",])
  ELEV_FL.BIF <- melt(ELEV_FL[ELEV_FL$Site.Code=="BIF",])
    BIF.Elev <- cast(ELEV_FL.BIF, Sampling.Unit.Name ~ variable)[,2]
    BIF.FL120 <- cast(ELEV_FL.BIF, Sampling.Unit.Name ~ variable)[,3]
    BIF.FG120 <- cast(ELEV_FL.BIF, Sampling.Unit.Name ~ variable)[,4] 
    BIF.Tmin <- as.data.frame(cast(CT.Temp.BIF, Sampling.Unit.Name ~ Year ~ variable)[,,1])
    BIF.Tmax <- as.data.frame(cast(CT.Temp.BIF, Sampling.Unit.Name ~ Year ~ variable)[,,2])
    BIF.Tvar <- round(as.data.frame(cast(CT.Temp.BIF, Sampling.Unit.Name ~ Year ~ variable)[,,3]),2)
    BIF.SD <- round(as.data.frame(cast(CT.Temp.BIF, Sampling.Unit.Name ~ Year ~ variable)[,,4]),2)
    BIF.Mean <- as.data.frame(cast(CT.Temp.BIF, Sampling.Unit.Name ~ Year ~ variable)[,,5])


      BIF.Elev <- scale(BIF.Elev)
      BIF.Elev[is.na(BIF.Elev)] <- 0

      BIF.FL120 <- scale(BIF.FL120)
      BIF.FL120[is.na(BIF.FL120)] <- 0

      BIF.FG120 <- scale(BIF.FG120)
      BIF.FG120[is.na(BIF.FG120)] <- 0

      BIF.Tmin <- BIF.Tmin[,1:6] # SUBSET SO THAT # YEARS MATCHES FOR SP & TEMP DATA
      BIF.Tmin <- scale(BIF.Tmin)
      BIF.Tmin[is.na(BIF.Tmin)] <- 0

      BIF.Tmax <- BIF.Tmax[,1:6]
      BIF.Tmax <- scale(BIF.Tmax)
      BIF.Tmax[is.na(BIF.Tmax)] <- 0

      BIF.Tvar <- BIF.Tvar[,1:6]
      BIF.Tvar <- scale(BIF.Tvar)
      BIF.Tvar[is.na(BIF.Tvar)] <- 0

      BIF.SD <- BIF.SD[,1:6]
      BIF.SD <- scale(BIF.SD)
      BIF.SD[is.na(BIF.SD)] <- 0

      BIF.Mean <- BIF.Mean[,1:6]
      BIF.Mean <- scale(BIF.Mean)
      BIF.Mean[is.na(BIF.Mean)] <- 0


          BIF_covs <- list(Elevation=BIF.Elev,
                ForestLossCT=BIF.FL120,
                ForestGainCT=BIF.FG120,
                Tmin=BIF.Tmin,
                Tmax=BIF.Tmax,
                Tvar=BIF.Tvar,
                Tsd=BIF.SD,
                Tmean=BIF.Mean)


## All Pasoh data (camera trap and covariates) need to be limited to Array 1 only
  CT.Temp.PSH <- melt(Alltemp500[Alltemp500$Site.Code=="PSH",])
  ELEV_FL.PSH <- melt(ELEV_FL[ELEV_FL$Site.Code=="PSH",])
    PSH.Elev <- cast(ELEV_FL.PSH, Sampling.Unit.Name ~ variable)[,2]
    PSH.FL120 <- cast(ELEV_FL.PSH, Sampling.Unit.Name ~ variable)[,3]
    PSH.FG120 <- cast(ELEV_FL.PSH, Sampling.Unit.Name ~ variable)[,4] 
    PSH.Tmin <- as.data.frame(cast(CT.Temp.PSH, Sampling.Unit.Name ~ Year ~ variable)[,,1])
    PSH.Tmax <- as.data.frame(cast(CT.Temp.PSH, Sampling.Unit.Name ~ Year ~ variable)[,,2])
    PSH.Tvar <- round(as.data.frame(cast(CT.Temp.PSH, Sampling.Unit.Name ~ Year ~ variable)[,,3]),2)
    PSH.SD <- round(as.data.frame(cast(CT.Temp.PSH, Sampling.Unit.Name ~ Year ~ variable)[,,4]),2)
    PSH.Mean <- as.data.frame(cast(CT.Temp.PSH, Sampling.Unit.Name ~ Year ~ variable)[,,5])


      PSH.Elev <- scale(PSH.Elev)
      PSH.Elev[is.na(PSH.Elev)] <- 0

      PSH.FL120 <- scale(PSH.FL120)
      PSH.FL120[is.na(PSH.FL120)] <- 0

      PSH.FG120 <- scale(PSH.FG120)
      PSH.FG120[is.na(PSH.FG120)] <- 0

      PSH.Tmin <- scale(PSH.Tmin)
      PSH.Tmin[is.na(PSH.Tmin)] <- 0

      PSH.Tmax <- scale(PSH.Tmax)
      PSH.Tmax[is.na(PSH.Tmax)] <- 0

      PSH.Tvar <- scale(PSH.Tvar)
      PSH.Tvar[is.na(PSH.Tvar)] <- 0

      PSH.SD <- scale(PSH.SD)
      PSH.SD[is.na(PSH.SD)] <- 0

      PSH.Mean <- scale(PSH.Mean)
      PSH.Mean[is.na(PSH.Mean)] <- 0

# Limit PSH covariate data to array 1 (rows 1:30)
          PSH_covs <- list(Elevation=PSH.Elev[1:30],
                ForestLossCT=PSH.FL120[1:30],
                ForestGainCT=PSH.FG120[1:30],
                Tmin=PSH.Tmin[1:30,],
                Tmax=PSH.Tmax[1:30,],
                Tvar=PSH.Tvar[1:30,],
                Tsd=PSH.SD[1:30,],
                Tmean=PSH.Mean[1:30,])

# Need to remove CT-YAN-2-16.1 and CT-YAN-1-04 from Elevation data 
  CT.Temp.YAN <- melt(Alltemp500[Alltemp500$Site.Code=="YAN",])
  ELEV_FL.YAN <- melt(ELEV_FL[ELEV_FL$Site.Code=="YAN",])
    YAN.Elev <- cast(ELEV_FL.YAN, Sampling.Unit.Name ~ variable)[,2]
    YAN.FL120 <- cast(ELEV_FL.YAN, Sampling.Unit.Name ~ variable)[,3]
    YAN.FG120 <- cast(ELEV_FL.YAN, Sampling.Unit.Name ~ variable)[,4] 
    YAN.Tmin <- as.data.frame(cast(CT.Temp.YAN, Sampling.Unit.Name ~ Year ~ variable)[,,1])
    YAN.Tmax <- as.data.frame(cast(CT.Temp.YAN, Sampling.Unit.Name ~ Year ~ variable)[,,2])
    YAN.Tvar <- round(as.data.frame(cast(CT.Temp.YAN, Sampling.Unit.Name ~ Year ~ variable)[,,3]),2)
    YAN.SD <- round(as.data.frame(cast(CT.Temp.YAN, Sampling.Unit.Name ~ Year ~ variable)[,,4]),2)
    YAN.Mean <- as.data.frame(cast(CT.Temp.YAN, Sampling.Unit.Name ~ Year ~ variable)[,,5])


      YAN.Elev <- scale(YAN.Elev)
      YAN.Elev[is.na(YAN.Elev)] <- 0

      YAN.FL120 <- scale(YAN.FL120)
      YAN.FL120[is.na(YAN.FL120)] <- 0

      YAN.FG120 <- scale(YAN.FG120)
      YAN.FG120[is.na(YAN.FG120)] <- 0

      YAN.Tmin <- scale(YAN.Tmin)
      YAN.Tmin[is.na(YAN.Tmin)] <- 0

      YAN.Tmax <- scale(YAN.Tmax)
      YAN.Tmax[is.na(YAN.Tmax)] <- 0

      YAN.Tvar <- scale(YAN.Tvar)
      YAN.Tvar[is.na(YAN.Tvar)] <- 0

      YAN.SD <- scale(YAN.SD)
      YAN.SD[is.na(YAN.SD)] <- 0

      YAN.Mean <- scale(YAN.Mean)
      YAN.Mean[is.na(YAN.Mean)] <- 0


          YAN_covs <- list(Elevation=YAN.Elev,
                ForestLossCT=YAN.FL120,
                ForestGainCT=YAN.FG120,
                Tmin=YAN.Tmin,
                Tmax=YAN.Tmax,
                Tvar=YAN.Tvar,
                Tsd=YAN.SD,
                Tmean=YAN.Mean)

  CT.Temp.NAK <- melt(Alltemp500[Alltemp500$Site.Code=="NAK",])
  ELEV_FL.NAK <- melt(ELEV_FL[ELEV_FL$Site.Code=="NAK",])
    NAK.Elev <- cast(ELEV_FL.NAK, Sampling.Unit.Name ~ variable)[,2]
    NAK.FL120 <- cast(ELEV_FL.NAK, Sampling.Unit.Name ~ variable)[,3]
    NAK.FG120 <- cast(ELEV_FL.NAK, Sampling.Unit.Name ~ variable)[,4] 
    NAK.Tmin <- as.data.frame(cast(CT.Temp.NAK, Sampling.Unit.Name ~ Year ~ variable)[,,1])
    NAK.Tmax <- as.data.frame(cast(CT.Temp.NAK, Sampling.Unit.Name ~ Year ~ variable)[,,2])
    NAK.Tvar <- round(as.data.frame(cast(CT.Temp.NAK, Sampling.Unit.Name ~ Year ~ variable)[,,3]),2)
    NAK.SD <- round(as.data.frame(cast(CT.Temp.NAK, Sampling.Unit.Name ~ Year ~ variable)[,,4]),2)
    NAK.Mean <- as.data.frame(cast(CT.Temp.NAK, Sampling.Unit.Name ~ Year ~ variable)[,,5])


      NAK.Elev <- scale(NAK.Elev)
      NAK.Elev[is.na(NAK.Elev)] <- 0

      NAK.FL120 <- scale(NAK.FL120)
      NAK.FL120[is.na(NAK.FL120)] <- 0

      NAK.FG120 <- scale(NAK.FG120)
      NAK.FG120[is.na(NAK.FG120)] <- 0

      NAK.Tmin <- scale(NAK.Tmin)
      NAK.Tmin[is.na(NAK.Tmin)] <- 0

      NAK.Tmax <- scale(NAK.Tmax)
      NAK.Tmax[is.na(NAK.Tmax)] <- 0

      NAK.Tvar <- scale(NAK.Tvar)
      NAK.Tvar[is.na(NAK.Tvar)] <- 0

      NAK.SD <- scale(NAK.SD)
      NAK.SD[is.na(NAK.SD)] <- 0

      NAK.Mean <- scale(NAK.Mean)
      NAK.Mean[is.na(NAK.Mean)] <- 0


          NAK_covs <- list(Elevation=NAK.Elev,
                ForestLossCT=NAK.FL120,
                ForestGainCT=NAK.FG120,
                Tmin=NAK.Tmin,
                Tmax=NAK.Tmax,
                Tvar=NAK.Tvar,
                Tsd=NAK.SD,
                Tmean=NAK.Mean)

  CT.Temp.RNF <- melt(Alltemp500[Alltemp500$Site.Code=="RNF",])
  ELEV_FL.RNF <- melt(ELEV_FL[ELEV_FL$Site.Code=="RNF",])
    RNF.Elev <- cast(ELEV_FL.RNF, Sampling.Unit.Name ~ variable)[,2]
    RNF.FL120 <- cast(ELEV_FL.RNF, Sampling.Unit.Name ~ variable)[,3]
    RNF.FG120 <- cast(ELEV_FL.RNF, Sampling.Unit.Name ~ variable)[,4] 
    RNF.Tmin <- as.data.frame(cast(CT.Temp.RNF, Sampling.Unit.Name ~ Year ~ variable)[,,1])
    RNF.Tmax <- as.data.frame(cast(CT.Temp.RNF, Sampling.Unit.Name ~ Year ~ variable)[,,2])
    RNF.Tvar <- round(as.data.frame(cast(CT.Temp.RNF, Sampling.Unit.Name ~ Year ~ variable)[,,3]),2)
    RNF.SD <- round(as.data.frame(cast(CT.Temp.RNF, Sampling.Unit.Name ~ Year ~ variable)[,,4]),2)
    RNF.Mean <- as.data.frame(cast(CT.Temp.RNF, Sampling.Unit.Name ~ Year ~ variable)[,,5])


      RNF.Elev <- scale(RNF.Elev)
      RNF.Elev[is.na(RNF.Elev)] <- 0

      RNF.FL120 <- scale(RNF.FL120)
      RNF.FL120[is.na(RNF.FL120)] <- 0

      RNF.FG120 <- scale(RNF.FG120)
      RNF.FG120[is.na(RNF.FG120)] <- 0

      RNF.Tmin <- RNF.Tmin[,1:6] # SUBSET SO THAT # YEARS MATCHES FOR SP & TEMP DATA
      RNF.Tmin <- scale(RNF.Tmin)
      RNF.Tmin[is.na(RNF.Tmin)] <- 0

      RNF.Tmax <- RNF.Tmax[,1:6]
      RNF.Tmax <- scale(RNF.Tmax)
      RNF.Tmax[is.na(RNF.Tmax)] <- 0

      RNF.Tvar <- RNF.Tvar[,1:6]
      RNF.Tvar <- scale(RNF.Tvar)
      RNF.Tvar[is.na(RNF.Tvar)] <- 0

      RNF.SD <- RNF.SD[,1:6]
      RNF.SD <- scale(RNF.SD)
      RNF.SD[is.na(RNF.SD)] <- 0

      RNF.Mean <- RNF.Mean[,1:6]
      RNF.Mean <- scale(RNF.Mean)
      RNF.Mean[is.na(RNF.Mean)] <- 0


          RNF_covs <- list(Elevation=RNF.Elev,
                ForestLossCT=RNF.FL120,
                ForestGainCT=RNF.FG120,
                Tmin=RNF.Tmin,
                Tmax=RNF.Tmax,
                Tvar=RNF.Tvar,
                Tsd=RNF.SD,
                Tmean=RNF.Mean)

# Remove X year from NAK to match available species data and temperature data
# NEED TO DO THIS TO MATCH UP ALL SITES SO THAT MODELS WILL CONVERGE! (e.g. chimpanzees at BIF won't converge w/out correction)
# BIF_covs[[6]] <- lapply(BIF_covs[[4:8]], "[", ,1:6)
          
All_covs <- list(VB__covs=VB__covs, #9 YEARS
                 UDZ_covs=UDZ_covs, #6 YEARS
                 BIF_covs=BIF_covs, #7 YEARS (sp only 6)
                 PSH_covs=PSH_covs, #6 YEARS
                 YAN_covs=YAN_covs, #5 YEARS
                 NAK_covs=NAK_covs, #5 YEARS
                 RNF_covs=RNF_covs) #7 YEARS (SP ONLY 6)

save(All_covs, file="All_covs.RData")
#save(All_covs, file="All_covs_unscaled.RData")

################ SAVE OBJECTS CONTAINING SITE LEVEL COVARIATES
#save(VB_covs, file="VB_covs.RData")
#save(UDZ_covs, file="UDZ_covs.RData")
#save(BIF_covs, file="BIF_covs.RData")
#save(PSH_covs, file="PSH_covs.RData")
#save(YAN_covs, file="YAN_covs.RData")
#save(NAK_covs, file="NAK_covs.RData")
#save(RNF_covs, file="RNF_covs.RData")


# Format EDI from Miguel for non-temporally varying covariate
# Need to create a vector for each species for the site CTs based on CT communities in Z.
load("Scaled_FPDist_0.5a.RData")
SitesBinary[[1]][,1] # Species index values for the site. To generalize, change to SitesBinary[[i]][,1]
SitesBinary[[1]][1,1] # Single species value to then extract elements from Z. To generalize, change to SitesBinary[[i]][j,1]
Z[[1]][1]

# the following loop produces output for a single species (#21) for all camera traps
test <- vector()
for(i in 1:length(Z[[1]])){
    test[i] <- Z[[1]][[i]]["21"]
  test[i] <- ifelse(is.na(test[i])==TRUE, 0, test[i])
}

# Now extend to all species for one TEAM site
test <- vector()
hold <-list()
#hold <- data.frame(Sampling.Unit.Name=colnames(SitesBinary[[1]])[2:length(colnames(SitesBinary[[1]]))])

  for(j in 1:length(SitesBinary[[1]][,1])){
    for(i in 1:length(Z[[1]])){
    test[i] <- Z[[1]][[i]][as.character(SitesBinary[[1]][,1])[j]]
    test[i] <- ifelse(is.na(test[i])==TRUE, 0, test[i]) 
    hold[[j]] <- test
  } 
}

names(hold) <- rownames(SitesBinary[[1]]) # Brings back species names, but could also name with UID index if needed later

# Now extend to all TEAM sites
load("Scaled_FPDist_0.5a.RData")
BIOTIC_pop <- vector()
BIOTIC_site <-list()
BIOTIC_all <- list()

for(k in 1:length(SitesBinary)){
  for(j in 1:length(SitesBinary[[k]][,1])){
    for(i in 1:length(Z[[k]])){
    BIOTIC_pop[i] <- Z[[k]][[i]][as.character(SitesBinary[[k]][,1])[j]]
    BIOTIC_pop[i] <- ifelse(is.na(BIOTIC_pop[i])==TRUE, 0, BIOTIC_pop[i]) 
    }   
    BIOTIC_site[[j]] <- BIOTIC_pop
    rm(BIOTIC_pop) # removes object and replaces is that the correct number of CT are used 
    BIOTIC_pop <- vector()
  }
  BIOTIC_all[[k]] <- BIOTIC_site
    names(BIOTIC_all[[k]]) <- rownames(SitesBinary[[k]])  # Brings back species names, but could also name with UID index if needed later
   rm(BIOTIC_site)
  BIOTIC_site <- list()
}
names(BIOTIC_all) <- Sitenames

# Create object that is a single list (rather than 7 lists of lists) to use as input for unmarked for easy indexing
BIOTIC_166 <- c(BIOTIC_all[[1]], BIOTIC_all[[2]], BIOTIC_all[[3]], BIOTIC_all[[4]], BIOTIC_all[[5]], BIOTIC_all[[6]], BIOTIC_all[[7]])
save(BIOTIC_166, file="BIOTIC_all.RData")

# Check to be sure that population level data is in the same order for y input and for BIOTIC covariate input
cbind(names(All_species7sites), names(BIOTIC_166))


# Format EDI from Miguel for temporally varying covariate
load("Scaled_FPDist_time_0.5a.RData")

BIOTIC_pop <- vector()
BIOTIC_year <- list()
BIOTIC_site <-list()
BIOTIC_all <- list()

for(k in 1:length(SitesBinaryAnnual)){
  for(m in 1:length(SitesBinaryAnnual[[k]])){ 
    for(j in 1:length(SitesBinaryAnnual[[k]][[m]][,1])){
      for(i in 1:length(Z[[k]][[m]])){
        BIOTIC_pop[i] <- Z[[k]][[m]][[i]][as.character(SitesBinaryAnnual[[k]][[m]][,1])[j]]
        BIOTIC_pop[i] <- ifelse(is.na(BIOTIC_pop[i])==TRUE, 0, BIOTIC_pop[i]) 
        }   
      BIOTIC_site[[j]] <- BIOTIC_pop
      rm(BIOTIC_pop) # removes object and replaces is that the correct number of CT are used 
      BIOTIC_pop <- list()
      }
    BIOTIC_year[[m]] <- BIOTIC_site
    rm(BIOTIC_site)
    BIOTIC_site <- list()
    }
  BIOTIC_all[[k]] <- BIOTIC_year
  #names(BIOTIC_all[[k]][[m]]) <- rownames(SitesBinaryAnnual[[k]][[m]])  # Brings back species names, but could also name with UID index if needed later
  rm(BIOTIC_year)
  BIOTIC_year <- list()
  }

names(BIOTIC_all) <- Sitenames

# Now convert BIOTIC_all output to a list of 166 dataframes (one for each population) with rows for each CT and columns for each year of data
VBannual <- melt(BIOTIC_all[[1]]) #L1 is nyears; L2 is species; L3 is CT
VBannual <- VBannual[61:dim(VBannual)[1],]
VBannual_cast <- cast(VBannual, L3 ~ L1 ~ L2)
VBannual_array <- alply(VBannual_cast,3)
VBannual_BIOTIC <- llply(VBannual_array, data.frame)

SiteList <- list()
for(i in 2:length(SitesBinaryAnnual)){
  test <- melt(BIOTIC_all[[i]]) #L1 is nyears; L2 is species; L3 is CT
  test2 <- cast(test, L3 ~ L1 ~ L2)
  test3 <- alply(test2,3)
  test4 <- llply(test3, data.frame)
  SiteList[[i]] <- test4
}

# Remove 5th year from NAK to match available CT temperature data
SiteList[[6]] <- lapply(SiteList[[6]], "[", ,1:4)

BIOTIC_ALL_YEARS <- c(VBannual_BIOTIC, SiteList[[2]], SiteList[[3]], SiteList[[4]], SiteList[[5]], SiteList[[6]], SiteList[[7]])

save(BIOTIC_ALL_YEARS, file="BIOTIC_ALL_YEARS.RData")






# Subset input data to exclude binomial cases (i.e rare populations) based on Cases$Full
# Note that this reduces the # of populations from 166 to 62 (35 constant; 27 simple)
Cases <- read.csv("Cases.csv")
Full <- paste(Cases$site, "species", Cases$site.sp, sep=".")
Cases <- cbind(Cases, Full)
Cases <- Cases[match(names(All_species7sites), Cases$Full),]
# NA value species are Helarctos malaynus at NAK and Tapirus terrestris at VB

# Subset the following objects based on Cases: All_species7sites, BIOTIC_166, BIOTIC_ALL_YEARS

ExcludeBinomial <- Cases[Cases$Case!="binomial",]
ExcludeBinomial <- na.omit(ExcludeBinomial)
IncludeIndex <- match(ExcludeBinomial$Full, names(All_species7sites))

Species7sites_Include <- All_species7sites[IncludeIndex]
BIOTIC_Include <- BIOTIC_166[IncludeIndex]
BIOTIC_ALL_YEARS_Include <- BIOTIC_ALL_YEARS[IncludeIndex]

save(Species7sites_Include, file="Species7sites_Include.RData")
save(BIOTIC_Include, file="BIOTIC_Include.RData")
save(BIOTIC_ALL_YEARS_Include, file="BIOTIC_ALL_YEARS_Include.RData")




######################## FUNCTION TO FORMAT TRINARY (1/0/NA) CT DATA MATRICES ###############
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
