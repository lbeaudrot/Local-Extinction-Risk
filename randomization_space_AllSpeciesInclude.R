############################
#Matrix randomization
############################
#By: Miguelito, Nov 2, 2016, 2:24pm 
#Modified for all 62 species in nms by: Lydia Beaudrot, Dec 16, 2016, 3:29 pm 

library(unmarked)

load('BIF.Cephalophus_nigrifrons.RData')
#bcn=getY(umf)
bcn <- umf@y
rm(umf)

load("Species7sites_Include.RData")
load("All_covs.RData")


Sp_Data <- Species7sites_Include

# Create loop that goes through Species_Data and creates a list with the umf for each population

umf_list <- list()

for(m in 1:length(Sp_Data)){
  print(m)
  
  # Define species for analysis and site using index value for list of all species
  sp.name <- names(Species_data)[m]
  species <- Species_data[[m]]
  site <- substr(names(Species_data)[[m]],1,3)
  
  # Define site covariates
  site_covs <- paste(site, "covs", sep="_")
  covs <- All_covs[names(All_covs)==site_covs]
  Elevation <- unlist(as.matrix(sapply(covs, "[", 1)))
  site.covs<-data.frame(Elevation)
  
  # Define number of primary periods
  yrs <- as.data.frame(sapply(covs, "[", 4))
  to=dim(yrs)[2]
  
  # Create object with data formatted for unmarked
  umf<-unmarkedMultFrame(y=species, siteCovs=site.covs, numPrimary=dim(yrs)[2])
  umf_list[[m]] <- umf@y
  
}





#datos=list(bcn,bpt,ntk,ucg,upl,vdn,vpt,yeb,yma)
datos <- umf_list
#reshuffled.col=vector("list",9) #reshuffling columns (i.e. time)
#reshuffled.all=vector("list",9) #reshuffling all values in the matrix
reshuffled.row=vector("list", length(umf_list)) #reshuffling rows(space)

for (i in 1:length(umf_list)){
  print(i)
  nr=dim(datos[[i]])[1]
  nc=dim(datos[[i]])[2]	
  for (j in 1:250){ #100 random matrices (it could be more or less depending on comp time)
    #reshuffled.col[[i]][[j]]=datos[[i]][,sample(1:nc)]
    #reshuffled.all[[i]][[j]]=matrix(sample(datos[[i]]),nrow=nr)
    reshuffled.row[[i]][[j]]=datos[[i]][sample(1:nr),]
  }
}

names(reshuffled.row) <- names(Sp_Data)

#save(reshuffled.col, file="reshuffled.col.RData")
#save(reshuffled.all, file="reshuffled.all.RData")
save(reshuffled.row, file="reshuffled.row.RData")
