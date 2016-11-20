# Follow up analysis based on Nature Reviewer #1 and Reviewer #4 comments
# Compare annual occupancy values to equilibrium occupancy
# Equilibrium occupancy = C / (C+E) where C is colonization and E is extinction 
# Lambda = (occupancy t + 1) / (occupancy t) directly estimates rate of change in site occupancy

# TOP MODELS FOR 9 ELEVATION SHIFT SPECIES IN NATURE SUBMISSION
# Double check model numbers for changes 
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


model <- fm2.1

# Calculate predicted annual occupancy values
timeseries <- as.numeric(smoothed(model)[2,])

elev.time <- smoothed(model, FALSE)
elev.time <- t(elev.time[2,,])

# Sort data so that it's in elevational order for display

elev_occ <- data.frame(BIF.Elev, elev.time)
elev_occ <- elev_occ[rev(do.call(order, elev_occ)),]
elev_occ <- as.matrix(elev_occ[,2:6])

# Plot changes in occupancy overtime per camera trap. Use colors to indicate elevation gradient. 
tsRainbow <- rainbow(ncol(as.zoo(t(elev.time))), start=0.5, end=1)
# Blue is highest elevations; Red is lowest elevations
# Plot with points
plot(elev_occ[1,] ~ c(1:5), las=1, xlab="Year", ylab="Occupancy Estimate")
for(i in 2:dim(elev_occ)[1]){
  points(elev_occ[i,] ~ c(1:5), col=tsRainbow[i])
}

# Plot with lines
library(zoo)
plot(as.zoo(t(elev.time)), col=tsRainbow, screens=1, las=1, ylab="Occupancy estimate", xlab="Year")

# Calculate the rate of change in site occupancy Lambda = (occupancy t + 1) / (occupancy t)
years <- length(timeseries)
lambda <- timeseries[2:years]/timeseries[1:(years-1)]

# Calculate estimates for colonization and extinction to estimate Equilibrium occupancy C/(C+E)
Cest <- coef(model, type = "col")[2]
Eest <- coef(model, type = "ext")

Peq <- Cest / (Cest + Eest)
