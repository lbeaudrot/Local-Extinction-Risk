# Convert interpolated temperature data into useable form for unmarked, etc.
# We need a matrix with each unique camera trap as a row and annual temperature variables as columns

library(reshape)
library(ggplot2)
# Read in forest loss data from Alex
traps_fc <- read.csv("traps_fc.csv")

# Read in output file from Interpolation script that combines min, max and var into a single cv file 
# VB
#temp <- read.csv("VB_InterpolatedTemperatures.csv") # Interpolated temperature data
#Orig.temp <- read.csv("VB_Covariate_Data.csv") # Original (non-interpolated) temperature data

temp <- read.csv("UDZ_InterpolatedTemperatures.csv") # Interpolated temperature data

# First year of temperature data are missing in interpolated data file
# Check with DSSG team about input data and/or correction and for the meantime use originals for these 

# Examine # of camera traps in each file
length(unique(temp$Sampling.Unit.Name))
length(unique(Orig.temp$Sampling.Unit))
head(temp)
dim(temp)

# Use cast to reshape data (data may already be melted....so may not need to do that, but explore)

temp.melt <- melt(temp, id.vars=c("Sampling.Unit.Name", "Time"))
temp.cast <- cast(temp.melt, Sampling.Unit.Name ~ variable + Time)
#write.csv(temp.cast, file="VB.temp.cast.csv")
#write.csv(temp.cast, file="UDZ.temp.cast.csv")



# Extract forest loss data from Alex's object "traps_fc"
VBtraps <- traps_fc[traps_fc$sitecode=="VB",]
UDZtraps <- traps_fc[traps_fc$sitecode=="UDZ" & traps_fc$buffer_m=="120m buffer",]
#write.csv(UDZtraps, file="UDZtraps.csv")

ELEV <- read.csv("CT_edgedist_elevation_final.txt")
#ELEVsub <- ELEV[match(temp.cast$Sampling.Unit.Name, ELEV$Sampling.Unit.Name),]

ftraps120 <- traps_fc[traps_fc$buffer_m=="120m buffer",]
ftraps120 <- cbind(ftraps120, ELEV[match(ftraps120$trap_ID, ELEV$Sampling.Unit.Name),])
ftraps60 <- traps_fc[traps_fc$buffer_m=="60m buffer",]
ftraps60 <- cbind(ftraps60, ELEV[match(ftraps60$trap_ID, ELEV$Sampling.Unit.Name),])
ftraps30 <- traps_fc[traps_fc$buffer_m=="30m buffer",]
ftraps30 <- cbind(ftraps30, ELEV[match(ftraps60$trap_ID, ELEV$Sampling.Unit.Name),])

ftrapsall <- rbind(ftraps120, ftraps60, ftraps30)

# Examine boxplots of camera traps by elevation
ggplot(ftrapsall, aes(sitecode, Elevation)) +
    geom_boxplot() +
    facet_grid(buffer_m~.) +
    xlab("Site") +
    ylab("Elevation")

# Examine the distribution of camera trap points by elevation
ggplot(ftraps120, aes(sitecode, Elevation)) +
    geom_point() +
   
    xlab("Site") +
    ylab("Elevation")

ftraps120.subset <- ftraps120[ftraps120$sitecode=="UDZ"|ftraps120$sitecode=="BIF"|ftraps120$sitecode=="YAN"|ftraps120$sitecode=="PSH"|ftraps120$sitecode=="NAK",]
# Examine the distribution of camera trap points by elevation
ggplot(ftraps120, aes(Elevation, fc_frac_loss)) + 
    ylim(0, 0.2) +
    geom_point() +
    facet_grid(sitecode~.) +
    xlab("Elevation") +
    ylab("Forest Loss")
    




