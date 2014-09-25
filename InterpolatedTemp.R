# Convert interpolated temperature data into useable form for unmarked, etc.
# We need a matrix with each unique camera trap as a row and annual temperature variables as columns

library(reshape)

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
write.csv(temp.cast, file="UDZ.temp.cast.csv")



# Extract forest loss data from Alex's object "traps_fc"
VBtraps <- traps_fc[traps_fc$sitecode=="VB",]
UDZtraps <- traps_fc[traps_fc$sitecode=="UDZ" & traps_fc$buffer_m=="120m buffer",]
write.csv(UDZtraps, file="UDZtraps.csv")
