# Convert interpolated temperature data into useable form for unmarked, etc.
# We need a matrix with each unique camera trap as a row and annual temperature variables as columns

library(reshape)

# Read in output file from Interpolation script that combines min, max and var into a single cv file 
temp <- read.csv("VB_InterpolatedTemperatures.csv") # Interpolated temperature data
Orig.temp <- read.csv("VB_Covariate_Data.csv") # Original (non-interpolated) temperature data

# Examine # of camera traps in each file
length(unique(temp$Sampling.Unit.Name))
length(unique(Orig.temp$Sampling.Unit))
head(temp)
dim(temp)

# Use cast to reshape data (data may already be melted....so may not need to do that, but explore)