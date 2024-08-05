#install.packages("drc")
library(drc)
#install.packages("xlsx")
library(readxl)
#install.package("dplyr")
library(dplyr)

# Dataset needs to have sample labels in first column, dilutions in second column, and ODs in third column. The first row should be the column labels.

### ---------------------
#### Standard curve
### ---------------------

# Import standard data
dat.standards <- read_excel("2023.07.21 HA ELISA STAPH#17 Skin/test.standards.xlsx")

# Calculate the average blank reading
blanks <- dat.standards[dat.standards$Label=="Blank",]
avg.blank <- mean(blanks$OD)

# Correct standard values by subtracting out average blank reading
corrected.standards <- dat.standards$OD-avg.blank

# Make a new dataset to replace the old readings with the corrected readings
dat.standards.corrected <- dat.standards
dat.standards.corrected$OD <- corrected.standards

# Remove blanks from dataset
dat.standards.corrected <- dat.standards.corrected[-c(15,16),]

# Fit 4PL model
st.curve <- drm(OD~Conc, fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")), data=dat.standards.corrected)

# TEMPORARY CORRECTION for slope having the wrong sign
st.curve$coefficients[1] <- abs(st.curve$coefficients[1])

# See 4PL parameter estimates
summary(st.curve)

# Plot 4PL standard curve
plot(st.curve)


### ---------------------
### Samples
### ---------------------

# Import sample data
dat.samples <- read_excel("2023.07.21 HA ELISA STAPH#17 Skin/test.samples.xlsx", col_types = c("text","numeric","numeric"))

# Convert dilution factor to concentration (leave dilution factor for later)
dat.samples$Conc <- 1/dat.samples$Dilution

# Correct standard values by subtracting out average blank reading
corrected.samples <- dat.samples$OD-avg.blank

# Make a new dataset to replace the old readings with the corrected readings
dat.samples.corrected <- dat.samples
dat.samples.corrected$OD <- round(corrected.samples, digits=3)

# Estimate the concentrations
DOSEx <- ED(st.curve, dat.samples.corrected$OD, type="absolute", display=F) #NaN values are lower than the curve

# Turn DOSEx into a dataframe
DOSEx <- data.frame(DOSEx)

# Change labels the proper labels
rownames(DOSEx) <- dat.samples$Label

# Scale estimates by dilution factors
DOSEx$Estimate <- DOSEx$Estimate*dat.samples$Dilution   # do not run this line multiple times; it will multiply again. If you do, restart from the first DOSEx line

# Add ODs back in for choosing which value to use
DOSEx$OD <- dat.samples$OD
DOSEx <- relocate(DOSEx, OD, Estimate)

# Remove values not on the curve
max.standard <- dat.standards[which(dat.standards$Conc==40),] # Isolate 40 ng standards
avg.max.standard <- mean(max.standard$OD)   # Calculate average OD of 40 ng standards

results <- subset(DOSEx, OD < avg.max.standard) # Filter the dataset to remove samples with larger OD than the average of the 40 ng standards

results <- na.omit(results) # Remove NaN values


# Print results!
round(results, digits=2)

### ---------------------
#### Sources
### ---------------------

# Using drc package
# https://katatrepsis.com/2018/09/24/elisa-analysis-in-r/

# Uinsg nlsLM package/manual method
# https://janalin.github.io/analyse-ELISA/calibration.html
