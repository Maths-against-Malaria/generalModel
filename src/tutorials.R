# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 27.10.23
# Last modified: 05.12.23

path <- "/Users/christian/Documents/phd/models/generalModel"

# Load external resources
source(paste0(path,'/src/model.R'))

# Install the necessary packages if necessary
#install.packages('openxlsx')   # Comment this line if openxlsx installed

# Loading libraries
library(openxlsx)

#################################
### Import Datasets
##################################
datasetNatural <- read.xlsx(paste0(path,'/exampleDatasets/dataset.xlsx'), 1)

# Transform the data to the standard format
datasetStandard <- datasetToStandard(datasetNatural, 2:ncol(datasetNatural))

#################################
### Estimate MLEs
##################################
## Choose markers of interests
markers <- 1:4
MLE(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE)

# Estimate MLEs with plugin estimate for MOI parameter
MLE(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, plugin = 1.0)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction
MLE(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, isBC = TRUE, replBC = 15000)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 15000 bootstrap replicates
MLE(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, isCI = TRUE, replCI = 15000, alpha = 0.10)

# Finding pairwise LD between two loci. The function outputs the LD measures D', r-squared.

pairwiseLD(datasetStandard, markersPair, idExists = FALSE)

# Finding pairwise LD between two loci (,i.e., at 1st and 4th column), using D' and r-squared with a 95% confidence interval
markersPair <- c(1,4)
pairwiseLD(datasetStandard, markersPair, idExists = FALSE, isCI=TRUE, replCI = 20000, alpha=0.10)
