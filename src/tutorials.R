# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 27.10.23
# Last modified: 31.11.23

path <- "/Users/christian/Library/CloudStorage/GoogleDrive-christian.tsoungui@aims-cameroon.org/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Christian/Models/FullGeneralModel"

# Load external resources
source(paste0(path,'/src/model.R'))#("/home/janedoe/Documents/src/STRmodel.R")

# Install the necessary packages if necessary
#install.packages('openxlsx')   # Comment this line if openxlsx installed

# Loading libraries
library(openxlsx)

#################################
### Import Datasets
##################################

## Import the dataset
datasetNaturalFormat <- read.xlsx(paste0(path,'/exampleDataset/exampleDatasetNaturalFormat.xlsx'), 1)

# Transform the data to the standard format
datasetStandard <- convertDatasetToStandardFormat(datasetNaturalFormat, 2:ncol(datasetNaturalFormat))

#################################
### Estimate MLEs
##################################

# Estimate MLEs
## Choose markers of interests
markers <- 1:4
calculateMaximumLikelihoodEstimatesWithAddOns(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE)

# Estimate MLEs with plugin estimate for MOI parameter
calculateMaximumLikelihoodEstimatesWithAddOns(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, pluginValueOfLambda = 1.0)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction ('Bootstrap')
calculateMaximumLikelihoodEstimatesWithAddOns(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, isBiasCorrection = TRUE, methodForBiasCorrection = "bootstrap", numberOfBootstrapReplicatesBiasCorrection = 15000)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction ('jackknife') with plugin
calculateMaximumLikelihoodEstimatesWithAddOns(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, pluginValueOfLambda = 1.0, isBiasCorrection = TRUE, methodForBiasCorrection = "jackknife")

# Finding MLEs (haplotype frequencies and MOI) using a 95% confidence interval
calculateMaximumLikelihoodEstimatesWithAddOns(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, isBiasCorrection = TRUE, methodForBiasCorrection = "bootstrap", numberOfBootstrapReplicatesBiasCorrection = 15000, isConfidenceInterval = TRUE, numberOfBootstrapReplicatesConfidenceInterval = 10000)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 20000 bootstrap samples
calculateMaximumLikelihoodEstimatesWithAddOns(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, isBiasCorrection = TRUE, methodForBiasCorrection = "bootstrap", numberOfBootstrapReplicatesBiasCorrection = 15000,  isConfidenceInterval = TRUE, numberOfBootstrapReplicatesConfidenceInterval = 20000, significanceLevel = 0.1)

# Finding pairwise LD between two loci. The function outputs the LD measures D', r-squared.
markersPair <- c(4,4)
calculatePairwiseLDWithAddons(datasetStandard,markersPair, idExists = FALSE)

# Finding pairwise LD between two loci. The function outputs the LD measures D', r-squared using a 95% confidence interval
calculatePairwiseLDWithAddons(datasetStandard,markersPair, idExists = FALSE, isConfidenceInterval=TRUE,numberOfBootstrapReplicatesConfidenceInterval=100, significanceLevel=0.05)
