# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 27.10.23
# Last modified: 23.12.24

# Load external resources
source('/home/johndoe/documents/src/model.R')

# Install the necessary packages if necessary
#install.packages('openxlsx')   # uncomment this line to install openxlsx

# Load libraries
library(openxlsx)

#################################
### Import Datasets
##################################
datasetNatural <- read.xlsx('/home/johndoe/documents/exampleDatasets/dataset.xlsx', 1)

# Transform the data to the standard format
dataset <- datasetToStandard(datasetNatural, 2:ncol(datasetNatural))

#################################
### Estimate MLEs
##################################
## Choose markers of interests
markers <- 1:4
MLE(dataset[[1]][,markers], dataset[[3]][markers], idExists = FALSE)

# Estimate MLEs with plugin estimate for MOI parameter
MLE(dataset[[1]][,markers], dataset[[3]][markers], idExists = FALSE, plugin = 1.0)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction
MLE(dataset[[1]][,markers], dataset[[3]][markers], idExists = FALSE, isBC = TRUE, replBC = 15000)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 15000 bootstrap replicates
MLE(dataset[[1]][,markers], dataset[[3]][markers], idExists = FALSE, isCI = TRUE, replCI = 15000, alpha = 0.10)

# Finding MLEs (haplotype frequencies and MOI) with bias-correction and a 90% confidence
MLE(datasetStandard[[1]][,markers], datasetStandard[[3]][markers], idExists = FALSE, isBC = TRUE, replBC = 20000, isCI = TRUE, replCI = 15000, alpha = 0.10)

# Finding pairwise LD between two loci (,i.e., at 1st and 4th column), using D' and r-squared with a 90% confidence interval
markersPair <- c(1,4)
pairwiseLD(dataset, markersPair, idExists = FALSE, isCI=TRUE, replCI = 20000, alpha=0.10)

#################################
### Estimate prevalence
##################################
## Estimating haplotype prevalence and MOI using a 90% confidence interval and 15000 bootstrap replicates
markers <- 1:2
PREV(dataset[[1]][,markers], dataset[[3]][markers], idExists = FALSE, isCI=TRUE, replCI = 15000, alpha = 0.10)

#################################
### Asymptotic variance of MLEs
##################################
## Estimate MLEs
mle <- MLE(dataset[[1]][,markers], dataset[[3]][markers], idExists = FALSE)

## Calculate covariance matrix for MOI parameter and frequencies
CRLB(mle, dataset[[3]][markers])

## Calculate covariance matrix for mean MOI and frequencies
CRLB(mle, dataset[[3]][markers], isPsi = TRUE)

## Calculate covariance matrix for mean MOI and prevalence
CRLB(mle, dataset[[3]][markers], isPsi = TRUE, isPrev = TRUE)
