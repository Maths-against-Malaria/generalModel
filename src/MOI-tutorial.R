# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 27.10.23
# Last modified: 06.01.25

# Load external resources
source('/home/johndoe/documents//src/MOI-MLE.R')

# Load libraries
library(openxlsx)

#################################
### Import Datasets
##################################
datasetNatural <- read.xlsx('/home/johndoe/documents/exampleDataset/dataset.xlsx', 1)

# Transform the data to the standard format
data <- datasetFormat(datasetNatural, 2:ncol(datasetNatural))

#################################
### Estimate MLEs
##################################
## Choose markers of interests
markers <- 1:4
MLE(dataset, markers, idExists = FALSE)

# Estimate MLEs with plugin estimate for MOI parameter
MLE(dataset, markers, idExists = FALSE, plugin = 1.0, allelesName = FALSE)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction
MLE(dataset, markers, idExists = FALSE, isBC = TRUE, replBC = 15000)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 15000 bootstrap replicates
MLE(dataset, markers, idExists = FALSE, isCI = TRUE, replCI = 15000, alpha = 0.10)

# Finding MLEs (haplotype frequencies and MOI) with bias-correction and a 90% confidence
MLE(dataset, markers, idExists = FALSE, isBC = TRUE, replBC = 20000, isCI = TRUE, replCI = 15000, alpha = 0.10)

# Finding pairwise LD between two loci (,i.e., at 1st and 4th column), using D' and r-squared with a 90% confidence interval
markersPair <- c(1,4)
pairwiseLD(dataset, markersPair, idExists = FALSE, isCI=TRUE, replCI = 20000, alpha=0.10)

#################################
### Estimate prevalence
##################################
## Estimating haplotype prevalence and MOI using a 90% confidence interval and 15000 bootstrap replicates
markers <- 1:2
PREV(dataset, markers, idExists = FALSE, isCI=TRUE, replCI = 15000, alpha = 0.10)

#################################
### Asymptotic variance of MLEs
##################################
## Calculate covariance matrix for MOI parameter and frequencies
markers <- 1:2
FI(dataset, markers, idExists = FALSE)

## Calculate covariance matrix for mean MOI and frequencies
FI(dataset, markers, idExists = FALSE, isPsi = TRUE)

## Calculate covariance matrix for mean MOI and prevalence
FI(dataset, markers, idExists = FALSE, isPsi = TRUE, isPrev = TRUE, allelesName = FALSE)
