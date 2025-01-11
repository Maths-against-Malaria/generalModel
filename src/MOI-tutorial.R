# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 27.10.23
# Last modified: 11.01.25

# Load external resources
source('/Users/christian/Documents/Phd/models/generalModel/src/MOI-MLE.R')

# Load libraries
library(openxlsx)

#################################
### Import Datasets
##################################
data <- read.xlsx('/Users/christian/Documents/Phd/models/generalModel/exampleDataset/dataset.xlsx', 1)

#################################
### Estimate MLEs
##################################
## Choose markers of interests
markers <- 1:4
MLE(data, markers)

# Estimate MLEs with plugin estimate for MOI parameter
MLE(data, markers, plugin = 1.0, allelesName = FALSE)

# Finding MLEs (haplotype frequencies and MOI) with bootstrap bias-correction
MLE(data, markers, isBC = TRUE, replBC = 15000)

# Finding MLEs (haplotype frequencies and MOI) using a 90% confidence interval and 15000 bootstrap replicates
MLE(data, markers, isCI = TRUE, replCI = 15000, alpha = 0.10)

# Finding MLEs (haplotype frequencies and MOI) with bias-correction and a 90% confidence
MLE(data, markers, isBC = TRUE, replBC = 20000, isCI = TRUE, replCI = 15000, alpha = 0.10)

# Finding pairwise LD between two loci (,i.e., at 1st and 4th column), using D' and r-squared with a 90% confidence interval
markersPair <- c(1,4)
pairwiseLD(data, markersPair, isCI = TRUE, replCI = 20000, alpha = 0.10)

#################################
### Estimate prevalence
##################################
## Estimating haplotype prevalence and MOI using a 90% confidence interval and 15000 bootstrap replicates
markers <- 1:2
PREV(data, markers, isCI = TRUE, replCI = 15000, alpha = 0.10)

#################################
### Asymptotic variance of MLEs
##################################
## Calculate covariance matrix for MOI parameter and frequencies
markers <- 1:2
FI(data, markers)

## Calculate covariance matrix for mean MOI and frequencies
FI(data, markers, isPsi = TRUE)

## Calculate covariance matrix for mean MOI and prevalence
FI(data, markers, isPsi = TRUE, isPrev = TRUE, allelesName = FALSE, isObserv = FALSE)
