path <- "/Users/christian/Library/CloudStorage/GoogleDrive-christian.tsoungui@aims-cameroon.org/.shortcut-targets-by-id/1Ulru-DjbFRaMVB7Vj9tJ4NfyzPDkhzOr/Maths against Malaria/Christian/Models/FullGeneralModel"

# Load external resources
source(paste0(path,'/src/model.R'))

calculateScoreMatrix <- function(maximumLikelihoodEstimates, numberOfAllelesAtEachMarker, sampleSize){

  lambda <<- maximumLikelihoodEstimates[[1]]
  haplotypeFrequencies <<- maximumLikelihoodEstimates[[2]]
  detectedHaplotypes <- maximumLikelihoodEstimates[[3]]
  scoreMatrix <- matrix(0 ,nrow = length(haplotypeFrequencies), ncol = length(haplotypeFrequencies))

  numberOfLoci <- length(numberOfAllelesAtEachMarker)
  mixRadixBaseCumulativeProduct <- c(1, cumprod(numberOfAllelesAtEachMarker)[1:(numberOfLoci-1)])
  detectedHaplotypesLabels <<- (detectedHaplotypes-1)%*%mixRadixBaseCumulativeProduct+1
  rownames(haplotypeFrequencies) <- detectedHaplotypesLabels
  haplotypePairs <- t(combn(detectedHaplotypesLabels[-which(detectedHaplotypesLabels==1)],2))

  possibleObservations <- buildAllPossibleObservations(numberOfAllelesAtEachMarker)
  observationsCompatibleWithFirstHaplotype <- findObservationsCompatibleWithFirstHaplotype(numberOfAllelesAtEachMarker,possibleObservations)
  subsetsFromObservations <<- buildAllSetsAndSubsets(possibleObservations, numberOfAllelesAtEachMarker)
  listOfSplittedSubobservationsSubsets <- splitSubsetsByCompabilityWithFirstHaplotype(numberOfAllelesAtEachMarker, subsetsFromObservations)

  subsetsFromObservationsCompatibleWithFirstHaplotype <- listOfSplittedSubobservationsSubsets[[1]]
  subsetsFromObservationsUncompatibleWithFirstHaplotype <- listOfSplittedSubobservationsSubsets[[2]]

  # Spsi,psi
  temp1 <- 0
  temp2 <- 0
  temp3 <- 0

  for(compatibleObservation in which(observationsCompatibleWithFirstHaplotype)){
    probabilityDensityFunctionPx <- calculatePDF(compatibleObservation, subsetsFromObservations[[1]]) # Px
    temporaryFirstDerivativePx <- calculateFirstDerivativeByLambdaPDF(compatibleObservation, subsetsFromObservations[[1]]) #dPx/dlambda

    subsetCompatibleObservation <- which(observationsCompatibleWithFirstHaplotype)[-which(which(observationsCompatibleWithFirstHaplotype)==compatibleObservation)]
    temporaryFirstDerivativePy <- 0
    for(subObservations in subsetCompatibleObservation){
      temporaryFirstDerivativePy <- temporaryFirstDerivativePy + calculateFirstDerivativeByLambdaPDF(subObservations, subsetsFromObservations[[1]]) # dPy/dlambda
    }

    temp1 <- temp1 + temporaryFirstDerivativePx*temporaryFirstDerivativePy
    temp2 <- temp2 + (sampleSize^2 - sampleSize + sampleSize/probabilityDensityFunctionPx)*temporaryFirstDerivativePx^2

    tempSubObservation <- 0
    for(subObservation in which(!observationsCompatibleWithFirstHaplotype)){
      probabilityDensityFunctionPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
      PDFUncompatibleWithFirstHaplotype <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py(2)
      temporaryFirstDerivativePy <- calculateFirstDerivativeByLambdaPDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # dPy(2)/dlamda
      tempSubObservation <- tempSubObservation + temporaryFirstDerivativePy*probabilityDensityFunctionPy/PDFUncompatibleWithFirstHaplotype
    }

    temp3 <- temp3 + temporaryFirstDerivativePx*tempSubObservation
  }
  temp4 <- 0
  temp5 <- 0
  for(uncompatibleObservation in which(!observationsCompatibleWithFirstHaplotype)){
    temporaryPx <- calculatePDF(uncompatibleObservation, subsetsFromObservations[[1]]) # Px
    PDFUncompatibleWithFirstHaplotype <- calculatePDF(uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Px2
    temporaryFirstDerivativePx <- calculateFirstDerivativeByLambdaPDF(uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dlambda
    temp <- temporaryFirstDerivativePx/PDFUncompatibleWithFirstHaplotype

    subsetUncompatibleObservation <- which(!observationsCompatibleWithFirstHaplotype)[-which(which(!observationsCompatibleWithFirstHaplotype)==uncompatibleObservation)]
    temporarySubObservation <- 0
    for(subObservation in subsetUncompatibleObservation){
      temporaryPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
      PDFUncompatibleWithFirstHaplotypePy <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py(2)
      temporaryFirstDerivativePy <- calculateFirstDerivativeByLambdaPDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # dPy(2)/dlambda
      temporarySubObservation <- temporarySubObservation + temporaryFirstDerivativePy*temporaryPy/PDFUncompatibleWithFirstHaplotypePy
    }
    temp4 <- temp4 + temp*temporarySubObservation*temporaryPx
    temp5 <- temp5 + ((temporaryPx*sampleSize)^2 - sampleSize*temporaryPx^2 + sampleSize*temporaryPx)*temp^2
  }

  scoreMatrix[1,1] <- (sampleSize*(sampleSize - 1)*(temp1 + 2*temp3 + temp4) + temp2 + temp5)/calculateFirstDerivativeMeanMOIFunction(lambda)^2

  # Spsi,pi
  for(haplotype in detectedHaplotypesLabels[-which(detectedHaplotypesLabels==1)]){

    temp1 <- 0
    temp2 <- 0
    temp3 <- 0
    temp4 <- 0
    temp5 <- 0
    temp6 <- 0

    for(compatibleObservation in which(observationsCompatibleWithFirstHaplotype)){
      #PDFCompatiblePx <- calculatePDFPartial(compatibleObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) # Px1
      #PDFUncompatiblePx <- calculatePDF(compatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Px2
      PDFPx <- calculatePDF(compatibleObservation, subsetsFromObservations[[1]]) # Px

      temporaryFirstDerivativeCompatiblePxi <- calculateFirstDerivativeByFrequencyPDFPartial(haplotype, compatibleObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) #dPx1/dPi
      temporaryFirstDerivativeUncompatiblePxi <- calculateFirstDerivativeByFrequencyPDF(haplotype, compatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPi
      temporaryFirstDerivativePxi <- temporaryFirstDerivativeCompatiblePxi + temporaryFirstDerivativeUncompatiblePxi # dPx/dpi

      temporaryFirstDerivativeCompatibleLambda <- calculateFirstDerivativeByLambdaPDFPartial(compatibleObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) #dPx1/dlambda
      temporaryFirstDerivativeUncompatibleLambda <- calculateFirstDerivativeByLambdaPDF(compatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dlambda
      temporaryFirstDerivativeLambda <- temporaryFirstDerivativeCompatibleLambda + temporaryFirstDerivativeUncompatibleLambda # dPx/dlambda

      subsetCompatibleObservation <- which(observationsCompatibleWithFirstHaplotype)[-which(which(observationsCompatibleWithFirstHaplotype)==compatibleObservation)]
      temporaryFirstDerivativePyi <- 0
      for(subObservation in subsetCompatibleObservation){
        temporaryFirstDerivativePy <- calculateFirstDerivativeByFrequencyPDFPartial(haplotype, subObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) + calculateFirstDerivativeByFrequencyPDF(haplotype, subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype)
        temporaryFirstDerivativePyi <- temporaryFirstDerivativePyi + temporaryFirstDerivativePy
      }
      temp1 <- temp1 + temporaryFirstDerivativeLambda*temporaryFirstDerivativePyi
      temp2 <- temp2 + (sampleSize^2 - sampleSize + sampleSize/PDFPx)*temporaryFirstDerivativePxi*temporaryFirstDerivativeLambda

      tempPi <- 0
      tempLambda <- 0
      for(subObservation in which(!observationsCompatibleWithFirstHaplotype)){
        PDFPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
        PDFUncompatiblePy <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py2
        temporaryFirstDerivativeUncompatiblePyi <- calculateFirstDerivativeByFrequencyPDF(haplotype, subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPY(2)/dPi
        temporaryFirstDerivativeUncompatibleLambda <- calculateFirstDerivativeByLambdaPDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPY(2)/dlambda
        tempPi <- tempPi + temporaryFirstDerivativeUncompatiblePyi*PDFPy/PDFUncompatiblePy
        tempLambda <- tempLambda + temporaryFirstDerivativeUncompatibleLambda*PDFPy/PDFUncompatiblePy
      }
      temp3 <- temp3 + tempPi*temporaryFirstDerivativeLambda
      temp4 <- temp4 + tempLambda*temporaryFirstDerivativePxi
    }

    for(uncompatibleObservation in which(!observationsCompatibleWithFirstHaplotype)){
      PDFPx <- calculatePDF(uncompatibleObservation, subsetsFromObservations[[1]]) # Px
      PDFUncompatiblePx <- calculatePDF(uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Px2
      temporaryFirstDerivativeUncompatiblePxi   <- calculateFirstDerivativeByFrequencyPDF(haplotype, uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPi
      temporaryFirstDerivativeUncompatibleLambda <- calculateFirstDerivativeByLambdaPDF(uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPY(2)/dlambda

      tempi <- temporaryFirstDerivativeUncompatibleLambda/PDFUncompatiblePx
      subsetUncompatibleObservation <- which(!observationsCompatibleWithFirstHaplotype)[-which(which(!observationsCompatibleWithFirstHaplotype)==uncompatibleObservation)]
      temporaryi <- 0
      for(subObservation in subsetUncompatibleObservation){
        PDFPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
        PDFUncompatiblePy <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py2
        temporaryFirstDerivativeUncompatiblePyi   <- calculateFirstDerivativeByFrequencyPDF(haplotype, subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPy2/dPi
        temporaryi <- temporaryi + PDFPy*temporaryFirstDerivativeUncompatiblePyi/PDFUncompatiblePy
      }
      temp5 <- temp5 + tempi*PDFPx*temporaryi
      temp6 <- temp6 + (sampleSize^2 - sampleSize + sampleSize/PDFPx)*temporaryFirstDerivativeUncompatiblePxi*temporaryFirstDerivativeUncompatibleLambda*(PDFPx/PDFUncompatiblePx)^2
    }

    scoreMatrix[1,as.numeric(haplotype)] <- (sampleSize*(sampleSize - 1)*(temp1 + temp3 + temp4 + temp5) + temp2 + temp6)/calculateFirstDerivativeMeanMOIFunction(lambda)
    scoreMatrix[as.numeric(haplotype),1] <- scoreMatrix[1,as.numeric(haplotype)]
  }
  # Spi,pi
  for(haplotype in detectedHaplotypesLabels[-which(detectedHaplotypesLabels==1)]){
    temp1 <- 0
    temp2 <- 0
    temp3 <- 0
    temp4 <- 0
    temp5 <- 0

    for(compatibleObservation in which(observationsCompatibleWithFirstHaplotype)){
      #PDFCompatiblePx <- calculatePDFPartial(compatibleObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) # Px1
      #PDFUncompatiblePx <- calculatePDF(compatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Px2
      PDFPx <- calculatePDF(compatibleObservation, subsetsFromObservations[[1]]) # Px

      temporaryFirstDerivativeCompatiblePxi <- calculateFirstDerivativeByFrequencyPDFPartial(haplotype, compatibleObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) #dPx1/dPi
      temporaryFirstDerivativeUncompatiblePxi <- calculateFirstDerivativeByFrequencyPDF(haplotype, compatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPi
      temporaryFirstDerivativePxi <- temporaryFirstDerivativeCompatiblePxi + temporaryFirstDerivativeUncompatiblePxi # dPx/dpi

      subsetCompatibleObservation <- which(observationsCompatibleWithFirstHaplotype)[-which(which(observationsCompatibleWithFirstHaplotype)==compatibleObservation)]
      temporaryFirstDerivativePyi <- 0
      for(subObservation in subsetCompatibleObservation){
        temporaryFirstDerivativePy <- calculateFirstDerivativeByFrequencyPDFPartial(haplotype, subObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) + calculateFirstDerivativeByFrequencyPDF(haplotype, subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype)
        temporaryFirstDerivativePyi <- temporaryFirstDerivativePyi + temporaryFirstDerivativePy
      }
      temp1 <- temp1 + temporaryFirstDerivativePxi*temporaryFirstDerivativePyi
      temp2 <- temp2 + (sampleSize^2 - sampleSize + sampleSize/PDFPx)*(temporaryFirstDerivativePxi)^2

      tempPi <- 0
      for(subObservation in which(!observationsCompatibleWithFirstHaplotype)){
        PDFPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
        PDFUncompatiblePy <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py2
        temporaryFirstDerivativeUncompatiblePyi <- calculateFirstDerivativeByFrequencyPDF(haplotype, subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPY(2)/dPi
        tempPi <- tempPi + temporaryFirstDerivativeUncompatiblePyi*PDFPy/PDFUncompatiblePy
      }
      temp3 <- temp3 + tempPi*temporaryFirstDerivativePxi
    }

    for(uncompatibleObservation in which(!observationsCompatibleWithFirstHaplotype)){
      PDFPx <- calculatePDF(uncompatibleObservation, subsetsFromObservations[[1]]) # Px
      PDFUncompatiblePx <- calculatePDF(uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Px2
      temporaryFirstDerivativeUncompatiblePxi   <- calculateFirstDerivativeByFrequencyPDF(haplotype, uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPi

      tempPi <- temporaryFirstDerivativeUncompatiblePxi*PDFPx/PDFUncompatiblePx
      subsetUncompatibleObservation <- which(!observationsCompatibleWithFirstHaplotype)[-which(which(!observationsCompatibleWithFirstHaplotype)==uncompatibleObservation)]
      temporaryi <- 0
      for(subObservation in subsetUncompatibleObservation){
        PDFPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
        PDFUncompatiblePy <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py2
        temporaryFirstDerivativeUncompatiblePyi   <- calculateFirstDerivativeByFrequencyPDF(haplotype, subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPy2/dPi
        temporaryi <- temporaryi + PDFPy*temporaryFirstDerivativeUncompatiblePyi/PDFUncompatiblePy
      }
      temp4 <- temp4 + tempPi*temporaryi
      temp5 <- temp5 + (sampleSize^2 - sampleSize + sampleSize/PDFPx)*tempPi^2
    }

    scoreMatrix[as.numeric(haplotype),as.numeric(haplotype)] <- sampleSize*(sampleSize - 1)*(temp1 + 2*temp3 + temp4) + temp2 + temp5
  }

  # Spi,pj
  for(haplotypePairIndex in seq_len(nrow(haplotypePairs))){
    haplotypePair <- haplotypePairs[haplotypePairIndex,]

    temp1 <- 0
    temp2 <- 0
    temp3 <- 0
    temp4 <- 0
    temp5 <- 0
    temp6 <- 0

    for(compatibleObservation in which(observationsCompatibleWithFirstHaplotype)){
      PDFPx <- calculatePDF(compatibleObservation, subsetsFromObservations[[1]]) # Px

      temporaryFirstDerivativeCompatiblePxi <- calculateFirstDerivativeByFrequencyPDFPartial(haplotypePair[1], compatibleObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) #dPx1/dPi
      temporaryFirstDerivativeUncompatiblePxi <- calculateFirstDerivativeByFrequencyPDF(haplotypePair[1], compatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPi
      temporaryFirstDerivativePxi <- temporaryFirstDerivativeCompatiblePxi + temporaryFirstDerivativeUncompatiblePxi # dPx/dpi

      temporaryFirstDerivativeCompatiblePxj <- calculateFirstDerivativeByFrequencyPDFPartial(haplotypePair[2], compatibleObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) #dPx1/dPj
      temporaryFirstDerivativeUncompatiblePxj <- calculateFirstDerivativeByFrequencyPDF(haplotypePair[2], compatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPj
      temporaryFirstDerivativePxj <- temporaryFirstDerivativeCompatiblePxj + temporaryFirstDerivativeUncompatiblePxj # dPx/dpj

      subsetCompatibleObservation <- which(observationsCompatibleWithFirstHaplotype)[-which(which(observationsCompatibleWithFirstHaplotype)==compatibleObservation)]
      temporaryFirstDerivativePyj <- 0
      for(subObservation in subsetCompatibleObservation){
        temporaryFirstDerivativePy <- calculateFirstDerivativeByFrequencyPDFPartial(haplotypePair[2], subObservation, subsetsFromObservationsCompatibleWithFirstHaplotype) + calculateFirstDerivativeByFrequencyPDF(haplotypePair[2], subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype)
        temporaryFirstDerivativePyj <- temporaryFirstDerivativePyj + temporaryFirstDerivativePy
      }
      temp1 <- temp1 + temporaryFirstDerivativePxi*temporaryFirstDerivativePyj

      temp2 <- temp2 + (sampleSize^2 - sampleSize + sampleSize/PDFPx)*temporaryFirstDerivativePxi*temporaryFirstDerivativePxj

      tempPi <- 0
      tempPj <- 0
      for(subObservation in which(!observationsCompatibleWithFirstHaplotype)){
        PDFPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
        PDFUncompatiblePy <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py2
        temporaryFirstDerivativeUncompatiblePyi <- calculateFirstDerivativeByFrequencyPDF(haplotypePair[1], subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPY2/dPi
        temporaryFirstDerivativeUncompatiblePyj <- calculateFirstDerivativeByFrequencyPDF(haplotypePair[2], subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPY2/dPj
        tempPi <- tempPi + temporaryFirstDerivativeUncompatiblePyi*PDFPy/PDFUncompatiblePy
        tempPj <- tempPj + temporaryFirstDerivativeUncompatiblePyj*PDFPy/PDFUncompatiblePy
      }
      temp3 <- temp3 + temporaryFirstDerivativePxi*tempPj
      temp4 <- temp4 + temporaryFirstDerivativePxj*tempPi
    }

    for(uncompatibleObservation in which(!observationsCompatibleWithFirstHaplotype)){
      PDFPx <- calculatePDF(uncompatibleObservation, subsetsFromObservations[[1]]) # Px
      PDFUncompatiblePx <- calculatePDF(uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Px2
      temporaryFirstDerivativeUncompatiblePxi   <- calculateFirstDerivativeByFrequencyPDF(haplotypePair[1], uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPi
      temporaryFirstDerivativeUncompatiblePxj   <- calculateFirstDerivativeByFrequencyPDF(haplotypePair[2], uncompatibleObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPx2/dPj

      tempPi <- temporaryFirstDerivativeUncompatiblePxi/PDFUncompatiblePx
      subsetUncompatibleObservation <- which(!observationsCompatibleWithFirstHaplotype)[-which(which(!observationsCompatibleWithFirstHaplotype)==uncompatibleObservation)]
      tempPj <- 0
      for(subObservation in subsetUncompatibleObservation){
        PDFPy <- calculatePDF(subObservation, subsetsFromObservations[[1]]) # Py
        PDFUncompatiblePy <- calculatePDF(subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) # Py2
        temporaryFirstDerivativeUncompatiblePyj   <- calculateFirstDerivativeByFrequencyPDF(haplotypePair[2], subObservation, subsetsFromObservationsUncompatibleWithFirstHaplotype) #dPy2/dPj
        tempPj <- tempPj + PDFPy*temporaryFirstDerivativeUncompatiblePyj/PDFUncompatiblePy
      }
      temp5 <- temp5 + tempPi*PDFPx*tempPj
      temp6 <- temp6 + (sampleSize^2 - sampleSize + sampleSize/PDFPx)*temporaryFirstDerivativeUncompatiblePxi*temporaryFirstDerivativeUncompatiblePxj*(PDFPx/PDFUncompatiblePx)^2
    }

    scoreMatrix[as.numeric(haplotypePair[1]),as.numeric(haplotypePair[2])] <- sampleSize*(sampleSize - 1)*(temp1 + temp3 + temp4 + temp5) + temp2 + temp6
    scoreMatrix[as.numeric(haplotypePair[2]),as.numeric(haplotypePair[1])] <- scoreMatrix[as.numeric(haplotypePair[1]),as.numeric(haplotypePair[2])]
  }

  list(scoreMatrix, solve(scoreMatrix))
}

calculatePDFPartial <- function(compatibleObservation,subsetsFromObservation){
  probabilityDensityFunction <- 0
  for(subObservation in seq_len(subsetsFromObservation[[compatibleObservation]][[1]])){
    uncompatibleHaplotypes <- detectedHaplotypesLabels[-which(detectedHaplotypesLabels %in%unique(c(1, subsetsFromObservation[[compatibleObservation]][[4]][[subObservation]])))]
    sumFrequenciesOfCompatibleHaplotypes <- 1 - sum(haplotypeFrequencies[as.numeric(uncompatibleHaplotypes),])
    inclusionExclusionCoefficient <- subsetsFromObservation[[compatibleObservation]][[3]][[subObservation]]
    probabilityDensityFunction <- probabilityDensityFunction + inclusionExclusionCoefficient*calculateGeneratingFunction(lambda, sumFrequenciesOfCompatibleHaplotypes)
  }
  probabilityDensityFunction
}

calculatePDF <- function(compatibleObservation,subsetsFromObservation){
  probabilityDensityFunction <- 0
  for(subObservation in seq_len(subsetsFromObservation[[compatibleObservation]][[1]])){
    sumFrequenciesOfCompatibleHaplotypes <- sum(haplotypeFrequencies[as.numeric(subsetsFromObservation[[compatibleObservation]][[4]][[subObservation]]),])
    inclusionExclusionCoefficient <- subsetsFromObservation[[compatibleObservation]][[3]][[subObservation]]
    probabilityDensityFunction <- probabilityDensityFunction + inclusionExclusionCoefficient*calculateGeneratingFunction(lambda, sumFrequenciesOfCompatibleHaplotypes)
  }
  probabilityDensityFunction
}

calculateFirstDerivativeByLambdaPDFPartial <- function(compatibleObservation,subsetsFromObservations){
  probabilityDensityFunction <- 0
  for(subObservation in seq_len(subsetsFromObservations[[compatibleObservation]][[1]])){
    uncompatibleHaplotypes <- detectedHaplotypesLabels[-which(detectedHaplotypesLabels %in%unique(c(1, subsetsFromObservations[[compatibleObservation]][[4]][[subObservation]])))]
    sumFrequenciesOfCompatibleHaplotypes <- 1 - sum(haplotypeFrequencies[as.numeric(uncompatibleHaplotypes),])
    inclusionExclusionCoefficient <- subsetsFromObservations[[compatibleObservation]][[3]][[subObservation]]
    probabilityDensityFunction <- probabilityDensityFunction + inclusionExclusionCoefficient*calculateFirstDerivativeGeneratingFunction(lambda, sumFrequenciesOfCompatibleHaplotypes)
  }
  probabilityDensityFunction
}

calculateFirstDerivativeByLambdaPDF <- function(compatibleObservation,subsetsOfObservations){
  probabilityDensityFunction <- 0
  for(subObservation in seq_len(subsetsOfObservations[[compatibleObservation]][[1]])){
    sumFrequenciesOfCompatibleHaplotypes <- sum(haplotypeFrequencies[as.numeric(subsetsOfObservations[[compatibleObservation]][[4]][[subObservation]]),])
    inclusionExclusionCoefficient <- subsetsOfObservations[[compatibleObservation]][[3]][[subObservation]]
    probabilityDensityFunction <- probabilityDensityFunction + inclusionExclusionCoefficient*calculateFirstDerivativeGeneratingFunction(lambda, sumFrequenciesOfCompatibleHaplotypes)
  }
  probabilityDensityFunction
}

calculateFirstDerivativeByFrequencyPDFPartial <- function(haplotype, compatibleObservation, subsetsFromObservations){
  probabilityDensityFunction <- 0
  for(subObservation in seq_len(subsetsFromObservations[[compatibleObservation]][[1]])){
    uncompatibleHaplotypes <- detectedHaplotypesLabels[-which(detectedHaplotypesLabels %in%unique(c(1, subsetsFromObservations[[compatibleObservation]][[4]][[subObservation]])))]
    sumFrequenciesOfUncompatibleHaplotypes <- 1 - sum(haplotypeFrequencies[as.numeric(uncompatibleHaplotypes),])
    haplotypeInSubset <- subsetsFromObservations[[compatibleObservation]][[4]][[subObservation]]
    isHaplotypeNotInSubset <- !(haplotype %in% haplotypeInSubset)
    isHaplotypeNotFirstHaplotype <- !(haplotype == 1)
    inclusionExclusionCoefficient <- subsetsFromObservations[[compatibleObservation]][[3]][[subObservation]]
    probabilityDensityFunction <- probabilityDensityFunction - inclusionExclusionCoefficient*calculateFirstDerivativeGeneratingFunctionByFrequency(lambda, sumFrequenciesOfUncompatibleHaplotypes)*isHaplotypeNotInSubset*isHaplotypeNotFirstHaplotype
  }
  probabilityDensityFunction
}

calculateFirstDerivativeByFrequencyPDF <- function(haplotype, compatibleObservation, subsetsFromObservations){
  probabilityDensityFunction <- 0
  for(subObservation in seq_len(subsetsFromObservations[[compatibleObservation]][[1]])){
    haplotypeInSubset <- subsetsFromObservations[[compatibleObservation]][[4]][[subObservation]]
    isHaplotypePresent <- haplotype %in% haplotypeInSubset
    sumFrequenciesOfCompatibleHaplotypes <- sum(haplotypeFrequencies[as.numeric(haplotypeInSubset),])
    inclusionExclusionCoefficient <- subsetsFromObservations[[compatibleObservation]][[3]][[subObservation]]
    probabilityDensityFunction <- probabilityDensityFunction + inclusionExclusionCoefficient*calculateFirstDerivativeGeneratingFunctionByFrequency(lambda, sumFrequenciesOfCompatibleHaplotypes)*isHaplotypePresent
  }
  probabilityDensityFunction
}

findObservationsCompatibleWithFirstHaplotype <- function(numberOfAllelesAtEachMarker, possibleObservations){
  alleleResearched <- 1
  numberOfLoci <- length(numberOfAllelesAtEachMarker)
  presence <- vector(,nrow(possibleObservations))
  for(locus in seq_len(numberOfLoci)){
    presence <- presence + mapply(function(x){as.numeric(intToBits(x)[1:numberOfAllelesAtEachMarker[locus]])[alleleResearched]},
           possibleObservations[,locus])
  }
  presence >= length(numberOfAllelesAtEachMarker)
}

splitSubsetsByCompabilityWithFirstHaplotype <- function(numberOfAllelesAtEachMarker, subsetsFromObservations){
  totalNumberOfPossibleObservations <- prod(2^numberOfAllelesAtEachMarker-1)
  compatibleSubObservations <- vector("list", totalNumberOfPossibleObservations)
  for(observation in seq_len(totalNumberOfPossibleObservations)){
    compatibleSubObservations[[observation]] <- findObservationsCompatibleWithFirstHaplotype(numberOfAllelesAtEachMarker, subsetsFromObservations[[1]][[observation]][[2]])
  }

  subsetsFromObservationsCompatibleWithFirstHaplotype <- vector("list", nrow(possibleObservations))
  subsetsFromObservationsUncompatibleWithFirstHaplotype <- vector("list", nrow(possibleObservations))
  for(observation in seq_len(nrow(possibleObservations))){
    temporaryCompatibleSubset <- subsetsFromObservations[[1]][[observation]]
    temporaryCompatibleSubset[[1]] <- sum(compatibleSubObservations[[observation]])
    temporaryCompatibleSubset[[2]] <- matrix(temporaryCompatibleSubset[[2]][compatibleSubObservations[[observation]],], nrow = temporaryCompatibleSubset[[1]])
    temporaryUncompatibleSubset <- subsetsFromObservations[[1]][[observation]]
    temporaryUncompatibleSubset[[1]] <- sum(!compatibleSubObservations[[observation]])
    temporaryUncompatibleSubset[[2]] <- temporaryUncompatibleSubset[[2]][!compatibleSubObservations[[observation]],]

    for(i in 3:4){
      temporaryCompatibleSubset[[i]] <- temporaryCompatibleSubset[[i]][compatibleSubObservations[[observation]]]
      temporaryUncompatibleSubset[[i]] <- temporaryUncompatibleSubset[[i]][!compatibleSubObservations[[observation]]]
    }
    subsetsFromObservationsCompatibleWithFirstHaplotype[[observation]] <- temporaryCompatibleSubset
    subsetsFromObservationsUncompatibleWithFirstHaplotype[[observation]] <- temporaryUncompatibleSubset
  }
  list(subsetsFromObservationsCompatibleWithFirstHaplotype, subsetsFromObservationsUncompatibleWithFirstHaplotype)
}

calculateProbabilityDesityFunction <- function(){
  for(subObservation in seq_len(subsetsFromObservationsCompatibleWithFirstHaplotype[[compatibleObservation]][[1]])){
    sump <- sum(pp[Ax[[u]][[4]][[k]],])
    vz   <- Ax[[u]][[3]][[k]]
    out1 <- out1 + vz*calculateGeneratingFunction(lambda,sump)   # calculateGeneratingFunction to build Px
  }
}

# Generating function
calculateGeneratingFunction <- function(lambda, haplotypeFrequencies){
  elp   <- exp(lambda*haplotypeFrequencies)
  elmo  <- exp(lambda) - 1
  (elp - 1)/elmo
}

# First derivative of generating function with respect to lambda
calculateFirstDerivativeGeneratingFunction <- function(lambda, haplotypeFrequencies){
  el <- exp(lambda)
  elmo <- exp(lambda)-1
  elp  <- exp(lambda*haplotypeFrequencies)

  haplotypeFrequencies*elp/elmo - el*(elp - 1)/(elmo^2)
}

calculateFirstDerivativeGeneratingFunctionByFrequency <- function(lambda, haplotypeFrequencies){
  elmo <- exp(lambda)-1
  lambda*exp(haplotypeFrequencies)/elmo
}

# Second derivative of generating function with respect to lambda
calculateSecondDerivativeGeneratingFunction <- function(lambda, haplotypeFrequencies){
  el   <- exp(lambda)
  elmo <- exp(lambda)-1
  elp  <- exp(lambda*haplotypeFrequencies)

  (elp*haplotypeFrequencies^2)/elmo - 2*haplotypeFrequencies*exp(lambda*(haplotypeFrequencies+1))/(elmo^2) - el*(elp-1)/(elmo^2) + 2*exp(2*lambda)*(elp-1)/(elmo^3)
}

# First derivative of mean MOI
calculateFirstDerivativeMeanMOIFunction <- function(lambda){
  eml <- 1-exp(-lambda)
  1/eml - (lambda*exp(-lambda))/(eml^2)
}

# Estimation of CRLB
cramerRaoLowerBound <- function(haplotypeFrequencies,lambda,sampleSize,numberOfAllelesAtEachMarker){
  #### Generate all possible observations given numberOfAllelesAtEachMarker
  detectedObservations    <- buildAllPossibleObservations(numberOfAllelesAtEachMarker)
  Nobs <- nrow(detectedObservations)

  # Find for each infection detectedObservations all components necessary for the computations
  Ax <- buildAllSetsAndSubsets(detectedObservations, numberOfAllelesAtEachMarker)[[1]]

  # Initialize the Fisher information matrix of degree d
  names(haplotypeFrequencies) <- seq_along(haplotypeFrequencies)
  pp   <- matrix(haplotypeFrequencies, ncol=1)
  hap  <- as.numeric(names(haplotypeFrequencies))
  #happ <- t(combn(hap,2))
  rownames(pp) <- hap

  d <- nrow(pp)+2
  diag.el <- diag(array(1:d^2,c(d,d)))[-1]
  diag.el <- diag.el[-length(diag.el)]

  I <- matrix(0, nrow=d, ncol=d) # d*(d+1)/2 dof
  rownames(I) <- 1:d
  colnames(I) <- 1:d

  elmo <- exp(lambda)-1

  # Ipsi,psi
  out <- 0
  for(u in 1:Nobs){ # For each observation x in ScrO
    out1 <- 0
    out2 <- 0
    out3 <- 0
    for(k in 1:(Ax[[u]][[1]][[1]])){  # For each observation y in the sub-observation ScrAx
      sump <- sum(pp[Ax[[u]][[4]][[k]],])
      vz   <- Ax[[u]][[3]][[k]]
      out1 <- out1 + vz*calculateGeneratingFunction(lambda,sump)   # calculateGeneratingFunction to build Px
      out2 <- out2 + vz*calculateFirstDerivativeGeneratingFunction(lambda,sump)  # dG/dl to build dPx/dlam
      out3 <- out3 + vz*calculateSecondDerivativeGeneratingFunction(lambda,sump) # d2G/dl to build d^2Px/dlam^2
    }
    out <- out + (out2^2)/out1 - out3 #Ill
  }
  I[1,1] <- sampleSize*out/(calculateFirstDerivativeMeanMOIFunction(lambda)^2)   # I_lam_lam

  # Ipi,pi
  out <- haplotypeFrequencies*0
  for(u in 1:Nobs){
    out1 <- 0
    out2 <- haplotypeFrequencies*0
    out3 <- haplotypeFrequencies*0
    for(k in 1:Ax[[u]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[Ax[[u]][[4]][[k]],])
      elp  <- exp(lambda*sump)
      vz   <- Ax[[u]][[3]][[k]]
      tmp  <- vz*(lambda*elp/elmo)
      out1 <- out1 + vz*calculateGeneratingFunction(lambda,sump)                                      # calculateGeneratingFunction
      out2[as.numeric(Ax[[u]][[4]][[k]])] <- out2[as.numeric(Ax[[u]][[4]][[k]])] + tmp          # dG/dPi * Indic
      out3[as.numeric(Ax[[u]][[4]][[k]])] <- out3[as.numeric(Ax[[u]][[4]][[k]])] + lambda*tmp       # d2G/dpi
    }
    out <- out + (out2^2)/out1 - out3
  }
  I[diag.el] <- sampleSize*out

  # Ipi,pj  new
  out <- 0 * (haplotypeFrequencies%*%t(haplotypeFrequencies))
  for(u in 1:Nobs){
    out1  <- 0
    out2 <- 0*haplotypeFrequencies
    out3  <- 0* haplotypeFrequencies %*% t(haplotypeFrequencies)
    rownames(out3) <- colnames(out3)
    for(k in 1:Ax[[u]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[Ax[[u]][[4]][[k]],])
      elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      tmp   <- vz*lambda*elp/elmo
      out1  <- out1 + vz*calculateGeneratingFunction(lambda,sump)
      out2[as.numeric(Ax[[u]][[4]][[k]])] <- out2[as.numeric(Ax[[u]][[4]][[k]])] + tmp          # dG/dpi * Indic
      out3[as.numeric(Ax[[u]][[4]][[k]]),as.numeric(Ax[[u]][[4]][[k]])]  <- out3[as.numeric(Ax[[u]][[4]][[k]]),as.numeric(Ax[[u]][[4]][[k]])] + lambda*tmp   # d2G/dpidpj
    }
    out <- out + (out2 %*% t(out2/out1)) - out3
  }
  I[2:(d-1),2:(d-1)] <- sampleSize*out

  # Ipsi,pi
  out <- 0*haplotypeFrequencies
  for(u in 1:Nobs){
    out1  <- 0
    out21 <- 0
    out22 <- 0*haplotypeFrequencies
    out3  <- 0*haplotypeFrequencies
    for(k in 1:Ax[[u]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[Ax[[u]][[4]][[k]],])
      elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1  + vz*calculateGeneratingFunction(lambda,sump)             # calculateGeneratingFunction
      out21 <- out21 + vz*calculateFirstDerivativeGeneratingFunction(lambda,sump)            # dG/dl
      out22[as.numeric(Ax[[u]][[4]][[k]])] <- out22[as.numeric(Ax[[u]][[4]][[k]])] + vz*(lambda*elp/elmo)    # dG/dpi
      out3[as.numeric(Ax[[u]][[4]][[k]])]  <- out3[as.numeric(Ax[[u]][[4]][[k]])]  + vz*((lambda*sump + 1)*elp/elmo - lambda*exp(lambda*(sump+1))/(elmo^2))  # d2G/dpidl
    }
    out <- out + (out21*out22)/out1 - out3
  }
  tmpI <- sampleSize*out/calculateFirstDerivativeMeanMOIFunction(lambda)
  I[1,2:(d-1)] <- tmpI
  I[2:(d-1),1] <- tmpI

  # Parameter space in Higher dimension
  dbeta <- c(0,rep(1,prod(numberOfAllelesAtEachMarker)),0)
  I[,d] <- dbeta
  I[d,] <- dbeta

  print("Inverting Now...")
  #list(I,matrix.inverse(I))
  out <- solve(I)
  out <- out[-d,]
  out <- out[,-d]
  I <- I[-d,]
  I <- I[,-d]
  list(I,out)
}

heterozygosity <- function(data1.tranf, snps.range, STR.markers, freq.threshold){
  # Estimating haplotype frequencies from transformed data
  DATA         <- data1.tranf[[1]][,snps.range]       # Data in range of markers of interest
  numberOfAllelesAtEachMarker         <- data1.tranf[[3]][snps.range]        # Genetic architecture in range of interest
  list.alleles <- data1.tranf[[2]][snps.range]        # Detected alleles
  data1 <- convertToListOfDatasets(DATA, idExists = FALSE)
  n <- sum(data1[[2]])

  estimate <- calculateMaximuLikelihoodEstimatesWithBiasCorrection(data1, numberOfAllelesAtEachMarker)
  freq.est <- estimate$p
  hapl.freqs <- nameHaplotypes(freq.est, list.alleles)

  # Heterozygosity among drug-resistant haplotypes with frequency > freq.threshold
  # Pick the right frequency threshold that excludes the wildtype.
  hapllist <- rownames(hapl.freqs)[hapl.freqs > freq.threshold]

  # Computing heterozygosity for each microsatelitte marker considering the haplotypes identified above
  out.het <- array(,c(length(hapllist),length(STR.markers)))  # matrix to store the heterozygosity values rows - haplotype background, columns are the microsatellite markers
  rownames(out.het) <- hapllist

  # Heterozygosity at each microsatellite marker
  for(kk in seq_along(STR.markers)){
    k            <- STR.markers[kk]                          # Pick microsatellite marker
    DATA         <- data1.tranf[[1]][,c(snps.range,k)]       # Data in range of markers of interest
    numberOfAllelesAtEachMarker         <- data1.tranf[[3]][c(snps.range,k)]        # Genetic architecture in range of interest
    list.alleles <- data1.tranf[[2]][c(snps.range,k)]        # Detected alleles
    data1 <- convertToListOfDatasets(DATA, idExists = FALSE)
    n <- sum(data1[[2]])

    estimate <- calculateMaximuLikelihoodEstimatesWithBiasCorrection(data1, numberOfAllelesAtEachMarker)
    freq.est <- estimate$p

    freq.est <- nameHaplotypes(freq.est,list.alleles) ## frequency estimates with proper haplotype names

    # next remove the microsatellite from the haplotype names
    names <- t(sapply(rownames(freq.est), function(x){
      y <- unlist(strsplit(x,"-"))
      c(paste(y[-length(y)],collapse="-"),y[length(y)])
    }
    ))

    freq.est1 <- as.data.frame(names)  # haplotype freqeuncy estimates
    freq.est1$freq <- freq.est

    # Heterozygosity at locus k is calculated from haplotype frequencies of
    # haplotypes containing the drug resistant alleles (conditional heterozygosity).
    for(h in hapllist){
      # Frequencies of alleles at locus k conditionned on drug-resistant haplotype h
      pp <- freq.est1[freq.est1[, 1]==h, 3]
      if(sum(pp) == 0){
        pp <- 0
      }else {
        pp <- pp/sum(pp)
      }
      out.het[h, kk] <- (1-sum(pp^2))*(n/(n-1))
    }
  }
  list(out.het, hapllist)
}

### function data.format outputs a list containing as first element  the data in new format (1 row per sample, 1 colum by marker) 
### entries are integerers. If transformed into binary numbers 0-1 vectors they indicate absence/presence of alleles 
### e.g., 0.... no allele present, 1... first allele present, 2... 2nd allele present, 3.... 1st and 2nd allele present ect. 
### The second element gives the order of the alleles pper locus
### The third element gives th number of alleles per locucs
data.format <- function(dat,markers){
    ### dat... is the input data set in standard format of package MLMOI, 1 st column contains smaple IDs
    ### markers ... vector of columms containing markers to be included

    ### list of alleles per marker###############
    allele.list <- apply(dat[,markers],2, function(x){ 
                                                        y=sort(unique(x))
                                                        y[!is.na(y)]
                                                    } 
                        )
    ###number of alleles per marker########
    allele.num <- unlist(lapply(allele.list,length))  
    #### split data by sample ID
    dat.split <- split(dat[,markers],dat[,1])
    #### Binary representation of allele being absent and present
    samples.coded <- t(sapply(dat.split, function(x){mapply(function(x,y,z){as.integer(is.element(y,x)) %*% 2^(0:(z-1))}, x , allele.list, allele.num)}))
    list(samples.coded,allele.list,allele.num)
}

### This function outputs a list
#### 1st gives the compact notation for the data
#### 2nd element base factos 1 g1 g1g2. ....,g1*...g(l-1)
#### 3rd element has the number of alleles in the prope order per locus
#### 4th element is the number of alleles per locus
dafa.formal.AL <- function(data){
        l <- length(data[[3]])   ## number of loci
        basefactors <- cumprod(c(1,2^(data[[3]][-l])))
        list(data[[1]] %*% basefactors, basefactors,data[[2]],data[[3]])
}

### Going back

div <- function(X,f){
    x <- NULL
    for(f in rev(dat[[2]])){
        #print(c(X,f))
        x <- c(X %/% f,x) 
        X <- X%%f
    }
    x
}

reconst.data <- function(data){
    newdat <- t(sapply(data[[1]], function(x) div(x,data[[2]]) ) )
    list(newdat,data[[3]],data[[4]])
}

# help functions 

varsets2 <- function(l){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  n <- length(l)
  B <- array(1,c(prod(l),n))
  B[1:l[1],1] <- 1:l[1]
  lkmo <- l[1]
  if(n>1){
    for(k in 2:n){
      if(l[k]>1){
        lk <- lkmo*l[k]
        pick1 <- (lkmo+1):lk
        B[pick1,] <- B[rep(1:lkmo,l[k]-1),]
        B[pick1,k] <- rep(2:l[k],each=lkmo)
        lkmo <- lk
      }
                              
    }
  }
  B
}

varsets1 <- function(l){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  n <- length(l)
  B <- array(0,c(prod(l),n))
  B[1:l[1],1] <- 0:(l[1]-1)
  lkmo <- l[1]
  if(n>1){
    for(k in 2:n){
      if(l[k]>1){
        lk <- lkmo*l[k]
        pick1 <- (lkmo+1):lk
        B[pick1,] <- B[rep(1:lkmo,l[k]-1),]
        B[pick1,k] <- rep(1:(l[k]-1),each=lkmo)
        lkmo <- lk   
      }
                           
    }
  }
  B
}

varsets <- function(l,n){   #calculate all var sets
  # n number of loci
  # l number of alleles per locus
  B <- array(0,c(l^n,n))
  B[1:l,1] <- 0:(l-1)
  lkmo <- l
  if(n>1){
    for(k in 2:n){
      lk <- lkmo*l
      pick1 <- (lkmo+1):lk
      B[pick1,] <- B[rep(1:lkmo,l-1),]
      B[pick1,k] <- rep(1:(l-1),each=lkmo)
      lkmo <- lk                        
    }
  }
  B
}

gead <- function(x,l,n){   ## calculates geadic expression of each element of vectorx 
  
  l <- rep(l,n)

  out <- array(0,c(length(x),n))
  div <- c(1,cumprod(l[1:(n-1)]))
  for(k in n:1){
    r <- x%%div[k]
    #print(c(r,div[k]))
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}

gead1 <- function(x,l){   ## calculates general geadic expression of each element of vector x
  n <- length(l)
  out <- array(0,c(length(x),n))
  div <- c(1,cumprod(l[1:(n-1)]))
  for(k in n:1){
    r <- x%%div[k]
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}

# main function
est <- function(X,Nx,l){
  eps <- 10^-8
  N <- sum(Nx)
  nn <- nrow(X) 
  n <- ncol(X)
  x <- X
  
  #_____________________________
  # this calculates the list to pick proper haplotype freuencies, i.e. sets Ay
  
  #allconf <- function(x,l,n){  # x array
  hapll <- list()
  if(length(l)==1){
    l <- rep(l,n)
  }else{
    n <- length(l)
  }
  ggead <- c(1,cumprod(l[1:(n-1)]))
  Hx <- list()
  Ax <- list()
  alnum <- list()
  bin2num <- list()
  for(k in 1:n){
    alnum[[k]] <- 1:l[[k]]
    bin2num[[k]] <- 2^(0:(l[[k]]-1))
  }
  alcnt <- array(0,n)
  for(u in 1:nn){
    Hx[[u]] <- list(array(0,n),list(),list(),list(),list(),array(0,n))
    Ax[[u]] <- list(list(),list(),list(),list())
    for(k in 1:n){
      temp <- gead(x[u,k],2,l[[k]])
      temp1 <- varsets1(temp+1)[-1,]
      Hx[[u]][[1]][k] <- sum(temp)
      Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
      Hx[[u]][[3]][[k]] <- temp1
      Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] # subsets of alleles
      Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,l[k]) # number of alleles in subset at locus k
      Hx[[u]][[6]][k] <- length(temp1%*%bin2num[[k]])
    }
    vz1 <- sum(Hx[[u]][[1]])
    temp2 <- prod(Hx[[u]][[6]])
    Ax[[u]][[1]] <- temp2
    Ax[[u]][[2]] <- varsets2(Hx[[u]][[6]])
    for(k in 1:n){
      Ax[[u]][[2]][,k] <- Hx[[u]][[4]][[k]][Ax[[u]][[2]][,k]]
    }
    for(j in 1:temp2){
      Ax[[u]][[3]][[j]] <- list() 
      for(k in 1:n){
        temp <- gead(Ax[[u]][[2]][j,k],2,l[[k]])
        temp1 <- (alnum[[k]])[temp*alnum[[k]]]
        alcnt[k] <- length(temp1)
        Ax[[u]][[3]][[j]][[k]] <- temp1
      }
      Ax[[u]][[4]][[j]] <- varsets2(alcnt)
      for(k in 1:n){
        Ax[[u]][[4]][[j]][,k] <- Ax[[u]][[3]][[j]][[k]][Ax[[u]][[4]][[j]][,k]]
      }    
      Ax[[u]][[4]][[j]] <- as.character((Ax[[u]][[4]][[j]]-1)%*%ggead+1)
      Ax[[u]][[3]][[j]] <- (-1)^(vz1+sum(alcnt))
    }
    hapll[[u]] <- Ax[[u]][[4]][[temp2]]
    
  }  
  
  hapl1 <- unique(unlist(hapll))
  
  #---------------------------------------
  # initialize parameters
  H <- length(hapl1)
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1
  
  #initial list#
  num0 <- pp*0
  cond1 <- 1  ## condition to stop EM alg! 
  lambda <- 2
  num <- num0
  rownames(num) <- hapl1
  Bcoeff <- num0
  rownames(Bcoeff) <- hapl1
  
  while(cond1>eps){
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    num <- num0  #reset numerator to 0 in next iteration
    for(u in 1:nn){
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[1]]){
        haplotypeFrequencies <- sum(pp[Ax[[u]][[4]][[k]],])
        vz <- Ax[[u]][[3]][[k]]
        lap <- lambda*haplotypeFrequencies
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz  ##   = (1-)^(Nx-Ny)*(Exp(lambda*sum haplotypeFrequencies)-1) = (1-)^(Nx-Ny)*calculateGeneratingFunction(sum haplotypeFrequencies)
        num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],]+ exlap#*pp[Ax[[u]][[1]][[k]],]
        ## exlap =  (1-)^(Nx-Ny) G'(sum haplotypeFrequencies)   --- denominator of generating functions cancels out!
        CC <- CC + exlap*haplotypeFrequencies
      }
      num <- num*pp
      denom <- Nx[u]/denom
      denom <- lambda*denom
      Ccoeff <- Ccoeff + CC*denom
      Bcoeff <- Bcoeff + num*denom
    }
    Ccoeff <- Ccoeff/N
    ppn <- Bcoeff/(sum(Bcoeff))
    
    ### Newton step
    cond2 <- 1
    xt <- Ccoeff   ### good initial condition
    while(cond2 > eps){
      ex <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      cond2 <- abs(xtn-xt)
      xt <- xtn
    }
    cond1 <- abs(xt-lambda)+sqrt(sum((pp-ppn)^2))
    lambda <- xt
    pp <- ppn
    #print("________")
    #print(cond1)
    #print(c(lambda,pp))
    
  }
  list(lambda,pp)
}

# likelihood function
likeGen <- function(pp,lambda,Nx,N,Ax){
  logli <- 0
  num0 <- 0
  Bcoeff <- num0 #reset B coefficients to 0 in next iteration
  num <- num0  #reset numerator to 0 in next iteration
  nn <- (length(Ax))
  for(u in 1:nn){
    denom <- 0
    num <- num0
    CC <- 0
    for(k in 1:Ax[[u]][[1]]){
      haplotypeFrequencies <- sum(pp[Ax[[u]][[4]][[k]],])
      vz <- Ax[[u]][[3]][[k]]
      lap <- lambda*haplotypeFrequencies
      exlap <- vz*exp(lap)
      denom <- denom + (exlap-vz) 
    }
    
    logli <-  logli+ Nx[u]*log(denom/(exp(lambda)-1))
  }
 logli 
}  

# labels output frequencies correctly
hapl.names <- function(pp,allist){ ## allist list with alleles per locos
  n <- length(allist)
  allnum <- unlist(lapply(allist,length))
  hapl <- gead1(as.numeric(rownames(pp))-1,allnum)+1
  hapl1 <- array(,dim(hapl))
  for(l in 1: ncol(hapl)){
    for(m in 1: nrow(hapl)){
      hapl1[m,l] <- allist[[l]][hapl[m,l]]
    }
  }
  hapl1 <- apply(hapl1,1,function(x) paste(x,sep="",collapse="-"))
  rownames(pp) <- hapl1
  pp
}

# plot of frequencies
freqplot <- function(out,out1,lv){
  al <- attr(out[[5]],"names") #alleles in 1st data set
  al1 <- attr(out1[[5]],"names") #alleles in 2st data set
  alleles <- unique(c(al,al1)) # alleles in both data sets
  
  pldata <- array(0,c(length(alleles),2)) # create array that contains frequencies
  
  rownames(pldata) <- alleles # the rows are labelled with the STR repet lengths
  colnames(pldata) <- lv  # the colums are the year intervals
  pldata[al,1] <- out[[3]]  # in the first column the alllele frequencies of the old data
  pldata[al1,2] <- out1[[3]]
  
  pldata <- melt(pldata,id.vars="all")  # reshapes the data 3 colums Var1, Var2, value
  
  # next transform plotdata into data frame
  pldata <- data.frame(allele=as.factor(pldata$Var1),gr=as.factor(pldata$Var2),freq=as.numeric(pldata$value))
  
  # Define colors for the plot - always be coulorbliend friendly
  #cbPalette <- c("#3ffdcc","#fa6403","#cdff00","#1e9864","#999999", "#E69F00")
  cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
  #plot  
  p <- ggplot(data=pldata, aes(x=allele, y=freq, fill=gr)) +
      geom_bar(stat="identity",position=position_dodge(), colour="black") 
  p <- p + theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(colour='black',fill=NA),
                 axis.text = element_text(size = rel(1.3),color='black'),
                 axis.text.x = element_text(angle=75,vjust=0.5),
                 axis.title = element_text(size = rel(1.3)),
                 plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
                 legend.text = element_text(size = rel(1.3)),
                 legend.title = element_text(size = rel(1.3)))
  p <- p + scale_fill_manual(values=cbPalette,name="Years") + ylim(0,1)
  p <- p + labs(x="alleles",y="frequencies",title=parse(text=out[[6]]))
  p
}
