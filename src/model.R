# Title        : Maximum likelihood method
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 15.08.23
# Last modified: 01.11.23

calculateMaximumLikelihoodEstimatesWithAddOns <- function(dataset, numberOfAllelesAtEachMarker, idExists=TRUE, pluginValueOfLambda=NULL, isConfidenceInterval=FALSE, isBiasCorrection=FALSE, methodForBiasCorrection="bootstrap", numberOfBootstrapReplicatesBiasCorrection=10000, numberOfBootstrapReplicatesConfidenceInterval=10000, significanceLevel=0.05){
  ### Dropping Missing data in dataset
  samplesWithMissingData <- rowSums(dataset == 0) > 0
  dataset <- dataset[!samplesWithMissingData,]

  ### Actual sample size
  sampleSizeWithNoMissingData <- nrow(dataset)
  listOfDatasets  <- convertToListOfDatasets(dataset, idExists=idExists)
  detectedObservations <<- listOfDatasets[[1]]
  numberOfEachDetectedObservations <- listOfDatasets[[2]]
  numberOfLoci <- ncol(detectedObservations)

  # MLEs
  maximumLikelihoodEstimates <- calculateMaximuLikelihoodEstimatesWithBiasCorrection(listOfDatasets, numberOfAllelesAtEachMarker, isBiasCorrection=isBiasCorrection, methodForBiasCorrection=methodForBiasCorrection, numberOfBootstrapReplicatesBiasCorrection=numberOfBootstrapReplicatesBiasCorrection, pluginValueOfLambda=pluginValueOfLambda)
  mixRadixBaseCumulativeProduct <- c(1, cumprod(numberOfAllelesAtEachMarker)[1:(numberOfLoci-1)])
  frequenciesEstimates <- maximumLikelihoodEstimates[[2]]
  lambdaEstimates <- maximumLikelihoodEstimates[[1]]
  rnames1 <- as.integer(rownames(frequenciesEstimates)) - 1
  rnames  <- rnames1
  nh <- length(rnames)
  mixRadixTableOfDetectedHaplotypes <- array(0,c(nh,numberOfLoci))
  for(k in 1:numberOfLoci){ #for each locus
    re <- rnames%%rev(mixRadixBaseCumulativeProduct)[k]
    mixRadixTableOfDetectedHaplotypes[,(numberOfLoci-k+1)] <- (rnames-re)/rev(mixRadixBaseCumulativeProduct)[k]
    rnames <- re
  }
  mixRadixTableOfDetectedHaplotypes <- mixRadixTableOfDetectedHaplotypes+1
  for(i in 1:nh){
    rnames[i] <- paste(mixRadixTableOfDetectedHaplotypes[i,], collapse = '')
  }
  rownames(frequenciesEstimates) <- rnames

  # Bootstrap CIs
  if(isConfidenceInterval){
    numberOfDetectedHaplotypes  <- length(frequenciesEstimates)
    sampleSize <<- sum(numberOfEachDetectedObservations)
    probabilityOfEachObservation  <<- numberOfEachDetectedObservations/sampleSize
    arrayOfBootstrappedEstimates <- array(0, dim = c((numberOfDetectedHaplotypes+1), numberOfBootstrapReplicatesConfidenceInterval=numberOfBootstrapReplicatesConfidenceInterval))
    rownames(arrayOfBootstrappedEstimates) <- c('lambda',(rnames1+1))
    observation <<- seq_along(numberOfEachDetectedObservations)
    for (bootstrapReplicate in 1:numberOfBootstrapReplicatesConfidenceInterval){
      bootstrappedListOfDatasets <- buildBootstrappedDataset()
      bootstrappedEstimates <- calculateMaximumLikelihoodEstimatesWithOrWithoutPlugin(bootstrappedListOfDatasets, numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda)
      detectedHaplotypes      <- as.integer(rownames(bootstrappedEstimates[[2]]))
      arrayOfBootstrappedEstimates[1,bootstrapReplicate]  <- unlist(bootstrappedEstimates[[1]])
      arrayOfBootstrappedEstimates[as.character(detectedHaplotypes),bootstrapReplicate] <- unlist(bootstrappedEstimates[[2]])
    }
    perc <- t(apply(arrayOfBootstrappedEstimates, 1, quantile, c(significanceLevel/2, (1-significanceLevel/2))))
    if(is.null(pluginValueOfLambda)){
      lambdaEstimates <- c(unlist(maximumLikelihoodEstimates[[1]]), perc[1,])
      names(lambdaEstimates) <- c('', paste0(as.character((significanceLevel/2)*100), '%'), paste0(as.character((1-significanceLevel/2)*100), '%'))
    }else{
      lambdaEstimates <- maximumLikelihoodEstimates[[1]]
      names(lambdaEstimates) <- ''
    }
    frequenciesEstimatesWithConfidenceInterval <- cbind(frequenciesEstimates,perc[2:(numberOfDetectedHaplotypes+1),])
    maximumLikelihoodEstimates <- list(lambdaEstimates, frequenciesEstimatesWithConfidenceInterval, mixRadixTableOfDetectedHaplotypes, sampleSizeWithNoMissingData)
    names(maximumLikelihoodEstimates) <- c('lambda', 'haplotypes_frequencies', 'detected_haplotypes', 'used_sample_size')
  }else{
    maximumLikelihoodEstimates <- list(lambdaEstimates, frequenciesEstimates, mixRadixTableOfDetectedHaplotypes, sampleSizeWithNoMissingData)
    names(maximumLikelihoodEstimates) <- c('lambda', 'haplotypes_frequencies', 'detected_haplotypes', 'used_sample_size')
  }
  maximumLikelihoodEstimates
}

calculateMaximuLikelihoodEstimatesWithBiasCorrection <- function(dataset, numberOfAllelesAtEachMarker, isBiasCorrection=FALSE, methodForBiasCorrection='bootstrap', numberOfBootstrapReplicatesBiasCorrection=10000, pluginValueOfLambda=NULL){
  maximumLikelihoodEstimates <- calculateMaximumLikelihoodEstimatesWithOrWithoutPlugin(dataset, numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda)
  detectedHaplotypes <- as.integer(rownames(maximumLikelihoodEstimates[[2]])) - 1
  numberOfHaplotypes <- length(maximumLikelihoodEstimates[[2]])
  detectedObservations <<- dataset[[1]]
  numberOfEachDetectedObservations  <<- dataset[[2]]
  isFrequenciesEstimatesCloseToZero <- round(maximumLikelihoodEstimates[[2]],2)==0
  if(isBiasCorrection){
    sampleSize <<- sum(numberOfEachDetectedObservations)
    if(methodForBiasCorrection == "bootstrap"){
      probabilityOfEachObservation <<- numberOfEachDetectedObservations/sampleSize
      arrayOfBootstrappedEstimates <- array(0, dim = c((numberOfHaplotypes+1), numberOfBootstrapReplicatesBiasCorrection))
      rownames(arrayOfBootstrappedEstimates) <- c('lambda',(detectedHaplotypes+1))
      observation <<- seq_along(numberOfEachDetectedObservations)
      for (bootstrapReplicate in seq(numberOfBootstrapReplicatesBiasCorrection)){bootstrappedListOfDatasets <- buildBootstrappedDataset()
        bootstrappedEstimates <- calculateMaximumLikelihoodEstimatesWithOrWithoutPlugin(bootstrappedListOfDatasets, numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda)
        detectedHaplotypes   <- as.integer(rownames(bootstrappedEstimates[[2]]))
        arrayOfBootstrappedEstimates[1,bootstrapReplicate] <- unlist(bootstrappedEstimates[[1]])
        arrayOfBootstrappedEstimates[as.character(detectedHaplotypes),bootstrapReplicate] <- unlist(bootstrappedEstimates[[2]])
      }
      meanValueOfEstimates <- rowSums(arrayOfBootstrappedEstimates)/numberOfBootstrapReplicatesBiasCorrection
      meanValueOfEstimates[-1][isFrequenciesEstimatesCloseToZero] <- 0
      biasCorrectedEstimateOfLambda      <- 2*maximumLikelihoodEstimates[[1]][1] - meanValueOfEstimates[1]
      biasCorrectedEstimateOfFrequencies <- 2*maximumLikelihoodEstimates[[2]] - meanValueOfEstimates[-1]
    }else{
      if(methodForBiasCorrection=="jackknife"){
        numberOfDistinctObservations <- length(numberOfEachDetectedObservations)
        arrayOfJackknifedEstimates <- array(0, dim = c((numberOfHaplotypes+1), numberOfDistinctObservations))
        rownames(arrayOfJackknifedEstimates) <- c('lambda',(detectedHaplotypes+1))
        for(observation in 1:numberOfDistinctObservations){
          jackknifedEstimates <- calculateJackknifedEstimates(observation,numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda)# calculateMaximumLikelihoodEstimatesWithOrWithoutPlugin(bootstrappedListOfDatasets, numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda)
          rnames <- as.integer(rownames(jackknifedEstimates[[2]]))
          arrayOfJackknifedEstimates[1,observation]  <- unlist(jackknifedEstimates[[1]])
          arrayOfJackknifedEstimates[as.character(rnames),observation] <- unlist(jackknifedEstimates[[2]])
        }
        bias  <- arrayOfJackknifedEstimates %*% numberOfEachDetectedObservations/sampleSize
        bias[-1][isFrequenciesEstimatesCloseToZero] <- 0
        biasCorrectedEstimateOfLambda       <- sampleSize*maximumLikelihoodEstimates[[1]][1] - (sampleSize-1)*bias[1]# maximumLikelihoodEstimates[[1]][1] - (sampleSize-1)*( bias[1] - maximumLikelihoodEstimates[[1]][1])
        biasCorrectedEstimateOfFrequencies  <- sampleSize*maximumLikelihoodEstimates[[2]] - (sampleSize-1)*bias[-1]#maximumLikelihoodEstimates[[2]] - (sampleSize-1)*( bias[-1] - maximumLikelihoodEstimates[[2]])
      }else{
        warning("method needs to be either bootstrap or jackknife")
      }
    }
    biasCorrectedEstimates <- list(biasCorrectedEstimateOfLambda, biasCorrectedEstimateOfFrequencies)
  }else{
    biasCorrectedEstimates <- maximumLikelihoodEstimates
  }
  biasCorrectedEstimates
}

calculateMaximumLikelihoodEstimatesWithOrWithoutPlugin <- function(dataset, numberOfAllelesAtEachMarker, pluginValueOfLambda=NULL){
  if(is.null(pluginValueOfLambda)){
    out <- calculateMaximumLikelihoodEstimates(dataset, numberOfAllelesAtEachMarker)               # calculates the uncorrected estimate
  }else{
    out <- calculateMaximumLikelihoodEstimatesWithPluginValueOfLambda(dataset, numberOfAllelesAtEachMarker, pluginValueOfLambda) # calculates the corrected estimate
  }
  out
}

calculateMaximumLikelihoodEstimatesWithPluginValueOfLambda <- function(dataset,numberOfAllelesAtEachMarker,lambda){
  tolerance <- 10^-8 # Error tolerance
  detectedObservations             <- dataset[[1]]
  numberOfEachDetectedObservations <- dataset[[2]]
  #numberOfLoci  <- ncol(detectedObservations)
  subsetsFromDataset <- buildAllSetsAndSubsets(detectedObservations, numberOfAllelesAtEachMarker)
  numberOfDetectedObservations <- nrow(detectedObservations)
  Ax <- subsetsFromDataset[[1]]
  hapll <- subsetsFromDataset[[2]]
  # x  <- detectedObservations
  #
  # # This calculates the list to pick proper haplotype freuencies, i.e. sets Ay
  # hapll <- list()
  # if(length(numberOfAllelesAtEachMarker)==1){
  #   numberOfAllelesAtEachMarker <- rep(numberOfAllelesAtEachMarker,numberOfLoci)
  # }else{
  #   numberOfLoci <- length(numberOfAllelesAtEachMarker)
  # }
  # ggead <- c(1,cumprod(numberOfAllelesAtEachMarker[1:(numberOfLoci-1)]))
  # Hx <- list()
  # Ax <- list()
  # alnum <- list()
  # bin2num <- list()
  # for(k in 1:numberOfLoci){
  #   alnum[[k]] <- 1:numberOfAllelesAtEachMarker[[k]]
  #   bin2num[[k]] <- 2^(0:(numberOfAllelesAtEachMarker[[k]]-1))
  # }
  # alcnt <- array(0,numberOfLoci)
  # for(u in 1:numberOfDetectedObservations){
  #   Hx[[u]] <- list(array(0,numberOfLoci),list(),list(),list(),list(),array(0,numberOfLoci))
  #   Ax[[u]] <- list(list(),list(),list(),list())
  #   for(k in 1:numberOfLoci){
  #     temp  <- gead(x[u,k],2,numberOfAllelesAtEachMarker[[k]])
  #     temp1 <- varsets1(temp)
  #     Hx[[u]][[6]][k] <- nrow(temp1)-1
  #     Hx[[u]][[1]][k]   <- sum(temp)
  #     temp1 <- temp1[-1,]
  #     Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
  #     Hx[[u]][[3]][[k]] <- temp1
  #     Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] # subsets of alleles
  #     Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,numberOfAllelesAtEachMarker[k]) # number of alleles in subset at locus k
  #   }
  #   vz1 <- sum(Hx[[u]][[1]])
  #   temp2 <- prod(Hx[[u]][[6]])
  #   Ax[[u]][[1]] <- temp2
  #   Ax[[u]][[2]] <- varsets2(Hx[[u]][[6]])
  #   for(k in 1:numberOfLoci){
  #     Ax[[u]][[2]][,k] <- Hx[[u]][[4]][[k]][Ax[[u]][[2]][,k]]
  #   }
  #   for(j in 1:temp2){
  #     Ax[[u]][[3]][[j]] <- list()
  #     for(k in 1:numberOfLoci){
  #       temp <- gead(Ax[[u]][[2]][j,k],2,numberOfAllelesAtEachMarker[[k]])
  #       temp1 <- (alnum[[k]])[temp*alnum[[k]]]
  #       alcnt[k] <- length(temp1)
  #       Ax[[u]][[3]][[j]][[k]] <- temp1
  #     }
  #     Ax[[u]][[4]][[j]] <- varsets2(alcnt)
  #     for(k in 1:numberOfLoci){
  #       Ax[[u]][[4]][[j]][,k] <- Ax[[u]][[3]][[j]][[k]][Ax[[u]][[4]][[j]][,k]]
  #     }
  #     Ax[[u]][[4]][[j]] <- as.character((Ax[[u]][[4]][[j]]-1)%*%ggead+1)
  #     Ax[[u]][[3]][[j]] <- (-1)^(vz1+sum(alcnt))
  #   }
  #   hapll[[u]] <- Ax[[u]][[4]][[temp2]]
  # }
  hapl1 <- unique(unlist(hapll))
  # initialize parameters
  H <- length(hapl1)
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1

  #initial list#
  num0 <- pp*0
  cond1 <- 1  ## condition to stop EM alg!
  la <- lambda
  num <- num0
  rownames(num) <- hapl1
  Bcoeff <- num0
  rownames(Bcoeff) <- hapl1

  while(cond1>tolerance){
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    for(u in 1:numberOfDetectedObservations){
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[1]]){
        p <- sum(pp[Ax[[u]][[4]][[k]],])
        vz <- Ax[[u]][[3]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz
        num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],]+ exlap
        CC <- CC + exlap*p
      }
      num <- num*pp
      denom <- numberOfEachDetectedObservations[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      Bcoeff <- Bcoeff + num*denom
    }
    ppn <- Bcoeff/(sum(Bcoeff))
    cond1 <- sqrt(sum((pp-ppn)^2))
    pp <- ppn
  }
  names(la) <- NULL
  out <- list(la,pp)
  names(out) <- c('lambda', 'p')
  out
}

calculateMaximumLikelihoodEstimates <- function(dataset,numberOfAllelesAtEachMarker){
  tolerance <- 10^-8 # Error tolerance
  detectedObservations             <- dataset[[1]]
  numberOfEachDetectedObservations <- dataset[[2]]
  sampleSize <- sum(numberOfEachDetectedObservations)
  subsetsFromDataset <- buildAllSetsAndSubsets(detectedObservations, numberOfAllelesAtEachMarker)
  numberOfDetectedObservations <- nrow(detectedObservations)
  Ax <- subsetsFromDataset[[1]]
  hapll <- subsetsFromDataset[[2]]
  hapl1 <- unique(unlist(hapll))

  # initialize parameters
  H <- length(hapl1)
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1

  #initial list
  num0  <- pp*0
  cond1 <- 1  ## condition to stop EM alg!
  la    <- 2
  num   <- num0
  rownames(num)    <- hapl1
  Bcoeff           <- num0
  rownames(Bcoeff) <- hapl1
  t <- 0
  globalNumberOfIterations <- 500
  numberOfIterationsLambda <- 300
  while(cond1>tolerance && t<globalNumberOfIterations){
    t <- t+1
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    #num <- num0  #reset numerator to 0 in next iteration
    for(u in 1:numberOfDetectedObservations){
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[1]]){
        p <- sum(pp[Ax[[u]][[4]][[k]],])
        vz <- Ax[[u]][[3]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz
        num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],]+ exlap # Sum with indicator function
        CC <- CC + exlap*p
      }
      num <- num*pp
      denom <- numberOfEachDetectedObservations[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      Bcoeff <- Bcoeff + num*denom
    }
    Ccoeff <- Ccoeff/sampleSize

    # Replacing NaN's in Ak by 0
    cnt <- sum(is.nan(Bcoeff))
    if(cnt > 0){
      break
    }else{
      ppn <- Bcoeff/(sum(Bcoeff))
    }

    ### Newton step
    cond2 <- 1
    xt    <- Ccoeff   ### good initial condition
    tau   <- 0
    while(cond2 > tolerance && tau < numberOfIterationsLambda){
      tau <- tau + 1
      ex  <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      if(is.nan(xtn) || (tau == (numberOfIterationsLambda-1)) || xtn < 0){
        xtn <- runif(1, 0.1, 2.5)
      }
      cond2 <- abs(xtn-xt)
      xt <- xtn
    }
    cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
    la <- xt
    pp <- ppn
  }
  names(la) <- NULL
  out <- list(la,pp)
  names(out) <- c('lambda', 'p')
  out
}

calculateJackknifedEstimates <- function(observation, numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda){
  numberOfEachDetectedObservations[observation]   <- numberOfEachDetectedObservations[observation]-1
  observationsStillPresent   <- numberOfEachDetectedObservations != 0
  bootstrappedListOfDatasets  <- list(detectedObservations[observationsStillPresent,], numberOfEachDetectedObservations[observationsStillPresent])
  calculateMaximumLikelihoodEstimatesWithOrWithoutPlugin(bootstrappedListOfDatasets, numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda)
}

buildBootstrappedDataset <- function(){
  bootstrappedSamples <- rmultinom(1, sampleSize, probabilityOfEachObservation)
  bootstrappedDataset <- detectedObservations[rep(observation,bootstrappedSamples),]
  bootstrappedListOfDatasets  <- convertToListOfDatasets(bootstrappedDataset, idExists=FALSE)
  bootstrappedListOfDatasets
  #calculateMaximumLikelihoodEstimatesWithOrWithoutPlugin(bootstrappedListOfDatasets, numberOfAllelesAtEachMarker, pluginValueOfLambda=pluginValueOfLambda)
}

buildAllSetsAndSubsets <- function(detectedObservations, numberOfAllelesAtEachMarker){
  numberOfDetectedObservations <- nrow(detectedObservations)
  n  <- ncol(detectedObservations)
  x  <- detectedObservations
  #_____________________________
  # this calculates the list to pick proper haplotype freuencies, i.e. sets Ay
  hapll <- list()
  if(length(numberOfAllelesAtEachMarker)==1){
    numberOfAllelesAtEachMarker <- rep(numberOfAllelesAtEachMarker,n)
  }else{
    n <- length(numberOfAllelesAtEachMarker)
  }
  ggead <- c(1,cumprod(numberOfAllelesAtEachMarker[1:(n-1)]))
  Hx <- list()
  Ax <- list()
  alnum <- list()
  bin2num <- list()
  for(k in 1:n){
    alnum[[k]] <- 1:numberOfAllelesAtEachMarker[[k]]
    bin2num[[k]] <- 2^(0:(numberOfAllelesAtEachMarker[[k]]-1))
  }
  alcnt <- array(0,n)
  for(u in 1:numberOfDetectedObservations){
    Hx[[u]] <- list(array(0,n),list(),list(),list(),list(),array(0,n))
    Ax[[u]] <- list(list(),list(),list(),list())
    for(k in 1:n){
      temp  <- gead(x[u,k],2,numberOfAllelesAtEachMarker[[k]])
      temp1 <- varsets1(temp) #varsets1(temp+1)[-1,]
      Hx[[u]][[6]][k] <- nrow(temp1)-1
      Hx[[u]][[1]][k]   <- sum(temp)
      temp1 <-temp1[-1,]
      Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
      Hx[[u]][[3]][[k]] <- temp1
      Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] # subsets of alleles
      Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,numberOfAllelesAtEachMarker[k]) # number of alleles in subset at locus k
      #Hx[[u]][[6]][k]   <- length(temp1%*%bin2num[[k]])
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
        temp <- gead(Ax[[u]][[2]][j,k],2,numberOfAllelesAtEachMarker[[k]])
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
  list(Ax,hapll)
}

convertToListOfDatasets <- function(dataset, idExists = TRUE){
  # Drop id column
  if(idExists){
    dataset <- dataset[,-1]
  }
  data.comp <- apply(dataset,1,function(x) paste(x,collapse="-"))
  numberOfEachDetectedObservations <- table(data.comp)
  numberOfEachDetectedObservations.names <- names(numberOfEachDetectedObservations)
  detectedObservations <- t(sapply(numberOfEachDetectedObservations.names, function(x) unlist(strsplit(x,"-")) ))
  rownames(detectedObservations) <- NULL
  detectedObservations <- array(as.numeric(detectedObservations),dim(detectedObservations))

  # Removing samples with missing information
  sel <- rowSums(detectedObservations==0)==0
  numberOfEachDetectedObservations <- numberOfEachDetectedObservations[sel]
  detectedObservations <- detectedObservations[sel,]

  list(detectedObservations, numberOfEachDetectedObservations)
}

buildAllPossibleHaplotypes <- function(numberOfAllelesAtEachMarker){
  numberOfLoci <- length(numberOfAllelesAtEachMarker)
  H <- array(0, c(prod(numberOfAllelesAtEachMarker),numberOfLoci))
  for(i in 1:(numberOfLoci-1)){
    H[,i] <- rep(0:(numberOfAllelesAtEachMarker[i]-1), each = prod(numberOfAllelesAtEachMarker[(i+1):numberOfLoci]))
  }
  H[,numberOfLoci] <- rep(0:(numberOfAllelesAtEachMarker[numberOfLoci]-1), each = 1)
  H
}

buildAllPossibleObservations <- function(numberOfAllelesAtEachMarker){
  numberOfAllelesAtEachMarker <- 2^numberOfAllelesAtEachMarker - 1

  numberOfLoci <- length(numberOfAllelesAtEachMarker)
  H <- array(0, c(prod(numberOfAllelesAtEachMarker),numberOfLoci))
  for(i in 1:(numberOfLoci-1)){
    H[,i] <- rep(0:(numberOfAllelesAtEachMarker[i]-1), each = prod(numberOfAllelesAtEachMarker[(i+1):numberOfLoci]))
  }
  H[,numberOfLoci] <- rep(0:(numberOfAllelesAtEachMarker[numberOfLoci]-1), each = 1)
  H + 1
}

sampleWithConditionalPoisson <- function(lambda,sampleSize){
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- NULL
  x <- runif(sampleSize,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1)
  for (i in 1:sampleSize){
    k <- 1
    while(x[i] > pvec[k]){
      k <- k+1
    }
    if(k >= m){ # if a m>=100 is drawn this is executed
      k <- k+1
      a <- dpois(k,lambda)*nc
      b <- pvec[m]+a
      while(x[i]>b){
        k <- k+1
        a <- a*lambda/k
        b <- b+a
      }
    }
    out <- c(out, k)
  }
  out
}

convertDatasetToStandardFormat <- function(dataset, markersOfInterest){
  ### dataset... is the input dataset in standard format of package MLMOI, 1 st column contains smaple IDs
  ### markers ... vector of columms containing markers to be included

  allele.list <- sapply(as.list(dataset[,markersOfInterest]), function(x){
    y <- sort(unique(x))
    y[!is.na(y)]
  }
  )

  ###number of alleles per marker########
  allele.num <- unlist(lapply(allele.list,length))
  #### split data b sample ID
  dataset.split <- split(dataset[,markersOfInterest],dataset[,1])

  #### Binary representation of allele being absent and present
  samples.coded <- t(sapply(dataset.split, function(x){
    mapply(function(x,y,z){
      as.integer(is.element(y,x)) %*% 2^(0:(z-1))
    },
           x,
           allele.list,
           allele.num
    )
  }
  )
  )
  list(samples.coded,allele.list,allele.num)
}

calculatePairwiseLDWithAddons <- function(dataset, markersPair, idExists = TRUE, isConfidenceInterval=FALSE,numberOfBootstrapReplicatesConfidenceInterval=10000, significanceLevel=0.05){
  listOfDatasets  <<- convertToListOfDatasets(dataset[[1]], idExists=idExists)
  numberOfAllelesAtEachMarker <<- dataset[[3]][markersPair]
  numberOfEachDetectedObservations <- listOfDatasets[[2]]
  listOfLDEstimates <- calculatePairwiseLD(listOfDatasets,dataset)

  # Bootstrap CIs
  if(isConfidenceInterval){
    sampleSize <<- sum(numberOfEachDetectedObservations)
    probabilityOfEachObservation <<- numberOfEachDetectedObservations/sampleSize
    arrayOfBootstrappedEstimates <- array(0, dim = c(2, numberOfBootstrapReplicatesConfidenceInterval))
    observation <<- seq_along(numberOfEachDetectedObservations)
    for (bootstrapReplicate in 1:numberOfBootstrapReplicatesConfidenceInterval){
      bootstrappedListOfDatasets <- buildBootstrappedDataset()
      bootstrappedListOfLDEstimates <- calculatePairwiseLD(bootstrappedListOfDatasets,dataset)
      arrayOfBootstrappedEstimates[,bootstrapReplicate]  <- unlist(bootstrappedListOfLDEstimates)
    }
    arrayOfBootstrappedEstimates <- arrayOfBootstrappedEstimates[ , colSums(is.na(arrayOfBootstrappedEstimates))==0]
    percentiles <- t(apply(arrayOfBootstrappedEstimates, 1, quantile, c(significanceLevel/2, (1-significanceLevel/2))))
    out <- cbind(unlist(listOfLDEstimates),percentiles)
    rownames(out) <- c("D'", bquote(r^2))
    colnames(out) <- c('', paste0(as.character((significanceLevel/2)*100), '%'), paste0(as.character((1-significanceLevel/2)*100), '%'))
  }else{
    out <- unlist(listOfLDEstimates)
    names(out) <- c("D'", expression(r^2))
  }
  out
}

calculatePairwiseLD <- function(listOfDatasets,dataset){
  ListOfMaximumLikelihoodEstimates <- calculateMaximumLikelihoodEstimates(listOfDatasets,numberOfAllelesAtEachMarker)
  frequencyEstimates <- ListOfMaximumLikelihoodEstimates$p
  frequencyEstimates <- labelFrequencyEstimates(dataset,frequencyEstimates,markersPair)
  frequencyEstimates <- as.data.frame(frequencyEstimates)
  listOfLDEstimates <- list(dPrime(frequencyEstimates), rSquare(frequencyEstimates))
  names(listOfLDEstimates) <- c("D'", "r^2")
  listOfLDEstimates
}

dPrime <- function(dataset){
  colnames(dataset) <- c("A","B","freq")
  dataset$freq <- round(dataset$freq, digits=32)
  dataset <- dataset[dataset$freq >10^-100, ]

  ftbl <- xtabs(freq ~ A + B, data=dataset)

  # Haplotypes frequency
  xij <- ftbl/sum(ftbl)

  # Alleles frequencies at locus A and B, respectively
  pi <- rowSums(xij)
  qj <- colSums(xij)

  # Dij
  piqj <- pi %*% t(qj)
  Dij <- xij - piqj

  ## Dmax
  if(length(pi)==0 || length(qj)==0){
    Dprime <- 0
  }else{
    Dmax <- array(do.call(pmin, list(piqj,(1-pi) %*% t(1-qj))),dim(Dij))*(Dij<0) + array(do.call(pmin, list(pi %*% t(1-qj),(1-pi) %*% t(qj))),dim(Dij))*(Dij>=0)
    Dijprime <- Dij/Dmax
    Dprime <- round(sum(abs(Dijprime) * piqj),2)
  }
  Dprime
}

rSquare <- function(dataset){ ### Calculates R^2 according to Hedrick 1987
  colnames(dataset) <- c("A","B","freq")
  dataset$freq <- round(dataset$freq, digits=32)
  dataset <- dataset[dataset$freq >10^-100, ]

  ftbl <- xtabs(freq ~ A + B, data=dataset)  # this creates tbale of haplotype frequencies

  # Haplotypes frequency
  xij <- ftbl/sum(ftbl)

  # Alleles frequencies at locus A and B, respectively
  pi <- rowSums(xij)
  qj <- colSums(xij)

  # Dij
  piqj <- pi %*% t(qj)
  Dij <- xij - piqj

  rsquare <- round(sum(Dij^2)/((1 - sum(pi^2))*(1 - sum(qj^2) )), 2)
  rsquare
}

labelFrequencyEstimates <- function(dataset,frequencyEstimates, markersPair){ ## allist list with alleles per locus
  numberOfLoci <- length(markersPair)
  mixRadixBaseCumulativeProduct <- c(1, cumprod(dataset[[3]][markersPair])[1:(numberOfLoci-1)])
  rnames1 <- as.integer(rownames(frequencyEstimates)) - 1
  rnames  <- rnames1
  numberOfHaplotypes <- length(rnames)
  mixRadixTableOfDetectedHaplotypes <- array(0,c(numberOfHaplotypes,numberOfLoci))
  for(k in 1:numberOfLoci){ #for each locus
    re <- rnames%%rev(mixRadixBaseCumulativeProduct)[k]
    mixRadixTableOfDetectedHaplotypes[,(numberOfLoci-k+1)] <- (rnames-re)/rev(mixRadixBaseCumulativeProduct)[k]
    rnames <- re
  }
  mixRadixTableOfDetectedHaplotypes <- mixRadixTableOfDetectedHaplotypes+1
  frequencyEstimates <- cbind(mixRadixTableOfDetectedHaplotypes, frequencyEstimates)
  frequencyEstimates
}

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

varsets1 <- function(l){   
  # calculate all haplotypes that could contribute to observation (Ax)
  ## n ... number of alleles
  ## l ... 0-1 vector indicating absence/presence of alleles

  n <- length(l)
  B <- array(0,c(2^sum(l),n))
  B[1:(l[1]+1),1] <- 0:(l[1])
  lkmo <- l[1]+1
  if(n>1){
    for(k in 2:n){
      if(l[k]>0){
        lk <- lkmo*(l[k]+1)
        pick1 <- (lkmo+1):lk
        B[pick1,] <- B[rep(1:lkmo,l[k]),]
        B[pick1,k] <- rep(1:(l[k]),each=lkmo)
        lkmo <- lk   
      }
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
    out[,k] <- (x-r)/div[k]
    x <- r
  }
  out
}
