# Title        : Maximum likelihood method
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 15.08.23
# Last modified : 10.10.24

MLE <- function(dataset, n_marker, idExists=TRUE, plugin=NULL, isCI=FALSE, isBC=FALSE, replBC=10000, replCI=10000, alpha=0.05){
  ### Dropping Missing data in dataset
  samplesWithMissingData <- rowSums(dataset == 0) > 0
  dataset <- dataset[!samplesWithMissingData,]

  ### Actual sample size
  sampleSizeWithNoMissingData <- nrow(dataset)
  listOfDatasets  <- datasetXNx(dataset, idExists=idExists)
  obs <<- listOfDatasets[[1]]
  nObsVec <- listOfDatasets[[2]]
  nLoci <- ncol(obs)

  # MLEs
  mle <- MLEBC(dataset, n_marker, isBC=isBC, replBC=replBC, plugin=plugin)
  mixRadixBaseCumulativeProduct <- c(1, cumprod(n_marker)[1:(nLoci-1)])
  frequenciesEstimates <- mle[[2]]
  lbdaEstim <- mle[[1]]
  rnames1 <- as.integer(rownames(frequenciesEstimates)) - 1
  rnames  <- rnames1
  nh <- length(rnames)
  mixRadixTableOfHap <- array(0,c(nh,nLoci))
  for(k in 1:nLoci){ #for each locus
    re <- rnames%%rev(mixRadixBaseCumulativeProduct)[k]
    mixRadixTableOfHap[,(nLoci-k+1)] <- (rnames-re)/rev(mixRadixBaseCumulativeProduct)[k]
    rnames <- re
  }
  mixRadixTableOfHap <- mixRadixTableOfHap+1
  for(i in 1:nh){
    rnames[i] <- paste(mixRadixTableOfHap[i,], collapse = '')
  }
  rownames(frequenciesEstimates) <- rnames

  # Bootstrap CIs
  if(isCI){
    nHap  <- length(frequenciesEstimates)
    N <<- sum(nObsVec)
    probaObs  <<- nObsVec/N
    arrayBootEstim <- array(0, dim = c((nHap+1), replCI=replCI))
    rownames(arrayBootEstim) <- c('lambda',(rnames1+1))
    observation <<- seq_along(nObsVec)
    for (bootRepl in 1:replCI){
      bootstrappedListOfDatasets <- bootstrapDataset()
      bootEstim <- MLEPluginChoice(bootstrappedListOfDatasets, n_marker, plugin=plugin)
      hap    <- as.integer(rownames(bootEstim[[2]]))
      arrayBootEstim[1,bootRepl]  <- unlist(bootEstim[[1]])
      arrayBootEstim[as.character(hap),bootRepl] <- unlist(bootEstim[[2]])
    }
    perc <- t(apply(arrayBootEstim, 1, quantile, c(alpha/2, (1-alpha/2))))
    if(is.null(plugin)){
      lbdaEstim <- c(unlist(mle[[1]]), perc[1,])
      names(lbdaEstim) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%'))
    }else{
      lbdaEstim <- mle[[1]]
      names(lbdaEstim) <- ''
    }
    freqEstimCI <- cbind(frequenciesEstimates,perc[2:(nHap+1),])
    colnames(freqEstimCI) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%'))
    mle <- list(lbdaEstim, freqEstimCI, mixRadixTableOfHap, sampleSizeWithNoMissingData)
    names(mle) <- c('lambda', 'haplotypes_frequencies', 'detected_haplotypes', 'used_sample_size')
  }else{
    mle <- list(lbdaEstim, frequenciesEstimates, mixRadixTableOfHap, sampleSizeWithNoMissingData)
    names(mle) <- c('lambda', 'haplotypes_frequencies', 'detected_haplotypes', 'used_sample_size')
  }
  mle
}

MLEBC <- function(datasetNaturalFormat, n_marker, isBC=FALSE, replBC=10000, plugin=NULL){
  isMissingData <- rowSums(datasetNaturalFormat == 0) == 0
  datasetNaturalFormat <<- datasetNaturalFormat[isMissingData,]
  NEff <- nrow(datasetNaturalFormat)
  dataset <- datasetXNx(datasetNaturalFormat, idExists = FALSE)
  mle <- MLEPluginChoice(dataset, n_marker, plugin=plugin)
  hap <- as.numeric(rownames(mle[[2]])) - 1
  nHaplotypes <- length(mle[[2]])
  obs <<- dataset[[1]]
  nObsVec  <<- dataset[[2]]
  isFreqZero <- round(mle[[2]],2)==0
  if(isBC){
    N <<- sum(nObsVec)
    probaObs <<- nObsVec/N
    arrayBootEstim <- array(0, dim = c((nHaplotypes+1), replBC))
    rownames(arrayBootEstim) <- c('lambda',(hap+1))
    for(bootRepl in seq(replBC)){
      bootDataset <- bootstrapDataset()
      bootEstim <- MLEPluginChoice(bootDataset, n_marker, plugin=plugin)
      hapl   <- as.numeric(rownames(bootEstim[[2]]))
      arrayBootEstim[1,bootRepl] <- unlist(bootEstim[[1]])
      tempBootFreqEstim <- unlist(bootEstim[[2]])
      arrayBootEstim[as.character(hapl),bootRepl] <- tempBootFreqEstim
    }
    meanValueOfEstimates <- rowSums(arrayBootEstim)/replBC
    meanValueOfEstimates[-1][isFreqZero] <- 0
    BCLambda      <- 2*mle[[1]][1] - meanValueOfEstimates[1]
    BCFrequencies <- 2*mle[[2]] - meanValueOfEstimates[-1]
    biasCorrectedEstimates <- list(BCLambda, BCFrequencies)
  }else{
    biasCorrectedEstimates <- mle
  }
  biasCorrectedEstimates
}

MLEPluginChoice <- function(dataset, n_marker, plugin=NULL){
  if(is.null(plugin)){
    out <- baseModel(dataset, n_marker)               # calculates the uncorrected estimate
  }else{
    out <- MLEPlugin(dataset, n_marker, plugin=plugin) # calculates the corrected estimate
  }
  out
}

MLEPlugin <- function(dataset,n_marker,plugin){
  tolerance <- 10^-8 # Error tolerance
  obs <- dataset[[1]]
  nObsVec <- dataset[[2]]
  subsetsFromDataset <- modelSubsets(obs, n_marker)
  nObs <- nrow(obs)
  Ax <- subsetsFromDataset[[1]]
  hapll <- subsetsFromDataset[[2]]
  hapl1 <- unique(unlist(hapll))
  # initialize parameters
  H <- length(hapl1)
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1

  #initial list#
  num0 <- pp*0
  cond1 <- 1  ## condition to stop EM alg!
  la <- plugin
  num <- num0
  rownames(num) <- hapl1
  Bcoeff <- num0
  rownames(Bcoeff) <- hapl1

  while(cond1>tolerance){
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    for(u in 1:nObs){
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
      denom <- nObsVec[u]/denom
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

baseModel <- function(dataset,n_marker){
  tolerance <- 10^-8 # Error tolerance
  obs     <- dataset[[1]]
  nObsVec <- dataset[[2]]
  N <- sum(nObsVec)
  subsetsFromDataset <- modelSubsets(obs, n_marker)
  nObs <- nrow(obs)
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
  nIterationsLambda <- 300
  while(cond1>tolerance && t<globalNumberOfIterations){
    t <- t+1
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    for(u in 1:nObs){
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
      denom <- nObsVec[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      Bcoeff <- Bcoeff + num*denom
    }
    Ccoeff <- Ccoeff/N

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
    while(cond2 > tolerance && tau < nIterationsLambda){
      tau <- tau + 1
      ex  <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      if(is.nan(xtn) || (tau == (nIterationsLambda-1)) || xtn < 0){
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

bootstrapDataset <- function(){
  bootSamples <- rmultinom(1, N, probaObs)
  bootDataset <- obs[rep(1:nrow(obs),bootSamples),]
  bootDataset  <- datasetXNx(bootDataset, idExists=FALSE)
  bootDataset
}

modelSubsets <- function(obs, n_marker){
  nObs <- nrow(obs)
  n  <- ncol(obs)
  x  <- obs
  #_____________________________
  # this calculates the list to pick proper haplotype freuencies, i.e. sets Ay
  hapll <- list()
  if(length(n_marker)==1){
    n_marker <- rep(n_marker,n)
  }else{
    n <- length(n_marker)
  }
  ggead <- c(1,cumprod(n_marker[1:(n-1)]))
  Hx <- list()
  Ax <- list()
  alnum <- list()
  bin2num <- list()
  for(k in 1:n){
    alnum[[k]] <- 1:n_marker[[k]]
    bin2num[[k]] <- 2^(0:(n_marker[[k]]-1))
  }
  alcnt <- array(0,n)
  for(u in 1:nObs){
    Hx[[u]] <- list(array(0,n),list(),list(),list(),list(),array(0,n))
    Ax[[u]] <- list(list(),list(),list(),list())
    for(k in 1:n){
      temp  <- gead(x[u,k],2,n_marker[[k]])
      temp1 <- varsets1(temp) #varsets1(temp+1)[-1,]
      Hx[[u]][[6]][k] <- nrow(temp1)-1
      Hx[[u]][[1]][k]   <- sum(temp)
      temp1 <-temp1[-1,]
      Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
      Hx[[u]][[3]][[k]] <- temp1
      Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] # subsets of alleles
      Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,n_marker[k]) # number of alleles in subset at locus k
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
        temp <- gead(Ax[[u]][[2]][j,k],2,n_marker[[k]])
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

datasetXNx <- function(dataset, idExists = TRUE){
  # Drop id column
  if(idExists){
    dataset <- dataset[,-1]
  }
  data.comp <- apply(dataset,1,function(x) paste(x,collapse="-"))
  nObsVec <- table(data.comp)
  nObsVec.names <- names(nObsVec)
  obs <- t(sapply(nObsVec.names, function(x) unlist(strsplit(x,"-")) ))
  rownames(obs) <- NULL
  obs <- array(as.numeric(obs),dim(obs))

  # Removing samples with missing information
  sel <- rowSums(obs==0)==0
  nObsVec <- nObsVec[sel]
  obs <- obs[sel,]

  list(obs, nObsVec)
}

allHaplotypes <- function(n_marker){
  nLoci <- length(n_marker)
  H <- array(0, c(prod(n_marker),nLoci))
  for(i in 1:(nLoci-1)){
    H[,i] <- rep(0:(n_marker[i]-1), each = prod(n_marker[(i+1):nLoci]))
  }
  H[,nLoci] <- rep(0:(n_marker[nLoci]-1), each = 1)
  H
}

allObservations <- function(n_marker){
  n_marker <- 2^n_marker - 1
  nLoci <- length(n_marker)
  H <- array(0, c(prod(n_marker),nLoci))
  for(i in 1:(nLoci-1)){
    H[,i] <- rep(0:(n_marker[i]-1), each = prod(n_marker[(i+1):nLoci]))
  }
  H[,nLoci] <- rep(0:(n_marker[nLoci]-1), each = 1)
  H + 1
}

conditionalPoissonSampling <- function(lambda,N){
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- NULL
  x <- runif(N,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1)
  for (i in 1:N){
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

datasetToStandard <- function(dataset, markersOfInterest){
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

pairwiseLD <- function(dataset, markersPair, idExists = TRUE, isCI=FALSE,replCI=10000, alpha=0.05){
  listOfDatasets  <<- datasetXNx(dataset[[1]], idExists=idExists)
  n_marker <<- dataset[[3]][markersPair]
  nObsVec <- listOfDatasets[[2]]
  listOfLDEstimates <- pairwiseLDBase(listOfDatasets,dataset)

  # Bootstrap CIs
  if(isCI){
    N <<- sum(nObsVec)
    probaObs <<- nObsVec/N
    arrayBootEstim <- array(0, dim = c(2, replCI))
    observation <<- seq_along(nObsVec)
    for (bootRepl in 1:replCI){
      bootDataset <- bootstrapDataset()
      bootEstim <- pairwiseLDBase(bootDataset,dataset)
      arrayBootEstim[,bootRepl]  <- unlist(bootEstim)
    }
    arrayBootEstim <- arrayBootEstim[ , colSums(is.na(arrayBootEstim))==0]
    perc <- t(apply(arrayBootEstim, 1, quantile, c(alpha/2, (1-alpha/2))))
    out <- cbind(unlist(listOfLDEstimates),perc)
    rownames(out) <- c("D'", bquote(r^2))
    colnames(out) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%'))
  }else{
    out <- unlist(listOfLDEstimates)
    names(out) <- c("D'", expression(r^2))
  }
  out
}

pairwiseLDBase <- function(listOfDatasets,dataset){
  ListOfMaximumLikelihoodEstimates <- baseModel(listOfDatasets,n_marker)
  freqEstim <- ListOfMaximumLikelihoodEstimates$p
  freqEstim <- labelFreqEstim(dataset,freqEstim,markersPair)
  freqEstim <- as.data.frame(freqEstim)
  listOfLDEstimates <- list(dPrime(freqEstim), rSquare(freqEstim))
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

labelFreqEstim <- function(dataset,freqEstim, markersPair){ ## allist list with alleles per locus
  nLoci <- length(markersPair)
  mixRadixBaseCumulativeProduct <- c(1, cumprod(dataset[[3]][markersPair])[1:(nLoci-1)])
  rnames1 <- as.integer(rownames(freqEstim)) - 1
  rnames  <- rnames1
  nHaplotypes <- length(rnames)
  mixRadixTableOfHap <- array(0,c(nHaplotypes,nLoci))
  for(k in 1:nLoci){ #for each locus
    re <- rnames%%rev(mixRadixBaseCumulativeProduct)[k]
    mixRadixTableOfHap[,(nLoci-k+1)] <- (rnames-re)/rev(mixRadixBaseCumulativeProduct)[k]
    rnames <- re
  }
  mixRadixTableOfHap <- mixRadixTableOfHap+1
  freqEstim <- cbind(mixRadixTableOfHap, freqEstim)
  freqEstim
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

baseModelSim <- function(dataset,n_marker){
  samplesWithMissingData <- rowSums(dataset == 0) > 0
  dataset <- dataset[!samplesWithMissingData,]

  ### Actual sample size
  sampleSizeWithNoMissingData <- nrow(dataset)
  listOfDatasets  <- datasetXNx(dataset, idExists=FALSE)
  tolerance <- 10^-8 # Error tolerance
  obs     <- listOfDatasets[[1]]
  nObsVec <- listOfDatasets[[2]]
  N <- sum(nObsVec)
  subsetsFromDataset <- modelSubsets(obs, n_marker)
  nObs <- nrow(obs)
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
  nIterationsLambda <- 300
  while(cond1>tolerance && t<globalNumberOfIterations){
    t <- t+1
    Ccoeff <- 0
    Bcoeff <- num0 #reset B coefficients to 0 in next iteration
    for(u in 1:nObs){
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
      denom <- nObsVec[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      Bcoeff <- Bcoeff + num*denom
    }
    Ccoeff <- Ccoeff/N

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
    while(cond2 > tolerance && tau < nIterationsLambda){
      tau <- tau + 1
      ex  <- exp(-xt)
      xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)
      if(is.nan(xtn) || (tau == (nIterationsLambda-1)) || xtn < 0){
        xtn <- runif(1, 0.1, 2.5)
      }
      cond2 <- abs(xtn-xt)
      xt <- xtn
    }
    cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
    la <- xt
    pp <- ppn
  }
  ## Ordering the frequencies
  pp <- pp[order(as.numeric(rownames(pp))), ]

  ## Setting the frequencies of the unobserved haplotypes to 0.0
  nhapl <- prod(n_marker)

  if(length(pp)<nhapl){
    out <- t(pp)
    name <- colnames(out)
    cnt <- 0
    for (i in 1:nhapl) {
      if (is.element(as.character(i), name)){
        cnt <- cnt + 1
      }else{
        pp <- append(pp, list(x = 0.0), i-1)
      }
    }
  }
  names(la) <- NULL
  out <- list(la,unlist(pp))
  names(out) <- c('lambda', 'p')
  out
}

CRLB <- function(freq,lambda,N,nloci){
  #### Generate all possible observations given nloci
  detectedObservations    <- allObservations(nloci)
  Nobs <- nrow(detectedObservations)

  # Find for each infection detectedObservations all components necessary for the computations
  Ax <- modelSubsets(detectedObservations, nloci)[[1]]

  # Initialize the Fisher information matrix of degree d
  names(freq) <- seq_along(freq)
  pp   <- matrix(freq, ncol=1)
  hap  <- as.numeric(names(freq))
  rownames(pp) <- hap

  d <- nrow(pp)+2
  diag.el <- diag(array(1:d^2,c(d,d)))[-1]
  diag.el <- diag.el[-length(diag.el)]

  I <- matrix(0, nrow=d, ncol=d) #dof
  rownames(I) <- 1:d
  colnames(I) <- 1:d

  elmo <- exp(lambda)-1

  # Ipsi,psi
  out <- 0
  for(u in 1:Nobs){ # For each observation x in ScrO
    out1 <- 0
    out2 <- 0
    for(k in 1:(Ax[[u]][[1]][[1]])){  # For each observation y in the sub-observation ScrAx
      sump <- sum(pp[Ax[[u]][[4]][[k]],])
      vz   <- Ax[[u]][[3]][[k]]
      out1 <- out1 + vz*GFunc(lambda,sump)   # GFunc to build Px
      out2 <- out2 + vz*dGFunc(lambda,sump)  # dG/dl to build dPx/dlam
    }
    out <- out + (out2^2)/out1 #Ill
  }
  I[1,1] <- N*out/(dPsi(lambda)^2)   # I_lam_lam

  # Ipi,pj  new
  out <- 0 * (freq%*%t(freq))
  for(u in 1:Nobs){
    out1  <- 0
    out2 <- 0*freq
    for(k in 1:Ax[[u]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[Ax[[u]][[4]][[k]],])
      elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1 + vz*GFunc(lambda,sump)
      out2[as.numeric(Ax[[u]][[4]][[k]])] <- out2[as.numeric(Ax[[u]][[4]][[k]])] + vz*lambda*elp/elmo   # dG/dpi * Indic
    }
    out <- out + (out2 %*% t(out2/out1))
  }
  I[2:(d-1),2:(d-1)] <- N*out

  # Ipsi,pi
  out <- 0*freq
  for(u in 1:Nobs){
    out1  <- 0
    out21 <- 0
    out22 <- 0*freq
    for(k in 1:Ax[[u]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[Ax[[u]][[4]][[k]],])
      elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1  + vz*GFunc(lambda,sump)             # GFunc
      out21 <- out21 + vz*dGFunc(lambda,sump)            # dG/dl
      out22[as.numeric(Ax[[u]][[4]][[k]])] <- out22[as.numeric(Ax[[u]][[4]][[k]])] + vz*(lambda*elp/elmo)    # dG/dpi
    }
    out <- out + (out21*out22)/out1
  }
  tmpI <- N*out/dPsi(lambda)
  I[1,2:(d-1)] <- tmpI
  I[2:(d-1),1] <- tmpI

  # Parameter space in Higher dimension
  dbeta <- c(0,rep(1,prod(nloci)),0)
  I[,d] <- dbeta
  I[d,] <- dbeta

  print("Inverting Now...")
  out <- matrix.inverse(I)
  out <- out[-d,]
  out <- out[,-d]
  I <- I[-d,]
  I <- I[,-d]
  list(I,zapsmall(out))
}

GFunc <- function(lambda, sumFreq){
  elp   <- exp(lambda*sumFreq)
  elmo  <- exp(lambda) - 1
  (elp - 1)/elmo
}

dGFunc <- function(lambda, freq){
  el <- exp(lambda)
  elmo <- exp(lambda)-1
  elp  <- exp(lambda*freq)

  freq*elp/elmo - el*(elp - 1)/(elmo^2)
}

dPsi <- function(lambda){
  eml <- 1-exp(-lambda)
  1/eml - (lambda*exp(-lambda))/(eml^2)
}
