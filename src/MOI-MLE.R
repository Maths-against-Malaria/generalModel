# Title        : Maximum likelihood method
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 15.08.23
# Last modified : 06.01.25

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

MLE <- function(data, markers, plugin=NULL, isCI=FALSE, isBC=FALSE, replBC=10000, replCI=10000, alpha=0.05, allelesName=TRUE){
  dat <- datasetFormat(data, 2:ncol(data))
  dataset <- dat[[1]][,markers]
  alleleList <- dat[[2]][markers]
  GA <- dat[[3]][markers]

  ### Dropping Missing data in dataset
  missData <- rowSums(dataset == 0) > 0
  dataset <- dataset[!missData,]

  ### Actual sample size
  Neff <- nrow(dataset)
  XNx  <- datasetXNx(dataset)
  X <- XNx[[1]]
  Nx <- XNx[[2]]
  nLoci <- ncol(X)

  # MLEs
  mle <- MLEBC(dataset, GA, isBC=isBC, replBC=replBC, plugin=plugin)
  mixRadCumProd <- c(1, cumprod(GA)[1:(nLoci-1)])
  freqEstim <- mle[[2]]
  lbdaEstim <- mle[[1]]
  rnames1 <- as.integer(rownames(freqEstim)) - 1
  rnames  <- rnames1
  nh <- length(rnames)
  mixRadixTableOfHap <- array(0,c(nh,nLoci))
  for(k in 1:nLoci){ #for each locus
    re <- rnames%%rev(mixRadCumProd)[k]
    mixRadixTableOfHap[,(nLoci-k+1)] <- as.character((rnames-re)/rev(mixRadCumProd)[k] + 1)
    rnames <- re
  }
  
  if(allelesName){
    freqEstim <- hapName(freqEstim,alleleList)
    for (i in 1:length(alleleList)){
      mixRadixTableOfHap[,i] <- alleleList[[i]][as.numeric(mixRadixTableOfHap[,i])]
    }
  }else {
    for(i in 1:nh){
      rnames[i] <- paste(mixRadixTableOfHap[i,], collapse = '.')
    }
    rownames(freqEstim) <- rnames
    mixRadixTableOfHap <- apply(mixRadixTableOfHap, 2, as.numeric)
  }
  
  colnames(mixRadixTableOfHap) <- names(GA)

  # Bootstrap CIs
  if(isCI){
    nHap  <- length(freqEstim)
    N <- sum(Nx)
    probaObs  <- Nx/N
    arrayBootEstim <- array(0, dim = c((nHap+1), replCI=replCI))
    rownames(arrayBootEstim) <- c('lambda',(rnames1+1))
    observation <- seq_along(Nx)
    for (bootRepl in 1:replCI){
      bootstrappedListOfDatasets <- bootstrapDataset(X, N, probaObs)
      bootEstim <- MLEPluginChoice(bootstrappedListOfDatasets, GA, plugin=plugin)
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
    freqEstimCI <- cbind(freqEstim,perc[2:(nHap+1),])
    colnames(freqEstimCI) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%'))
    mle <- list(lbdaEstim, freqEstimCI, mixRadixTableOfHap, Neff)
    names(mle) <- c('lambda', 'haplotypes_frequencies', 'detected_haplotypes', 'used_sample_size')
  }else{
    mle <- list(lbdaEstim, freqEstim, mixRadixTableOfHap, Neff)
    names(mle) <- c('lambda', 'haplotypes_frequencies', 'detected_haplotypes', 'used_sample_size')
  }
  rownames(mle[[2]]) <- paste("p", rownames(mle[[2]]), sep="") 
  mle
}

MLEBC <- function(data, GA, isBC=FALSE, replBC=10000, plugin=NULL){
  isMissingData <- rowSums(data == 0) == 0
  data <- data[isMissingData,]
  NEff <- nrow(data)
  dataset <- datasetXNx(data)
  mle <- MLEPluginChoice(dataset, GA, plugin=plugin)
  hap <- as.numeric(rownames(mle[[2]])) - 1
  nHaplotypes <- length(mle[[2]])
  X <- dataset[[1]]
  Nx  <- dataset[[2]]
  isFreqZero <- round(mle[[2]],2)==0
  if(isBC){
    N <- sum(Nx)
    probaObs <- Nx/N
    arrayBootEstim <- array(0, dim = c((nHaplotypes+1), replBC))
    rownames(arrayBootEstim) <- c('lambda',(hap+1))
    for(bootRepl in seq(replBC)){
      bootDataset <- bootstrapDataset(X, N, probaObs)
      bootEstim <- MLEPluginChoice(bootDataset, GA, plugin=plugin)
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

MLEPluginChoice <- function(data, GA, plugin=NULL){
  if(is.null(plugin)){
    out <- baseModel(data, GA)               # calculates the uncorrected estimate
  }else{
    out <- MLEPlugin(data, GA, plugin=plugin) # calculates the corrected estimate
  }
  out
}

MLEPlugin <- function(data, GA, plugin){
  tol <- 10^-8 # Error tol
  X <- data[[1]]
  Nx <- data[[2]]
  allSubs <- modelSubs(X, GA)
  nObs <- nrow(X)
  Ax <- allSubs[[1]]
  hapll <- allSubs[[2]]
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

  while(cond1>tol){
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
      denom <- Nx[u]/denom
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

baseModel <- function(data, GA){
  eps <- 10^-8
  X <- data[[1]]
  Nx <- data[[2]]
  N <- sum(Nx)
  allSubs <- modelSubs(X, GA)
  nn <- nrow(X)
  Ax <- allSubs[[1]]
  hapll <- allSubs[[2]]
  hapl1 <- unique(unlist(hapll))
  
  #_____________________________
  # this calculates the list to pick proper haplotype freuencies, i.e. sets Ay
  
  hapll <- list()
  H <- length(hapl1)
  attemp <- 0
  pp <- array(rep(1/H,H),c(H,1))
  rownames(pp) <- hapl1
  
  #initial list#
  num0 <- pp*0
  cond1 <- 1  ## condition to stop EM alg! 
  la <- 2
  num <- num0
  rownames(num) <- hapl1
  Bcoeff <- num0
  rownames(Bcoeff) <- hapl1
  
  ## If there is no sign of super-infections the estimate is degenerate p =Nx/N, and it is treated separately in the if-else statement
  if(any(lapply(Ax,function(x) length(x[[3]]))>1)){ # generic case
    stp <- 0 # this is to try different initial conditions if convergence is slow, if stp reaches 120 step a new initial condition is used,
    # at most 50 different initial conditions are used
    
    while(cond1>eps){
      stp <- stp + 1 
      Ccoeff <- 0
      Bcoeff <- num0 #reset B coefficients to 0 in next iteration
      num <- num0  #reset numerator to 0 in next iteration
      for(u in 1:nn){
        denom <- 0
        num <- num0
        CC <- 0
        for(k in 1:Ax[[u]][[1]]){
          p <- sum(pp[Ax[[u]][[4]][[k]],])
          vz <- Ax[[u]][[3]][[k]]
          lap <- la*p
          exlap <- vz*exp(lap)
          denom <- denom + exlap-vz  ##   = (1-)^(Nx-Ny)*(Exp(lambda*sum p)-1) = (1-)^(Nx-Ny)*G(sum p)
          num[Ax[[u]][[4]][[k]],] <- num[Ax[[u]][[4]][[k]],]+ exlap#*pp[Ax[[u]][[1]][[k]],]
          ## exlap =  (1-)^(Nx-Ny) G'(sum p)   --- denominator of generating functions cancels out!
          CC <- CC + exlap*p
        }
        num <- num*pp
        denom <- Nx[u]/denom
        denom <- la*denom
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
      cond1 <- abs(xt-la) + sqrt(sum((pp-ppn)^2))
      la <- xt
      pp <- ppn
      if(stp == 120 & attemp < 50 ){ # if algorithm dod not converge within 120 steps try new initial condition
        #print(attemp)
        attemp <- attemp + 1
        stp <- 0
        # new initial conditions
        
        tmp <- table(unlist(rep(lapply(Ax, function(x) x[[4]][[length(x[[4]])]]),Nx)))
        #pp1 <- pp
        
        pp[names(tmp),] <- tmp/sum(tmp)
        num0 <- pp*0
        cond1 <- 1  ## condition to stop EM alg! 
        la <- 2 + attemp
        num <- num0
        rownames(num) <- hapl1
        Bcoeff <- num0
        rownames(Bcoeff) <- hapl1
      }
    }
  }else{
    la <- 0
    pp[,1] <- Nx/N 
  }
  out <- list(unname(la),pp)
  names(out) <- c("MOI.est","Freq.est")
  out
}

bootstrapDataset <- function(X, N, probaObs){
  bootSamples <- rmultinom(1, N, probaObs)
  bootDataset <- X[rep(1:nrow(X),bootSamples),]
  bootDataset  <- datasetXNx(bootDataset)
  bootDataset
}

hapName <- function(pp, allist){ ## allist list with alleles per locos
  allist <- lapply(allist,as.character)  # this makes sure factor-levels are used
  n <- length(allist)
  allnum <- unlist(lapply(allist,length))
  hapl <- geadrepr(as.numeric(rownames(pp))-1,allnum)+1
  hapl1 <- array(NA,dim(hapl))
  for(l in 1: ncol(hapl)){
    for(m in 1: nrow(hapl)){
      hapl1[m,l] <- allist[[l]][hapl[m,l]]
    }
  }
  hapl1 <- apply(hapl1,1,function(x) paste(x,sep="",collapse="."))
  rownames(pp) <- hapl1
  pp
}

modelSubs <- function(X, GA){
  nObs <- nrow(X)
  n  <- ncol(X)
  x  <- X
  #_____________________________
  # this calculates the list to pick proper haplotype freuencies, i.e. sets Ay
  hapll <- list()
  if(length(GA)==1){
    GA <- rep(GA,n)
  }else{
    n <- length(GA)
  }
  ggead <- c(1,cumprod(GA[1:(n-1)]))
  Hx <- list()
  Ax <- list()
  alnum <- list()
  bin2num <- list()
  for(k in 1:n){
    alnum[[k]] <- 1:GA[[k]]
    bin2num[[k]] <- 2^(0:(GA[[k]]-1))
  }
  alcnt <- array(0,n)
  for(u in 1:nObs){
    Hx[[u]] <- list(array(0,n),list(),list(),list(),list(),array(0,n))
    Ax[[u]] <- list(list(),list(),list(),list())
    for(k in 1:n){
      temp  <- gead(x[u,k],2,GA[[k]])
      temp1 <- varsets1(temp) #varsets1(temp+1)[-1,]
      Hx[[u]][[6]][k] <- nrow(temp1)-1
      Hx[[u]][[1]][k]   <- sum(temp)
      temp1 <-temp1[-1,]
      Hx[[u]][[2]][[k]] <- (alnum[[k]])[temp*alnum[[k]]]
      Hx[[u]][[3]][[k]] <- temp1
      Hx[[u]][[4]][[k]] <- temp1%*%bin2num[[k]] # subsets of alleles
      Hx[[u]][[5]][[k]] <- Hx[[u]][[3]][[k]]%*%rep(1,GA[k]) # number of alleles in subset at locus k
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
        temp <- gead(Ax[[u]][[2]][j,k],2,GA[[k]])
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

datasetXNx <- function(data){
  data.comp <- apply(data,1,function(x) paste(x,collapse="-"))
  Nx <- table(data.comp)
  Nx.names <- names(Nx)
  X <- t(sapply(Nx.names, function(x) unlist(strsplit(x,"-")) ))
  rownames(X) <- NULL
  X <- array(as.numeric(X),dim(X))

  # Removing samples with missing information
  sel <- rowSums(X==0)==0
  Nx <- Nx[sel]
  X <- X[sel,]

  list(X, Nx)
}

allHap <- function(GA){
  nLoci <- length(GA)
  H <- array(0, c(prod(GA),nLoci))
  for(i in 1:(nLoci-1)){
    H[,i] <- rep(0:(GA[i]-1), each = prod(GA[(i+1):nLoci]))
  }
  H[,nLoci] <- rep(0:(GA[nLoci]-1), each = 1)
  H
}

allObs <- function(GA){
  GA <- 2^GA - 1
  nLoci <- length(GA)
  H <- array(0, c(prod(GA),nLoci))
  for(i in 1:(nLoci-1)){
    H[,i] <- rep(0:(GA[i]-1), each = prod(GA[(i+1):nLoci]))
  }
  H[,nLoci] <- rep(0:(GA[nLoci]-1), each = 1)
  H + 1
}

cPoiss <- function(lambda,N){
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

datasetFormat <- function(data, markers){
  ### dataset... is the input data in standard format of package MLMOI, 1 st column contains smaple IDs
  ### markers ... vector of columms containing markers to be included

  allele.list <- sapply(as.list(data[,markers]), function(x){
    y <- sort(unique(x))
    y[!is.na(y)]
  }
  )

  ###number of alleles per marker########
  allele.num <- unlist(lapply(allele.list,length))

  #### split data b sample ID
  dataset.split <- split(data[,markers],data[,1])

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
  #names(allele.num) <- paste0('m',1:length(allele.num))
  out <- list(samples.coded,allele.list,allele.num)
  names(out) <- c('dataset', 'alleles', 'genetic_architecture')
  out
}

pairwiseLD <- function(data, markersPair, isCI=FALSE,replCI=10000, alpha=0.05){
  data <- datasetFormat(data, 2:ncol(data))

  XNx  <- datasetXNx(data[[1]][,markersPair])
  X <- XNx[[1]]
  GA <- data[[3]][markersPair]
  Nx <- XNx[[2]]
  ldEstim <- pairwiseLDBase(XNx,data,markersPair, GA)

  # Bootstrap CIs
  if(isCI){
    N <- sum(Nx)
    probaObs <- Nx/N
    arrayBootEstim <- array(0, dim = c(2, replCI))
    observation <- seq_along(Nx)
    for (bootRepl in 1:replCI){
      bootDataset <- bootstrapDataset(X, N, probaObs)
      bootEstim <- pairwiseLDBase(bootDataset,data,markersPair,GA)
      arrayBootEstim[,bootRepl]  <- unlist(bootEstim)
    }
    arrayBootEstim <- arrayBootEstim[ , colSums(is.na(arrayBootEstim))==0]
    perc <- t(apply(arrayBootEstim, 1, quantile, c(alpha/2, (1-alpha/2))))
    out <- cbind(unlist(ldEstim),perc)
    rownames(out) <- c("D'", bquote(r^2))
    colnames(out) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%'))
  }else{
    out <- unlist(ldEstim)
    names(out) <- c("D'", expression(r^2))
  }
  out
}

pairwiseLDBase <- function(XNx, data, markersPair, GA){
  mle <- MLEPluginChoice(XNx, GA, plugin=NULL)
  freqEstim <- mle[[2]]
  freqEstim <- labelFreqEstim(data, freqEstim, markersPair)
  freqEstim <- as.data.frame(freqEstim)
  ldEstim <- list(dPrime(freqEstim), rSquare(freqEstim))
  names(ldEstim) <- c("D'", "r^2")
  ldEstim
}

dPrime <- function(data){
  colnames(data) <- c("A","B","freq")
  data$freq <- round(data$freq, digits=32)
  data <- data[data$freq >10^-100, ]

  ftbl <- xtabs(freq ~ A + B, data=data)

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

rSquare <- function(data){ ### Calculates R^2 according to Hedrick 1987
  colnames(data) <- c("A","B","freq")
  data$freq <- round(data$freq, digits=32)
  data <- data[data$freq >10^-100, ]

  ftbl <- xtabs(freq ~ A + B, data=data)  # this creates tbale of haplotype frequencies

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

labelFreqEstim <- function(data,freqEstim, markersPair){ ## allist list with alleles per locus
  nLoci <- length(markersPair)
  mixRadCumProd <- c(1, cumprod(data[[3]][markersPair])[1:(nLoci-1)])
  rnames1 <- as.integer(rownames(freqEstim)) - 1
  rnames  <- rnames1
  nHaplotypes <- length(rnames)
  mixRadixTableOfHap <- array(0,c(nHaplotypes,nLoci))
  for(k in 1:nLoci){ #for each locus
    re <- rnames%%rev(mixRadCumProd)[k]
    mixRadixTableOfHap[,(nLoci-k+1)] <- (rnames-re)/rev(mixRadCumProd)[k]
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

geadrepr <- function(x,l){   ## calculates general geadic expression of each element of vector x
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

baseModelSim <- function(data,GA){
  dat <- rowSums(data == 0) > 0
  data <- data[!dat,]

  ### Actual sample size
  #Neff <- nrow(data)
  XNx  <- datasetXNx(data)
  tol <- 10^-8 # Error tolerance
  X     <- XNx[[1]]
  Nx <- XNx[[2]]
  N <- sum(Nx)
  allSubs <- modelSubs(X, GA)
  nObs <- nrow(X)
  Ax <- allSubs[[1]]
  hapll <- allSubs[[2]]
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
  while(cond1>tol && t<globalNumberOfIterations){
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
      denom <- Nx[u]/denom
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
    while(cond2 > tol && tau < nIterationsLambda){
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
  nhapl <- prod(GA)

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

FI <- function(data, markers, isPsi = FALSE, isPrev = FALSE, allelesName=TRUE, isObserv = FALSE){
  dat <- datasetFormat(data, 2:ncol(data))
  alleleList <- dat[[2]][markers]
  GA <- dat[[3]][markers]

  mle <- MLE(data, markers, allelesName = FALSE)

  if(isPsi){
    if(isPrev){
      out <- FIPsiPrev(data, markers, isObserv=isObserv)
      var <- diag(out)
      if(allelesName){
        alleles <- strsplit(gsub("[p]", "", rownames(mle[[2]])), "[.]")
        rownames(mle[[2]]) <- sapply(alleles, function(x) rank(as.numeric(x), GA))
        hapname <- paste("q",rownames(hapName(mle[[2]], alleleList)), sep = "")
      }else{
        hapname <- rownames(mle[[2]])
      }
      nameMean <- c('Mean_MOI', hapname)
      names(var) <- nameMean
      rownames(out) <- nameMean
      colnames(out) <- nameMean
      out <- list(out, var)
      names(out) <- c('Covariance matrix', 'Variance')
      out
    }else {
      out <- FIPsi(data, markers, isObserv=isObserv)
      var <- diag(out)
      if(allelesName){
        alleles <- strsplit(gsub("[p]", "", rownames(mle[[2]])), "[.]")
        rownames(mle[[2]]) <- sapply(alleles, function(x) rank(as.numeric(x), GA))
        hapname <- paste("p",rownames(hapName(mle[[2]], alleleList)), sep = "")
      }else{
        hapname <- rownames(mle[[2]])
      }
      nameMean <- c('Mean_MOI', hapname)
      names(var) <- nameMean
      rownames(out) <- nameMean
      colnames(out) <- nameMean
      out <- list(out, var)
      names(out) <- c('Covariance matrix', 'Variance')
      out
    }
  }else {
    if(isObserv){
      out <- round(solve(FIMObs(data, markers)),5)
    }else{
      out <- round(solve(FIM(mle, GA)),5)
    }
    out <- out[-nrow(out), -ncol(out)]
    var <- diag(out)
    if(allelesName){
      alleles <- strsplit(gsub("[p]", "", rownames(mle[[2]])), "[.]")
      rownames(mle[[2]]) <- sapply(alleles, function(x) rank(as.numeric(x), GA))
      hapname <- paste("p",rownames(hapName(mle[[2]], alleleList)), sep = "")
    }else{
      hapname <- rownames(mle[[2]])
    }
    names(var) <- c('Lambda', hapname)
    rownames(out) <- c('Lambda', hapname)
    colnames(out) <- c('Lambda', hapname)
    out <- list(out, var)
    names(out) <- c('Covariance matrix', 'Variance')
    out
  }
  out
}

FIPsiPrev <- function(data, markers, isObserv = FALSE){

  dat <- datasetFormat(data, 2:ncol(data))
  GA <- dat[[3]][markers]
  mle <- MLE(data, markers, allelesName = FALSE)
  if(isObserv){
    I <- FIMObs(data, markers)
  }else{
    I <- FIM(mle, GA)
  }
  
  lambda <- mle[[1]]
  dgdl <- c(dGdL(lambda,mle[[2]]))
  dgdpi <- c(dGdPi(lambda,mle[[2]]))
  tmpvec <- c(1/dPsi(lambda),dgdpi,1)
  J <- matrix(diag(tmpvec), ncol=(nrow(mle[[2]])+2))
  J[2:(nrow(J)-1),1] <- dgdl

  invJ <- solve(J)
  
  out <- invJ%*%solve(I)%*%t(invJ)
  round(out[-nrow(out), -ncol(out)],5)
}

FIPsi <- function(data, markers, isObserv = FALSE){
  dat <- datasetFormat(data, 2:ncol(data))
  GA <- dat[[3]][markers]
  mle <- MLE(data, markers, allelesName = FALSE)
  if(isObserv){
    I <- FIMObs(data, markers)
  }else{
    I <- FIM(mle, GA)
  }

  lambda <- mle[[1]]
  M <- diag(nrow(mle[[2]])+2)
  M[1,1] <- dPsi(lambda)
  out <- M%*%solve(I)%*%M
  round(out[-nrow(out), -ncol(out)],5)
}

FIM <- function(mle, GA){
  lambda <- mle[[1]]
  freq <- mle[[2]]
  hapnames <- mle[[3]]
  N <- mle[[4]]

  #### Generate all possible observations given GA
  X  <- allObs(GA)
  Nobs <- nrow(X)

  # Find for each infection X all components necessary for the computations
  Ax <- modelSubs(X, GA)[[1]]

  # Initialize the Fisher information matrix of degree d
  hapl <- rank(hapnames, GA)
  rownames(freq) <- hapl
  pp <- matrix(0, ncol=1, nrow = prod(GA)) 
  pp[as.numeric(rownames(freq)),] <- freq

  d <- nrow(freq)+2
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
        sump <- sum(pp[as.numeric(Ax[[u]][[4]][[k]]),])
        vz   <- Ax[[u]][[3]][[k]]
        out1 <- out1 + vz*GFunc(lambda,sump)   # GFunc to build Px
        out2 <- out2 + vz*dGFunc(lambda,sump)  # dG/dl to build dPx/dlam
    }
    if(out1 != 0){
      out <- out + (out2^2)/out1 #Ill
    }
  }
  I[1,1] <- N*out#/(dPsi(lambda)^2)   # I_lam_lam

  # Ipi,pj  new
  out <- 0 * (freq%*%t(freq))
  for(u in 1:Nobs){
    out1 <- 0
    out2 <- 0*pp
    for(k in 1:Ax[[u]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[as.numeric(Ax[[u]][[4]][[k]]),])
      elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1 + vz*GFunc(lambda,sump)
      out2[as.numeric(Ax[[u]][[4]][[k]]),1] <- out2[as.numeric(Ax[[u]][[4]][[k]]),1] + vz*lambda*elp/elmo   # dG/dpi * Indic
    }
    if(out1 != 0){
    out <- out + (out2 %*% t(out2/out1))[c(hapl), c(hapl)]
    }
  }
  I[2:(d-1),2:(d-1)] <- N*out

  # Ipsi,pi
  out <- 0*freq
  for(u in 1:Nobs){
    out1  <- 0
    out21 <- 0
    out22 <- 0*pp
    for(k in 1:Ax[[u]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[as.numeric(Ax[[u]][[4]][[k]]),])
      elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1  + vz*GFunc(lambda,sump)             # GFunc
      out21 <- out21 + vz*dGFunc(lambda,sump)            # dG/dl
      out22[as.numeric(Ax[[u]][[4]][[k]])] <- out22[as.numeric(Ax[[u]][[4]][[k]])] + vz*(lambda*elp/elmo)    # dG/dpi
    }
    if(out1 != 0){
      out <- out + ((out21*out22)/out1)[c(hapl)]
    }
  }
  tmpI <- N*out#/dPsi(lambda)
  I[1,2:(d-1)] <- tmpI
  I[2:(d-1),1] <- tmpI

  # Parameter space in Higher dimension
  dbeta <- c(0,rep(1,nrow(hapl)),0)
  I[,d] <- dbeta
  I[d,] <- dbeta

  I
}

FIMObs <- function(data, markers){
  dat <- datasetFormat(data, 2:ncol(data))
  dataset <- dat[[1]][,markers]
  GA <- dat[[3]][markers]

  ### Dropping Missing data in dataset
  missData <- rowSums(dataset == 0) > 0
  dataset <- dataset[!missData,]

  XNx  <- datasetXNx(dataset)
  X <- XNx[[1]]
  Nx <- XNx[[2]]

  # MLEs
  mle <- MLEBC(dataset, GA, isBC=FALSE, replBC=15000, plugin=NULL)
  lambda <- mle[[1]]
  pp <- mle[[2]]

  # Find for each infection X all components necessary for the computations
  Ax <- modelSubs(X, GA)[[1]]

  # Initialize the Fisher information matrix of degree d
  hapl <- as.numeric(rownames(mle[[2]]))

  d <- nrow(pp)+2
  diag.el <- diag(array(1:d^2,c(d,d)))[-1]
  diag.el <- diag.el[-length(diag.el)]

  I <- matrix(0, nrow=d, ncol=d) #dof
  rownames(I) <- 1:d
  colnames(I) <- 1:d

  # Ilam,lam
  out <- 0
  for(u in seq(length(Nx))){ # For each sampled observation in X
    out1 <- 0
    out2 <- 0
    out3 <- 0
    for(k in 1:Ax[[u]][[1]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump <- sum(pp[Ax[[u]][[4]][[k]],])
      vz   <- Ax[[u]][[3]][[k]]
      out1 <- out1 + vz*GFunc(lambda,sump)   # GFunc to build Px
      out2 <- out2 + vz*dGdL(lambda,sump)    # dG/dl to build dPx/dlam
      out3 <- out3 + vz*d2GdL(lambda,sump)   # d2G/dl2 to build d2Px/dlam2
    }
    if(out1 != 0){
      out <- out + Nx[u]*((out2^2)/out1 - out3)/out1 #Ill
    }
  }
  I[1,1] <- out   # I_lam_lam

  # Ipi,pj  new
  out <- 0 * (pp%*%t(pp))
  for(u in seq(length(Nx))){
    out1 <- 0
    out2 <- 0*pp
    out3 <- 0*pp %*% t(pp)
    for(k in 1:Ax[[u]][[1]][[1]]){  # For each observation y in the sub-observation ScrAx
        sump  <- sum(pp[Ax[[u]][[4]][[k]],])
        vz    <- Ax[[u]][[3]][[k]]
        out1  <- out1 + vz*GFunc(lambda,sump)
        out2[Ax[[u]][[4]][[k]],1]                 <- out2[Ax[[u]][[4]][[k]],1] + vz*dGdPi(lambda, sump)   # dG/dpi * Indici
        out3[Ax[[u]][[4]][[k]],Ax[[u]][[4]][[k]]] <- out3[Ax[[u]][[4]][[k]],Ax[[u]][[4]][[k]]] + vz*d2GdPidPj(lambda,sump)  # d2G/dpidpj * Indici*indicj
    }
    if(out1 != 0){
    out <- out + Nx[u]*((out2 %*% t(out2/out1)) - out3)/out1
    }
  }
  I[2:(d-1),2:(d-1)] <- out

  # Ilam,pi
  out <- 0*pp
  for(u in seq(length(Nx))){
    out1  <- 0
    out21 <- 0
    out22 <- 0*pp
    out3 <- 0*pp
    for(k in 1:Ax[[u]][[1]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[Ax[[u]][[4]][[k]],])
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1  + vz*GFunc(lambda,sump)           # GFunc
      out21 <- out21 + vz*dGdL(lambda,sump)            # dG/dl
      out22[Ax[[u]][[4]][[k]],1] <- out22[Ax[[u]][[4]][[k]],1] + vz*dGdPi(lambda, sump)    # dG/dpi
      out3[Ax[[u]][[4]][[k]],1]  <- out3[Ax[[u]][[4]][[k]],1]  + vz*d2GdLPi(lambda, sump)  # d2GdLPi
    }
    if(out1 != 0){
      out <- out + Nx[u]*(out21*out22/out1 - out3)/out1
    }
  }
  
  I[1,2:(d-1)] <- out
  I[2:(d-1),1] <- out

  # Parameter space in Higher dimension
  dbeta <- c(0,rep(1,length(hapl)),0)
  I[,d] <- dbeta
  I[d,] <- dbeta

  I
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

rank <- function(hap, nloci){
hap <- hap-1
gk <- cumprod(nloci)
rk <- hap%*%c(1,gk)[-(length(gk)+1)] + 1
rk
}

prev0 <- function(mle){
  lambda <- mle[[1]]
  pp <- mle[[2]]
  lp <- lambda*pp
  mle[[2]][,1] <- (1-exp(-lp))/(1-exp(-lambda))
  mle
}

PREV <- function(data, markers, idExists=TRUE, plugin=NULL, isCI=FALSE, replCI=10000, alpha=0.05, allelesName=TRUE){
  data <- datasetFormat(data, 2:ncol(data))
  dataset <- data[[1]][,markers]
  alleleList <- data[[2]][markers]
  GA <- data[[3]][markers]
  
  ### Dropping Missing data in dataset
  missData <- rowSums(dataset == 0) > 0
  dataset <- dataset[!missData,]

  ### Actual sample size
  Neff <- nrow(dataset)
  XNx  <- datasetXNx(dataset)

  ### Dropping Missing data in dataset
  missData <- rowSums(dataset == 0) > 0
  dataset <- dataset[!missData,]
  X <- XNx[[1]]
  Nx <- XNx[[2]]
  nLoci <- ncol(X)

  mle <- MLEPluginChoice(XNx, GA, plugin=NULL)
  mle <- prev0(mle)

  mixRadCumProd <- c(1, cumprod(GA)[1:(nLoci-1)])
  prevEstimates <- mle[[2]]
  lbdaEstim <- mle[[1]]
  rnames1 <- as.integer(rownames(prevEstimates)) - 1
  rnames  <- rnames1
  nh <- length(rnames)
  mixRadixTableOfHap <- array(0,c(nh,nLoci))
  for(k in 1:nLoci){ #for each locus
    re <- rnames%%rev(mixRadCumProd)[k]
    mixRadixTableOfHap[,(nLoci-k+1)] <- as.character((rnames-re)/rev(mixRadCumProd)[k]+1)
    rnames <- re
  }
  #mixRadixTableOfHap <- mixRadixTableOfHap+1
  if(allelesName){
    prevEstimates <- hapName(prevEstimates,alleleList)
    for (i in 1:length(alleleList)){
      mixRadixTableOfHap[,i] <- alleleList[[i]][as.numeric(mixRadixTableOfHap[,i])]
    }
  }else {
    for(i in 1:nh){
      rnames[i] <- paste(mixRadixTableOfHap[i,], collapse = '.')
    }
    rownames(prevEstimates) <- rnames
    mixRadixTableOfHap <- apply(mixRadixTableOfHap, 2, as.numeric)
  }
  colnames(mixRadixTableOfHap) <- names(GA)

  # Bootstrap CIs
  if(isCI){
    nHap  <- length(prevEstimates)
    N <- sum(Nx)
    probaObs  <- Nx/N
    arrayBootEstim <- array(0, dim = c((nHap+1), replCI=replCI))
    rownames(arrayBootEstim) <- c('lambda',(rnames1+1))
    #observation <- seq_along(Nx)
    for (bootRepl in 1:replCI){
      bootstrappedListOfDatasets <- bootstrapDataset(X, N, probaObs)
      bootEstim <- prev0(MLEPluginChoice(bootstrappedListOfDatasets, GA, plugin=plugin))
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
    prevEstimCI <- cbind(prevEstimates,perc[2:(nHap+1),])
    colnames(prevEstimCI) <- c('', paste0(as.character((alpha/2)*100), '%'), paste0(as.character((1-alpha/2)*100), '%'))
    mle <- list(lbdaEstim, prevEstimCI, mixRadixTableOfHap, Neff)
    names(mle) <- c('lambda', 'haplotypes_prevalence', 'detected_haplotypes', 'used_sample_size')
  }else{
    mle <- list(lbdaEstim, prevEstimates, mixRadixTableOfHap, Neff)
    names(mle) <- c('lambda', 'haplotypes_prevalence', 'detected_haplotypes', 'used_sample_size')
  }
  rownames(mle[[2]]) <- paste("q", rownames(mle[[2]]), sep="") 
  mle
}

dGdPi <- function(lambda, pp){
  elmo <- exp(lambda) - 1
  lambda*exp(lambda*pp)/elmo
}

dGdL <- function(lambda, pp){
  pp*exp(-lambda*pp)/(1-exp(-lambda))-exp(lambda)*(1-exp(-lambda*pp))/(1-exp(-lambda))^2
}

d2GdL <- function(lambda, pp){
  elmo <- exp(lambda) - 1
  pp^2*exp(lambda*pp)/elmo - 2*pp*exp(lambda*(pp+1))/elmo^2 - exp(lambda)*(exp(lambda*pp)-1)/elmo^2 + 2*exp(2*lambda)*(exp(lambda*pp)-1)/elmo^3
}

d2GdPidPj <- function(lambda, pp){
  lambda^2*exp(lambda*pp)/(exp(lambda)-1)
}

d2GdLPi <- function(lambda, pp){
  elmo <- exp(lambda) - 1
  ((lambda*pp + 1)*exp(lambda*pp)/elmo - lambda*exp(lambda*(pp+1))/(elmo^2))
}

dataGen <- function(P,lambda, N, GA){ 
  H <- allHap(GA)                   
  out <- matrix(0, nrow=N, ncol=length(GA))
  m <- cPoiss(lambda, N)              
  for(j in 1:N){                     
    s <- rmultinom(1, m[j], P) 
    out[j,] <- observ(H[s!=0,], GA) 
  }
  out
}

observ <- function(M, GA){
  M <- data.frame(matrix(M, ncol=length(GA)))
  alleles <- lapply(GA, function(x) seq(0, x-1)) 
  xx <- lapply(M,function(x) sort(unique(x)))
  detectAllel <- Map(function(x,y) is.element(y,x), xx ,alleles) 
  yy <- lapply(GA, function(x) 2^(0:(x-1))) 

  mapply('%*%',detectAllel,yy)
}
