# Title        : Simulation study using simulated data
# Objectives   : Implement the EM-algorithm on simulated data and save the estimates
# Created by   : christian Tsoungui Obama, Kristan A. Schneider
# Created on   : 10.08.23
# Last modified: 04.10.23

# Import library
#library(matrixcalc)

psi_estim <- function(la){
  la/(1-exp(-la))
}

baseModelSim <- function(data,GA){
  
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
    pp[,1] <- Nx/N 
    la <- 0
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
  out <- list(la,matrix(unlist(pp),ncol = 1))
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

baseModelSim2 <- function(XNx,GA){
  
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

FIPsiSim <- function(XNx, GA, isObserv = FALSE){
  mle <- baseModelSim(XNx, GA)
  mle[[3]] <- as.numeric(names(mle[[2]]))
  mle[[4]] <- sum(XNx[[2]])
  if(isObserv){
    I <- FIMObs(data, markers)
  }else{
    I <- FIMSim(mle, GA)
  }

  lambda <- mle[[1]]
  M <- diag(length(mle[[2]])+2)
  M[1,1] <- dPsi(lambda)
  out <- M%*%solve(I)%*%M
  out <- out[-nrow(out), -ncol(out)]
  diag(out)
}

FIMSim <- function(mle, GA){
  lambda <- mle[[1]]
  freq <- matrix(mle[[2]],ncol=1)
  hapl <- mle[[3]]
  N <- mle[[4]]

  #### Generate all possible observations given GA
  X  <- allObs(GA)
  Nobs <- nrow(X)

  # Find for each infection X all components necessary for the computations
  Ax <- modelSubs(X, GA)[[1]]

  # Initialize the Fisher information matrix of degree d
  #hapl <- rank(hapnames, GA)
  rownames(freq) <- hapl
  pp <- matrix(0, ncol=1, nrow = prod(GA)) 
  pp[as.numeric(rownames(freq)),] <- freq

  d <- nrow(freq)+2
  diag.el <- diag(array(1:d^2,c(d,d)))[-1]
  diag.el <- diag.el[-length(diag.el)]

  I <- matrix(0, nrow=d, ncol=d) #dof
  rownames(I) <- 1:d
  colnames(I) <- 1:d

  #elmo <- exp(lambda)-1

  # Ipsi,psi
  out <- 0
  for(u in 1:Nobs){ # For each observation x in ScrO
    out1 <- 0
    out2 <- 0
    for(k in 1:Ax[[u]][[1]][[1]]){  # For each observation y in the sub-observation ScrAx
        sump <- sum(pp[as.numeric(Ax[[u]][[4]][[k]]),])
        vz   <- Ax[[u]][[3]][[k]]
        out1 <- out1 + vz*G(lambda,sump)   # GFunc to build Px
        out2 <- out2 + vz*dGdL(lambda,sump)  # dG/dl to build dPx/dlam
    }
    if(!is.nan(out1)){
      if(out1 != 0){
        out <- out + (out2^2)/out1 #Ill
      }
    }
  }
  I[1,1] <- N*out#/(dPsi(lambda)^2)   # I_lam_lam

  # Ipi,pj  new
  out <- 0 * (freq%*%t(freq))
  for(u in 1:Nobs){
    out1 <- 0
    out2 <- 0*pp
    for(k in 1:Ax[[u]][[1]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[as.numeric(Ax[[u]][[4]][[k]]),])
      #elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1 + vz*G(lambda,sump)
      out2[as.numeric(Ax[[u]][[4]][[k]]),1] <- out2[as.numeric(Ax[[u]][[4]][[k]]),1] + vz*dGdPi(lambda, sump)# lambda*elp/elmo   # dG/dpi * Indic
    }
    if(!is.nan(out1)){
      if(out1 != 0){
      out <- out + (out2 %*% t(out2/out1))#[c(hapl), c(hapl)]
      }
    }
  }
  I[2:(d-1),2:(d-1)] <- N*out

  # Ipsi,pi
  out <- 0*freq
  for(u in 1:Nobs){
    out1  <- 0
    out21 <- 0
    out22 <- 0*pp
    for(k in 1:Ax[[u]][[1]][[1]]){  # For each observation y in the sub-observation ScrAx
      sump  <- sum(pp[as.numeric(Ax[[u]][[4]][[k]]),])
     # elp   <- exp(lambda*sump)
      vz    <- Ax[[u]][[3]][[k]]
      out1  <- out1  + vz*G(lambda,sump)             # G
      out21 <- out21 + vz*dGdL(lambda,sump)            # dG/dl
      out22[as.numeric(Ax[[u]][[4]][[k]])] <- out22[as.numeric(Ax[[u]][[4]][[k]])] + vz*dGdPi(lambda, sump)#(lambda*elp/elmo)    # dG/dpi
    }
    if(!is.nan(out1)){
      if(out1 != 0){
        out <- out + ((out21*out22)/out1)#[c(hapl)]
      }
    }
  }
  tmpI <- N*out#/dPsi(lambda)
  I[1,2:(d-1)] <- tmpI
  I[2:(d-1),1] <- tmpI

  # Parameter space in Higher dimension
  dbeta <- c(0,rep(1,length(hapl)),0)
  I[,d] <- dbeta
  I[d,] <- dbeta

  I
}

# Relative path
path <- "/Users/christian/Documents/phd/models/generalModel/"

# Loading external ressources
source(paste0(path, "src/MOI-MLE.R"))          ## Loading Model
#source(paste0(path, "forSimulation/dataGenerator.R"))  

# True Poisson parameter
lbdavec <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5)
NLbd <- length(lbdavec)

# Genetic architecture
GA <- matrix(c(2,2,4,7), ncol = 2, byrow = TRUE)
Nn <- nrow(GA)
#nvec <- matrix(NumbLoci, nrow = 1, ncol = Nn)

# Number of possible haplotypes
Hvec <-  apply(GA,1, prod) 
NH <- length(Hvec)
Hvecpo <- Hvec + 1

# Sample sizes considered
Nvec <- c(50, 100, 150, 200, 500) #c(500)
NN   <- length(Nvec)

# Number of estimates generated in the simulation
NEst <- 100000

# Number of frequencies set per (number of loci) case
NFreq <- 2

# Extra parameters
ParExtra <- list(NLbd, Nn, Hvec, NN, NEst, NFreq, GA)

# True haplotype frequencies
Pvec <- vector(mode="list", length=Nn)
for (i in 1:Nn){
  Pvec[[i]] <- array(1/Hvec[i], c(1, Hvec[i]))
  Pvec[[i]] <- rbind(Pvec[[i]], c(0.70, rep(0.30/(Hvec[i]-1), (Hvec[i]-1))))
}

# True parameter
True_param <- list(Pvec, lbdavec, Nvec)

# Simulation
out  <- vector(mode = "list", length = Nn)

for (i in 1:2){
  print(paste0("processing frequency distributions for m=", GA[i,1], " and n=",  GA[i,2], " alleles, respectively."))
  sizelist <- vector(mode = "list", length = NN)
  for (j in 1:length(Nvec)){                               ## For each value of the sample size
    lbdalist <- vector(mode = "list", length = NLbd)
    for (k in 1:NLbd){                                     ## For each value of the lambda parameter
      Estim <- array(0, dim = c(Hvecpo[i], NEst, NFreq))
      print(c(Nvec[j], lbdavec[k]))
      for (cnt in 1:NFreq){
        for (l in 1:NEst){
          infct   <- dataGen(unlist(Pvec[[i]][cnt,]) ,unlist(lbdavec[k]) ,Nvec[j], GA[i,])      ## Generating data for the simulation
          XNx <- datasetXNx(infct)
          #est.tmp <- baseModel(XNx, GA[i,])
          I     <- FIPsiSim(XNx, GA[i,])
          tmp   <- (sqrt(I)/c(psi_estim(lbdavec[k]), Pvec[[i]][cnt,]))*100
          Estim[as.numeric(colnames(I)),l,cnt] <- tmp      ## Coefficient of variation in %
        }
      }
      lbdalist[[k]]  <- Estim
    }
    sizelist[[j]] <- lbdalist
    #saveRDS(lbdalist, file = paste0(path, "dataset/full_modelEstimates_", i,"_Sample_", Nvec[j], ".rds"))
  }
  saveRDS(sizelist, file = paste0(path, "dataset/full_modelEstimates_", i, ".rds"))
  # End of simulation warning
  print(paste0("Finished simulation frequency distributions for m=", GA[i,1], " and n=",  GA[i,2], " alleles, respectively."))
}

# Saving the true parameters
saveRDS(True_param, file = paste0(path, "dataset/full_true_Parameters.rds"))

# Saving the extra parameters
saveRDS(ParExtra, file = paste0(path, "dataset/full_extra_Parameters.rds"))