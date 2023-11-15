
# varsets <- function(l,n){   #calculate all var sets
#   # n number of loci
#   # l number of alleles per locus
#   B <- array(0,c(l^n,n))
#   B[1:l,1] <- 0:(l-1)
#   lkmo <- l
#   if(n>1){
#     for(k in 2:n){
#       lk <- lkmo*l
#       pick1 <- (lkmo+1):lk
#       B[pick1,] <- B[rep(1:lkmo,l-1),]
#       B[pick1,k] <- rep(1:(l-1),each=lkmo)
#       lkmo <- lk
#     }
#   }
#   B
# }


# gead1 <- function(x,l){   ## calculates general geadic expression of each element of vector x
#   n <- length(l)
#   out <- array(0,c(length(x),n))
#   div <- c(1,cumprod(l[1:(n-1)]))
#   for(k in n:1){
#     r <- x%%div[k]
#     out[,k] <- (x-r)/div[k]
#     x <- r
#   }
#   out
# }

# Base model

# List of components of each observation x necessary for computations

# Model with plugin value for lambda


#################################
# The function calculateMaximuLikelihoodEstimatesWithBiasCorrection(dataset, numberOfAllelesAtEachMarker) implements the EM algorithm with the option for bias correction (isBiasCorrection) using either
# a "bootstrap", or a "Jacknife" method, and returns the MLEs, i.e., estimates of haplotype frequencies and Poisson parameter.
#################################




#################################
# The function mle() wraps the reform(X1,id) and either calculateMaximuLikelihoodEstimatesWithBiasCorrection(dataset, numberOfAllelesAtEachMarker) or calculateMaximumLikelihoodEstimatesWithPluginValueOfLambda(dataset, numberOfAllelesAtEachMarker, pluginValueOfLambda) to find the MLEs
# with or without the Poisson parameter as a plug-in estimate, respectively. Moreover, the option to ouput the bias corrected (isBiasCorrection) estimates with
# confidence intervals (isConfidenceInterval) is available. The function outputs the estimates for haplotype frequencies, Poisson parameters, and a matrix of detected haplotypes.
#################################

#################################
# Function buildAllPossibleHaplotypes(numberOfAllelesAtEachMarker) takes as input a vector of the number of alleles per locus and outputs
#  a matrix of all possible haplotypes in geadic representation
#################################

#################################
# Function buildAllPossibleObservations(numberOfAllelesAtEachMarker) takes as input a vector of the number of alleles per locus and outputs
#  a matrix of all possible observations in geadic representation
#################################

#################################
# Function sampleWithConditionalPoisson(lambda,numberOfLoci) outputs n randomly drawn integer from a condidtional Poisson distribution
# with parameter lambda
#################################

### function convertDatasetToStandardFormat outputs a list containing as first element  the data in new format (1 row per sample, 1 colum by marker)
### entries are integerers. If transformed into binary numbers 0-1 vectors they indicate absence/presence of alleles
### e.g., 0.... no allele present, 1... first allele present, 2... 2nd allele present, 3.... 1st and 2nd allele present ect.
### The second element gives the order of the alleles pper locus
### The third element gives th number of alleles per locucs

# labels output frequencies correctly


# #################################
# # The function obs(M) gives a representation for an infection with haplotypes given by the matrix M.
# # The input M is  k x n matrix with entries 0 and 1, where each row  is a haplotype corresponding to
# # a 0-1 vector. Obs returns the corresponding vector representation of the buildAllPossibleObservations
# #################################
# obs <- function(M, numberOfAllelesAtEachMarker){
#   M <- matrix(M, ncol = 2)
#   out <- array(0,2)
#   arch1 <- numberOfAllelesAtEachMarker - 1
#   x <- list()
#   x[[1]] <- seq(0,arch1[1])
#   x[[2]] <- seq(0,arch1[2])

#   for(i in 1:2){
#     bin <- 2^x[[i]]
#     binx <- is.element(x[[i]], M[,i])
#     out[i] <- bin %*% binx
#   }
#   out
# }

# #################################
# # Function datasetgen(P,lambda,sampleSize,n) is used to simulate data. It generates sampleSize observations assuming n biallelic loci with haplotype distribution P
# # wich must be a vector of length 2n a sampleSize x n matrix of observations sampled using the multinomial and Poisson distribution
# # of parameters (m, P) and lambda respectively. m is the MOI for the corresponding sample.
# #################################
# datasetgen <- function(P,lambda,sampleSize,numberOfAllelesAtEachMarker){
#   # This function simulates the data as a Nx2 matrix of infections

#   H <- buildAllPossibleHaplotypes(numberOfAllelesAtEachMarker)                    # Set of possible haplotypes
#   out <- matrix(0,nrow=sampleSize, ncol=2)
#   m <- sampleWithConditionalPoisson(lambda,sampleSize)              # MOI values for each sample following CPoiss(lambda)
#   for(j in 1:sampleSize){
#     s <- rmultinom(1, m[j], P) #multinomially select M[j] haplotypes from the haplotype pool
#     out[j,] <- obs(H[s!=0,], numberOfAllelesAtEachMarker) #Summing up the mixed-radix representation of a number representing the infection
#   } #vector of infections
#   out
# }

# datagen <- function(P,lambda,sampleSize,numberOfAllelesAtEachMarker){
#   # This function generates the data as a list containing the observations detected (detectedObservations)
#   # and the number of times each infection is detected (numberOfEachDetectedObservations
#   # Output: list(detectedObservations, numberOfEachDetectedObservations)
#   out <- datasetgen(P,lambda,sampleSize,numberOfAllelesAtEachMarker)
#   reform(out, numberOfAllelesAtEachMarker, idExists = FALSE)
# }

### This function outputs a list
#### 1st gives the compact notation for the data
#### 2nd element base factos 1 g1 g1g2. ....,g1*...g(l-1)
#### 3rd element has the number of alleles in the prope order per locus
#### 4th element is the number of alleles per locus
#dafa.formal.AL <- function(data){
#        l <- length(data[[3]])   ## number of loci
#        basefactors <- cumprod(c(1,2^(data[[3]][-l])))
#        list(data[[1]] %*% basefactors, basefactors,data[[2]],data[[3]])
#}

### Going back

#div <- function(detectedObservations,f){
#    x <- NULL
#    for(f in rev(dataset[[2]])){
#        #print(c(detectedObservations,f))
#        x <- c(detectedObservations %/% f,x)
#        detectedObservations <- detectedObservations%%f
#    }
#    x
#}

#reconst.data <- function(data){
#    newdat <- t(sapply(data[[1]], function(x) div(x,data[[2]]) ) )
#    list(newdat,data[[3]],data[[4]])
#}

#_____________________________
# labels output frequencies correctly
#labelFrequencyEstimates <- function(pp,allist){ ## allist list with alleles per locos
#  n <- length(allist)
#  allnum <- unlist(lapply(allist,length))
#  hapl <- gead1(as.numeric(rownames(pp))-1,allnum)+1
#  hapl1 <- array(,dim(hapl))
#  for(l in 1: ncol(hapl)){
#    for(m in 1: nrow(hapl)){
#      hapl1[m,l] <- allist[[l]][hapl[m,l]]
#    }
#  }
#  hapl1 <- apply(hapl1,1,function(x) paste(x,sep="",collapse="-"))
#  rownames(pp) <- hapl1
#  pp
#}

#_____________________________
# plot of frequencies
# freqplot <- function(out,out1,lv){
#   al <- attr(out[[5]],"names") #alleles in 1st data set
#   al1 <- attr(out1[[5]],"names") #alleles in 2st data set
#   alleles <- unique(c(al,al1)) # alleles in both data sets
#
#   pldata <- array(0,c(length(alleles),2)) # create array that contains frequencies
#
#   rownames(pldata) <- alleles # the rows are labelled with the STR repet lengths
#   colnames(pldata) <- lv  # the colums are the year intervals
#   pldata[al,1] <- out[[3]]  # in the first column the alllele frequencies of the old data
#   pldata[al1,2] <- out1[[3]]
#
#   pldata <- melt(pldata,id.vars="all")  # reshapes the data 3 colums Var1, Var2, value
#
#   # next transform plotdata into data frame
#   pldata <- data.frame(allele=as.factor(pldata$Var1),gr=as.factor(pldata$Var2),freq=as.numeric(pldata$value))
#
#   # Define colors for the plot - always be coulorbliend friendly
#   #cbPalette <- c("#3ffdcc","#fa6403","#cdff00","#1e9864","#999999", "#E69F00")
#   cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
#   #plot
#   p <- ggplot(data=pldata, aes(x=allele, y=freq, fill=gr)) +
#       geom_bar(stat="identity",position=position_dodge(), colour="black")
#   p <- p + theme(panel.background = element_blank(),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  panel.border = element_rect(colour='black',fill=NA),
#                  axis.text = element_text(size = rel(1.3),color='black'),
#                  axis.text.x = element_text(angle=75,vjust=0.5),
#                  axis.title = element_text(size = rel(1.3)),
#                  plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
#                  legend.text = element_text(size = rel(1.3)),
#                  legend.title = element_text(size = rel(1.3)))
#   p <- p + scale_fill_manual(values=cbPalette,name="Years") + ylim(0,1)
#   p <- p + labs(x="alleles",y="frequencies",title=parse(text=out[[6]]))
#   p
# }


# varsets <- function(l,n){   #calculate all var sets
#   # n number of loci
#   # l number of alleles per locus
#   B <- array(0,c(l^n,n))
#   B[1:l,1] <- 0:(l-1)
#   lkmo <- l
#   if(n>1){
#     for(k in 2:n){
#       lk <- lkmo*l
#       pick1 <- (lkmo+1):lk
#       B[pick1,] <- B[rep(1:lkmo,l-1),]
#       B[pick1,k] <- rep(1:(l-1),each=lkmo)
#       lkmo <- lk
#     }
#   }
#   B
# }



# gead1 <- function(x,l){   ## calculates general geadic expression of each element of vector x
#   n <- length(l)
#   out <- array(0,c(length(x),n))
#   div <- c(1,cumprod(l[1:(n-1)]))
#   for(k in n:1){
#     r <- x%%div[k]
#     out[,k] <- (x-r)/div[k]
#     x <- r
#   }
#   out
# }

# Base model

# List of components of each observation x necessary for computations

# Model with plugin value for lambda


#################################
# The function calculateMaximuLikelihoodEstimatesWithBiasCorrection(dataset, numberOfAllelesAtEachMarker) implements the EM algorithm with the option for bias correction (isBiasCorrection) using either
# a "bootstrap", or a "Jacknife" method, and returns the MLEs, i.e., estimates of haplotype frequencies and Poisson parameter.
#################################




#################################
# The function mle() wraps the reform(X1,id) and either calculateMaximuLikelihoodEstimatesWithBiasCorrection(dataset, numberOfAllelesAtEachMarker) or calculateMaximumLikelihoodEstimatesWithPluginValueOfLambda(dataset, numberOfAllelesAtEachMarker, pluginValueOfLambda) to find the MLEs
# with or without the Poisson parameter as a plug-in estimate, respectively. Moreover, the option to ouput the bias corrected (isBiasCorrection) estimates with
# confidence intervals (isConfidenceInterval) is available. The function outputs the estimates for haplotype frequencies, Poisson parameters, and a matrix of detected haplotypes.
#################################

#################################
# Function buildAllPossibleHaplotypes(numberOfAllelesAtEachMarker) takes as input a vector of the number of alleles per locus and outputs
#  a matrix of all possible haplotypes in geadic representation
#################################

#################################
# Function buildAllPossibleObservations(numberOfAllelesAtEachMarker) takes as input a vector of the number of alleles per locus and outputs
#  a matrix of all possible observations in geadic representation
#################################


#################################
# Function sampleWithConditionalPoisson(lambda,numberOfLoci) outputs n randomly drawn integer from a condidtional Poisson distribution
# with parameter lambda
#################################

### function convertDatasetToStandardFormat outputs a list containing as first element  the data in new format (1 row per sample, 1 colum by marker)
### entries are integerers. If transformed into binary numbers 0-1 vectors they indicate absence/presence of alleles
### e.g., 0.... no allele present, 1... first allele present, 2... 2nd allele present, 3.... 1st and 2nd allele present ect.
### The second element gives the order of the alleles pper locus
### The third element gives th number of alleles per locucs

# labels output frequencies correctly


# #################################
# # The function obs(M) gives a representation for an infection with haplotypes given by the matrix M.
# # The input M is  k x n matrix with entries 0 and 1, where each row  is a haplotype corresponding to
# # a 0-1 vector. Obs returns the corresponding vector representation of the buildAllPossibleObservations
# #################################
# obs <- function(M, numberOfAllelesAtEachMarker){
#   M <- matrix(M, ncol = 2)
#   out <- array(0,2)
#   arch1 <- numberOfAllelesAtEachMarker - 1
#   x <- list()
#   x[[1]] <- seq(0,arch1[1])
#   x[[2]] <- seq(0,arch1[2])

#   for(i in 1:2){
#     bin <- 2^x[[i]]
#     binx <- is.element(x[[i]], M[,i])
#     out[i] <- bin %*% binx
#   }
#   out
# }

# #################################
# # Function datasetgen(P,lambda,sampleSize,n) is used to simulate data. It generates sampleSize observations assuming n biallelic loci with haplotype distribution P
# # wich must be a vector of length 2n a sampleSize x n matrix of observations sampled using the multinomial and Poisson distribution
# # of parameters (m, P) and lambda respectively. m is the MOI for the corresponding sample.
# #################################
# datasetgen <- function(P,lambda,sampleSize,numberOfAllelesAtEachMarker){
#   # This function simulates the data as a Nx2 matrix of infections

#   H <- buildAllPossibleHaplotypes(numberOfAllelesAtEachMarker)                    # Set of possible haplotypes
#   out <- matrix(0,nrow=sampleSize, ncol=2)
#   m <- sampleWithConditionalPoisson(lambda,sampleSize)              # MOI values for each sample following CPoiss(lambda)
#   for(j in 1:sampleSize){
#     s <- rmultinom(1, m[j], P) #multinomially select M[j] haplotypes from the haplotype pool
#     out[j,] <- obs(H[s!=0,], numberOfAllelesAtEachMarker) #Summing up the mixed-radix representation of a number representing the infection
#   } #vector of infections
#   out
# }

# datagen <- function(P,lambda,sampleSize,numberOfAllelesAtEachMarker){
#   # This function generates the data as a list containing the observations detected (detectedObservations)
#   # and the number of times each infection is detected (numberOfEachDetectedObservations
#   # Output: list(detectedObservations, numberOfEachDetectedObservations)
#   out <- datasetgen(P,lambda,sampleSize,numberOfAllelesAtEachMarker)
#   reform(out, numberOfAllelesAtEachMarker, idExists = FALSE)
# }

### This function outputs a list
#### 1st gives the compact notation for the data
#### 2nd element base factos 1 g1 g1g2. ....,g1*...g(l-1)
#### 3rd element has the number of alleles in the prope order per locus
#### 4th element is the number of alleles per locus
#dafa.formal.AL <- function(data){
#        l <- length(data[[3]])   ## number of loci
#        basefactors <- cumprod(c(1,2^(data[[3]][-l])))
#        list(data[[1]] %*% basefactors, basefactors,data[[2]],data[[3]])
#}

### Going back

#div <- function(detectedObservations,f){
#    x <- NULL
#    for(f in rev(dataset[[2]])){
#        #print(c(detectedObservations,f))
#        x <- c(detectedObservations %/% f,x)
#        detectedObservations <- detectedObservations%%f
#    }
#    x
#}

#reconst.data <- function(data){
#    newdat <- t(sapply(data[[1]], function(x) div(x,data[[2]]) ) )
#    list(newdat,data[[3]],data[[4]])
#}

#_____________________________
# labels output frequencies correctly
#labelFrequencyEstimates <- function(pp,allist){ ## allist list with alleles per locos
#  n <- length(allist)
#  allnum <- unlist(lapply(allist,length))
#  hapl <- gead1(as.numeric(rownames(pp))-1,allnum)+1
#  hapl1 <- array(,dim(hapl))
#  for(l in 1: ncol(hapl)){
#    for(m in 1: nrow(hapl)){
#      hapl1[m,l] <- allist[[l]][hapl[m,l]]
#    }
#  }
#  hapl1 <- apply(hapl1,1,function(x) paste(x,sep="",collapse="-"))
#  rownames(pp) <- hapl1
#  pp
#}

#_____________________________
# plot of frequencies
# freqplot <- function(out,out1,lv){
#   al <- attr(out[[5]],"names") #alleles in 1st data set
#   al1 <- attr(out1[[5]],"names") #alleles in 2st data set
#   alleles <- unique(c(al,al1)) # alleles in both data sets
#
#   pldata <- array(0,c(length(alleles),2)) # create array that contains frequencies
#
#   rownames(pldata) <- alleles # the rows are labelled with the STR repet lengths
#   colnames(pldata) <- lv  # the colums are the year intervals
#   pldata[al,1] <- out[[3]]  # in the first column the alllele frequencies of the old data
#   pldata[al1,2] <- out1[[3]]
#
#   pldata <- melt(pldata,id.vars="all")  # reshapes the data 3 colums Var1, Var2, value
#
#   # next transform plotdata into data frame
#   pldata <- data.frame(allele=as.factor(pldata$Var1),gr=as.factor(pldata$Var2),freq=as.numeric(pldata$value))
#
#   # Define colors for the plot - always be coulorbliend friendly
#   #cbPalette <- c("#3ffdcc","#fa6403","#cdff00","#1e9864","#999999", "#E69F00")
#   cbPalette <- c( "#0072B2", "#E69F00" , "#009E73","#56B4F9", "#CC79A7")
#   #plot
#   p <- ggplot(data=pldata, aes(x=allele, y=freq, fill=gr)) +
#       geom_bar(stat="identity",position=position_dodge(), colour="black")
#   p <- p + theme(panel.background = element_blank(),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  panel.border = element_rect(colour='black',fill=NA),
#                  axis.text = element_text(size = rel(1.3),color='black'),
#                  axis.text.x = element_text(angle=75,vjust=0.5),
#                  axis.title = element_text(size = rel(1.3)),
#                  plot.title = element_text(size = rel(1.4),color='black',hjust=0.5),
#                  legend.text = element_text(size = rel(1.3)),
#                  legend.title = element_text(size = rel(1.3)))
#   p <- p + scale_fill_manual(values=cbPalette,name="Years") + ylim(0,1)
#   p <- p + labs(x="alleles",y="frequencies",title=parse(text=out[[6]]))
#   p
# }
