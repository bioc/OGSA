#'copaStat
#'
#'Calculates outlier statistics by the Tibshirani-Hastie method
#'@usage copaStat (data, phenotype, tail='right', perms=100, permType='array')
#'@param data A matrix of nGene by nSample
#'@param phenotype A vector of 0s and 1s of length nSample, where 1 = case, 
#'0 = control
#'@param tail Indicates whether outliers are up (right) or down (left) outliers
#'@param perms The number of permutations
#'@param permType By all on array or by gene, if by gene increase perms 
#'significantly and plan on lots of time; in theory array should be fine as 
#'genes are rescaled
#'@return A vector with outlier counts by gene
#'@examples
#'
#'data(ExampleData)
#'
#' #Set up phenotype
#' phenotype <- pheno  
#' names(phenotype) <- colnames(cnv)
#'
#' #set up values for expr-meth-cnv in that order
#' tailLRL <- c('left', 'right', 'left')   
#' 
#' #setup dataList
#' dataSet <- list(expr, meth, cnv)
#'
#' data <- dataSet[[1]]
#' 
#' tibL <- copaStat(data, phenotype, tail='right', perms=100, permType='array')
#' 
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@export

copaStat <- function(data, phenotype, tail='right', perms=100,
										permType='array') {
    
    ##function input check
    if (all(phenotype < 0 & phenotype > 1 & phenotype%%1==0)) {
        stop("elements in phenotype must be 0 or 1")
    }
    
    nG <- dim(data)[1]
    nS <- length(phenotype)
    nT <- sum(phenotype)

    nData <- data[, phenotype==0]

    dataC <- matrix(nrow=nG, ncol=nS)

    copaValue <- vector(length=nG)
    copaStat <- vector(length=nG)
    copaPerm <- matrix(nrow=nG, ncol=perms)

    # generate permutation COPA values
    for (j in 1:perms) {
        phenoRand <- rep(0, nS)
        phenoRand[sample(1:nS, nT)] <- 1

        # generate Tibshirani eqn 2.2
        # putting genes on same scale
        nData <- data[, phenoRand==0]
        medG <- apply(data, 1, median)
        madGN <- apply(nData, 1, mad)
        
        for (i in 1:nG) {
            temp <- (data[i, ] - medG[i])/madGN[i]
            dataC[i, ] <- temp
        }

        # calc outlier stat on scaled data
        tData <- dataC[, phenoRand==1]
        nData <- dataC[, phenoRand==0]

        iqrG <- apply(dataC, 1, IQR)
        if (tail == 'right') {
            quantG <- apply(dataC, 1, quantile, probs = 0.75)
        } else if (tail == 'left') {
            quantG <- apply(dataC, 1, quantile, probs = 0.25)
        }

        if (tail == 'right') {
            for (i in 1:nG) {
                Ind <- tData[i, ] > quantG[i] + iqrG[i]
                copaValue[i] <- sum(tData[i, ] * Ind)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                Ind <- tData[i, ] < quantG[i] - iqrG[i]
                copaValue[i] <- sum(tData[i, ] * Ind)
            }
        }
        copaPerm[, j] <- copaValue
    }
    nPerm <- dim(copaPerm)[1] * dim(copaPerm)[2]

    # Now for real data
    # generate Tibshirani eqn 2.2
    # putting genes on same scale
    medG <- apply(data, 1, median)
    madGN <- apply(nData, 1, mad)

    for (i in 1:nG) {
        temp <- (data[i, ] - medG[i])/madGN[i]
        data[i, ] <- temp
    }
    
    # generate stat on scaled data
    tData <- data[, phenotype==1]
    nData <- data[, phenotype==0]
    iqrG <- apply(data, 1, IQR)
    if (tail == 'right') {
        quantG <- apply(data, 1, quantile, probs = 0.75)
    } else if (tail == 'left') {
        quantG <- apply(data, 1, quantile, probs = 0.25)
    }
    
    if (tail == 'right') {
        for (i in 1:nG) {
            Ind <- tData[i, ] > quantG[i] + iqrG[i]
            copaValue[i] <- sum(tData[i, ] * Ind)
        }
    } else if (tail == 'left') {
        for (i in 1:nG) {
            Ind <- tData[i, ] < quantG[i] - iqrG[i]
            copaValue[i] <- sum(tData[i, ] * Ind)
        }
    }

    # generate comparison of COPA values to permutations
    if(permType == 'array') {
        nPerm <- dim(copaPerm)[1] * dim(copaPerm)[2]
        if (tail == 'right') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm > copaValue[i])+1)/(nPerm+1)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm < copaValue[i])+1)/(nPerm+1)
            }
        }
    } else if (permType == 'gene') {
        nPerm <- dim(copaPerm)[2]
        if (tail == 'right') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm[i, ] > copaValue[i])+1)/(nPerm+1)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                copaStat[i] <- (sum(copaPerm < copaValue[i])+1)/(nPerm+1)
            }
        }
    } else {
        print("Unknown Permutation Comparison Type")
    }

    names(copaStat) <- rownames(data)

    return(copaStat)
}