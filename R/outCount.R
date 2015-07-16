#'outCount
#'
#'Counts outliers by the Tibshirani and Hastie method. Adds the ability to 
#'subtract for outliers in the normals using corr = TRUE
#'@usage outCount (data, phenotype, tail='right', corr=FALSE)
#'@param data A matrix of nGene by nSample
#'@param phenotype A vector of 0s and 1s of length nSample, where 1 = case, 
#'0 = control
#'@param tail Indicates whether outliers are up (right) or down (left) outliers
#'@param corr Whether to correct for normal outliers
#'@return A vector with outlier counts by gene
#'@examples
#'
#'data(ExampleData)
#' # Set up Phenotype
#' phenotype <- pheno  
#' names(phenotype) <- colnames(cnv)
#'
#' #set up datalist
#' dataSet <- list(expr,meth,cnv)
#' 
#' # set up values for expr-meth-cnv in that order
#' tailLRL <- c('left', 'right', 'left')   
#'
#'outTibLRL <- outCallTib(dataSet, phenotype=pheno, 
#'                          names=c('Expr', 'Meth', 'CNV'), tail=tailLRL)
#'
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@export

outCount <- function(data, phenotype, tail='right', corr=FALSE) {
    
    nG <- dim(data)[1]
    nS <- length(phenotype)

    outCount <- vector(length = nG)

    nData <- data[, phenotype==0]

    # generate Tibshirani eqn 2.2
    # putting genes on same scale
    medG <- apply(data, 1, median)
    madGN <- apply(nData, 1, mad)

    for (i in 1:nG) {
        temp <- (data[i, ] - medG[i]) / madGN[i]
        data[i, ] <- temp
    }

    # count outliers relative to scaled data
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
            outCount[i] <- sum(Ind)
        }
    } else if (tail == 'left') {
        for (i in 1:nG) {
            Ind <- tData[i, ] < quantG[i] - iqrG[i]
            outCount[i] <- sum(Ind)
        }
    }
    
    # if correcting, subtract number of normals
    # that are called outliers
    if (corr) {
        if (tail == 'right') {
            for (i in 1:nG) {
                Ind <- nData[i, ] > quantG[i] + iqrG[i]
                outCount[i] <- outCount[i] - sum(Ind)
            }
        } else if (tail == 'left') {
            for (i in 1:nG) {
                Ind <- nData[i, ] < quantG[i] - iqrG[i]
                outCount[i] <- outCount[i] - sum(Ind)
            }
        }
    }


    names(outCount) <- rownames(data)
    return(outCount)
}