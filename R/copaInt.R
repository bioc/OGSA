#'copaInt
#'
#'Counts outliers by Tibshirani-Hastie method by calling outCount after setting 
#'up list or by rank outlier method by calling outRank
#'@usage copaInt(dataSet, phenotype, tails, thres = 0.05, method='Tibshirani', 
#' corr=FALSE, offsets=NULL)
#'@param dataSet Set of matrices of molecular data
#'@param phenotype Vector of 1 for case, 0 for control
#'@param tails Vector equal to number of matrices with values left or right for 
#'where to find outliers
#'@param thres alpha value
#'@param method Tibshirani , Rank
#'@param corr Whether to correct for normal outliers
#'@param offsets A vector equal to the number of matrices which sets the minimum
#' value relative to normal to call outlier (corrected rank only)
#'@return  A vector with outlier counts by gene
#'@examples
#'
#' data(ExampleData)
#'
#' #Set up phenotype
#' phenotype <- pheno  
#' names(phenotype) <- colnames(cnv)
#'
#' #set up values for expr-meth-cnv in that order
#' tailLRL <- c('left', 'right', 'left')   
#' 
#' #setup dataSet
#' dataSet <- list(expr, meth, cnv)
#' 
#' tibLRL <- copaInt(dataSet, phenotype, tails=tailLRL)
#' 
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@export

copaInt <- function(dataSet, phenotype, tails, thres = 0.05, 
        method='Tibshirani', corr=FALSE, offsets=NULL) {
    tempDataSet <- dataSet

    ##function input checks
    if (length(tails) != length(dataSet)) {
        stop("tails must correspond to the input dataSet")
    }
    if (all(phenotype < 0 & phenotype > 1 & phenotype%%1==0)) {
        stop("elements in phenotype must be 0 or 1")
    }
    
    nType <- length(dataSet)
    nGene <- dim(dataSet[[1]])[1]
    
    # reorder the columns in all data matrices
    # to match
    for (i in 1:nType) {
        tempMat <- dataSet[[i]]
        temp <- match(names(phenotype), colnames(tempMat))
        tempMat <- tempMat[, temp]
        dataSet[[i]] <- tempMat
    }
    
    # reorder the rows in all data matrices
    # to match
    geneNames <- rownames(dataSet[[1]])
    missingData <- rep(FALSE, length(geneNames))
    for (i in 1:nType) {
        tempMat <- dataSet[[i]]
        temp <- match(geneNames, rownames(tempMat))
        tempMat <- tempMat[temp, ]
        dataSet[[i]] <- tempMat
        missingData[is.na(tempMat[, 1])] <- TRUE
    }


    # remove rows that are not in all matrices
    # for counting, all genes must be in each matrix
    for (i in 1:nType) {
        tempMat <- dataSet[[i]]
        tempMat <- tempMat[!missingData, ]
        dataSet[[i]] <- tempMat
    }
    nGene <- dim(tempMat)[1]
    outCts <- rep(0, nGene)
    names(outCts) <- rownames(tempMat)

    # remove missing data and run Tibshirani method
    # counting outliers
    if (method == 'Tibshirani') {
        for (i in 1:nType) {
            dataMat <- dataSet[[i]]
            noData <- is.na(dataMat[1, ])
            dataMat <- dataMat[, !noData]
            phenoMat <- phenotype[!noData]
            temp <- outCount(dataMat, phenoMat, tail=tails[i], corr=corr)
            outCts <- outCts + temp
        }
    } else if (method == 'Rank') {
        for (i in 1:nType) {
            dataMat <- dataSet[[i]]
            noData <- is.na(dataMat[1, ])
            dataMat <- dataMat[, !noData]
            phenoMat <- phenotype[!noData]
            dataSet[i] <- list(exprs = dataMat)
            phenotype <- phenoMat
        }
        outCts <- outRank(dataSet, phenotype, thres = thres, tail=tails, 
                          corr=corr, offsets=offsets)
    } else {
        print('Unknown Outlier Counting Method; Returning Zeros')
    }
    return(outCts)
}