#'outRank
#'
#'Counts outliers by the Ghosh method.
#'@usage outRank (dataSet, phenotype, thres= 0.05, tail='right', corr=FALSE,
#' offsets=NULL)
#'@param dataSet Set of matrices of molecular data
#'@param phenotype Vector of 1 for case, 0 for control
#'@param thres Alpha value
#'@param tail Vector equal to number of matrices with values 'left' or 'right' 
#'for where to find outliers
#'@param corr Whether to correct for normal outliers
#'@param offsets Vector equal to number of matrices which sets minimum value 
#'relative to normal to call outlier (corrected rank only)
#'@return A vector with outlier counts by gene
#'@examples
#'
#'data(ExampleData)
#'
#' # Set up Phenotype
#' phenotype <- pheno  
#' names(phenotype) <- colnames(cnv)
#' 
#' #set up dataSet
#' dataSet <- list(expr, meth, cnv)
#' 
#' # set up values for expr-meth-cnv in that order
#' tailLRL <- c('left', 'right', 'left')   
#'
#' outRankLRL <- outRank(dataSet, phenotype, thres= 0.05, tail=tailLRL, 
#'                              corr=FALSE, offsets=NULL)
#'
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@references D. Ghosh. (2010). Discrete Nonparametric Algorithms for Outlier 
#'Detection with Genomic Data. J. Biopharmaceutical Statistics, 20(2), 193-208.
#'@export

outRank <- function(dataSet, phenotype, thres= 0.05, tail='right', corr=FALSE, 
												offsets=NULL) {
    print(dim(dataSet[[1]]))
    temp <- dataSet[[1]]
    print(dim(temp))
#    temp <- temp[[1]]
    nG <- dim(temp)[1]
    print(nG)

    if (!corr) {
        offsets <- rep(0.0,length(phenotype))
    } else if (is.null(offsets)) {
        print('No Offsets Set with Correction Requested')
        return()
    }

    outP <- matrix(nrow=nG, ncol=2)
    outCount <- rep(0, nG)

    for (d in 1:length(dataSet)) {
        data <- dataSet[[d]]
#        phenotype <- data[[2]]
#        data <- data[[1]]
        nS <- length(phenotype)
        thisTail <- tail[d]
        
        adjust <- offsets[d]
        nData <- data[, phenotype==0]
        tData <- data[, phenotype==1]
        nT <- dim(tData)[2]
        
        # generate empicical pValues as the
        # number of sum(normals{<, >}tumor)/nN
        empirP <- matrix(nrow=nG, ncol=nT)
        if (thisTail == 'right') {
            for (i in 1:nG) {
                tumor <- tData[i, ]
                baseline <- nData[i, ]
                result <- sapply(1:length(tumor), function(j)
                    sum((baseline+adjust) > tumor[j]))
                empirP[i, ] <- result / length(baseline)
            }
        } else if (thisTail == 'left') {
            for (i in 1:nG) {
                tumor <- tData[i, ]
                baseline <- nData[i, ]
                result <- sapply(1:length(tumor), function(j)
                    sum((baseline-adjust) < tumor[j]))
                empirP[i, ] <- result / length(baseline)
            }
        }
        outP <- cbind(outP, empirP)
    }
    yP <- dim(outP)[2]
    outP <- outP[, 3:yP]
    
    # count the number of genes that have outliers based
    # on desired threshold
    for (i in 1:nG) {
        outCount[i] <- sum(outP[i, ] < thres)
    }

    names(outCount) <- rownames(data)
    return(outCount)
}