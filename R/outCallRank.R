#'outCallRank
#'
#'Counts outliers by the Ghosh method and generates list objects with all 
#'outliers noted
#'@usage outCallRank (dataSet, phenotype, thres= 0.05, tail='right', corr=FALSE,
#' offsets=NULL, names=NULL)
#'@param dataSet Set of matrices of molecular data
#'@param phenotype A vector of 0s and 1s of length nSample, where 1 = case, 
#'0 = control
#'@param thres Alpha value
#'@param tail A vector equal to the number of matrices with values left or right
#' for where to find outliers
#'@param corr Whether to correct for normal outliers
#'@param offsets A vector equal to the number of matrices which sets the minimum
#' value relative to normal to call outlier (corrected rank only)
#'@param names A vector equal to the number of matrices to name molecular type 
#'of data (e.g., CNV)
#'@return A list with all specific outlier calls for each molecular type in each
#' case sample
#'@examples
#'
#'data(ExampleData)
#' 
#' #set up dataSet
#' dataSet <- list(expr, meth,cnv)
#' 
#' # Set up Phenotype
#' phenotype <- pheno  
#' names(phenotype) <- colnames(cnv)
#' 
#' # set up values for expr-meth-cnv in that order
#' tailLRL <- c('left', 'right', 'left')   
#'
#'outRankLRL <- outCallRank(dataSet, phenotype, names=c('Expr',
#'                              'Meth', 'CNV'), tail=tailLRL)
#'
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@references D. Ghosh. (2010). Discrete Nonparametric Algorithms for Outlier 
#'Detection with Genomic Data. J. Biopharmaceutical Statistics, 20(2), 193-208.
#'@export

outCallRank <- function(dataSet, phenotype, thres= 0.05, tail='right',
								corr=FALSE, offsets=NULL, names=NULL) {
    
    ##function input checks
    if (all(thres <= 0 & thres > 1)) {
        stop("the thres alpha value must be between 0 and 1")
    }
    
    if (all(tail == "right" || tail == "left")) {
        
    }else {stop("values in 'tail' must be 'right' or 'left'")}
    
    if (length(dataSet) != length(tail)){
        stop("length of 'tail' must equal length of 'dataSet'")
    }
    
    
    if (is.null(names)) {
        names <- vector(length=length(dataSet), mode='character')
        for (d in 1:length(dataSet)) {
            names[d] <- paste('Data', d)
        }
    
    }
    temp <- dataSet[[1]]
#    temp <- temp[[1]]
    nG <- dim(temp)[1]
    outList <- list()

    if (!corr) {
        offsets <- rep(0.0, dim(temp)[2])
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
                    sum((baseline + adjust) > tumor[j]))
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
        empirP <- empirP < thres
        rownames(empirP) <- rownames(data)
        colnames(empirP) <- colnames(tData)
        outList[[d]] <- empirP
    }
    names(outList) <- names
    return(outList)
}