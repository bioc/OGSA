#'outCallTibE
#'
#'Counts outliers by the Tibshirani and Hastie method and generates a list 
#'object with all outliers noted
#'@usage outCallTibE (expressionSet, tail='right', corr=FALSE, names=NULL)
#'@param expressionSet ExpressionSet object containing sets of data and 
#'phenotype information
#'@param tail Vector equal to number of matrices with values 'left' or 'right' 
#'for where to find outliers
#'@param corr whether to correct for normal outliers ONLY for compatibility, 
#'since method does not allow determining specific changes in cases, it will 
#'just print message if corr = TRUE
#'@param names Vector equal to number of matrices to name molecular type of data
#' (e.g., 'CNV').
#'@import Biobase
#'@return A list with all specific outlier calls for each molecular type in each
#' case sample
#'@examples
#'
#'data(ExampleData)
#'data('KEGG_BC_GS')
#'
#'  library(Biobase)
#' # building the Annotated Data Frame
#'  phenoData <- AnnotatedDataFrame(
#'      data.frame(
#'         type = factor(x = pheno, labels = c("Control", "Case")),
#'          row.names = colnames(expr)
#'      )
#'  )
#' # build environment
#'  inputData <- list2env(list(exprs = expr, meth = meth, cnv = cnv))
#'
#' # build expressionSet - other information can be added here
#'  expressionSet <- ExpressionSet(inputData, phenoData)   
#'  
#' # set up values for for the tails in the order that they are exported, for example:
#' tailLRL <- c('left', 'right', 'left') 
#'  
#'  
#' outTibLRL <- outCallTibE(expressionSet, names=c('Expr', 'Meth', 'CNV'), tail=tailLRL)
#' 
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@references D. Ghosh. (2010). Discrete Nonparametric Algorithms for Outlier 
#'Detection with Genomic Data. J. Biopharmaceutical Statistics, 20(2), 193-208.
#'@export

outCallTibE <- function(expressionSet, tail='right', corr=FALSE, names=NULL) {
    
    dataSet <- expressionSetDataSet(expressionSet)    
    phenotype <- expressionSetPheno(expressionSet)
    
    ##function input checks
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
    
    if (corr) {
        print("No Correction Applicable in Tibshirani Indicators")
    }
    temp <- dataSet[[1]]
    nG <- dim(temp)[1]
    nT <- dim(temp)[2]
    outList <- list()
    
    for (d in 1:length(dataSet)) {
        data <- dataSet[[d]]

        nS <- length(phenotype)
        thisTail <- tail[d]
        
        nData <- data[, phenotype==0]
        
        # generate Tibshirani eqn 2.2
        # putting genes on same scale
        medG <- apply(data, 1, median)
        madGN <- apply(nData, 1, mad)
        
        for (i in 1:nG) {
            temp <- (data[i, ] - medG[i]) / madGN[i]
            data[i, ] <- temp
        }
        tData <- data[, phenotype==1]
        nData <- data[, phenotype==0]
        nT <- dim(tData)[2]
        
        # count outliers relative to scaled data
        iqrG <- apply(data, 1, IQR)
        if (thisTail == 'right') {
            quantG <- apply(data, 1, quantile, probs = 0.75)
        } else if (thisTail == 'left') {
            quantG <- apply(data, 1, quantile, probs = 0.25)
        }
        
        # generate indicators and place in output matrix
        empirP <- matrix(nrow=nG, ncol=nT)
        if (thisTail == 'right') {
            for (i in 1:nG) {
                Ind <- tData[i, ] > (quantG[i] + iqrG[i])
                empirP[i, ] <- Ind
            }
        } else if (thisTail == 'left') {
            for (i in 1:nG) {
                Ind <- tData[i, ] < (quantG[i] - iqrG[i])
                empirP[i, ] <- Ind
            }
        }
        rownames(empirP) <- rownames(data)
        colnames(empirP) <- colnames(tData)
        outList[[d]] <- empirP
    }
    
    names(outList) <- names
    return(outList)
}