#'testGScogps
#'
#'Performs gene set test on outlier counts
#'@usage testGScogps (outlierCts, geneSets)
#'@param outlierCts Vector with gene names and outlier counts
#'@param geneSets List of gene sets
#'@import limma
#'@return A vector with rank sum gene set statistics
#'@examples
#'\dontrun{
#'data(ExampleData)
#'data('_BC_GS')
#'
#' #Set up your phenotype
#' phenotype <- rep(0, 69)  
#' phenotype[annot[, 3] == 'Event'] <- 1
#' names(phenotype) <- rownames(annot)
#' 
#' # set up values for expr-meth-cnv in that order
#' tailLRL <- c('left', 'right', 'left')   
#' 
#' dataSet <- list(expr, meth, cnv)
#' 
#' tibLRLcorr <- copaInt(dataSet, phenotype, tails=tailLRL, corr=TRUE)
#' gsTibLRLcorr <- testGScogps(tibLRLcorr, pathGS)
#' }
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@export

testGScogps <- function(outlierCts, geneSets) {
        
    N <- length(geneSets)
    pGS <- vector(length=N)

    duplGenes <- duplicated(names(outlierCts))
    outlierCts <- outlierCts[!duplGenes]

    for (i in 1:N) {
        genes <- geneSets[[i]]
        genes <- unique(genes)
        genes <- names(outlierCts) %in% genes
        if (sum(genes) > 4) {
            pGS[i] <- wilcoxGST(genes, outlierCts, alternative='mixed')
        } else {
            pGS[i] <- 1.0
        }
    }
    names(pGS) <- names(geneSets)
    return(pGS)
}