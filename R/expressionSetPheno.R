#'expressionSetPheno
#'
#'Returns pheno info
#'@usage expressionSetPheno(expressionSet)
#'@param expressionSet ExpressionSet object containing sets of data and 
#'phenotype information
#'@keywords internal
#'@return A list with all specific outlier calls for each molecular type in each
#' case sample
#'@examples
#'
#' data(ExampleData)
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
#'# retrieve phenotype data
#'  phenotype <- expressionSetPheno(expressionSet) 
#' 
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@export


expressionSetPheno <- function(expressionSet) {
    if (!(class(expressionSet)[1] == "ExpressionSet")) {
        stop("The object must be an ExpressionSet")
    }
    
    #This puts the phenotype data into the form needed
    pheno.v <- relevel(expressionSet@phenoData@data$type, "Case")
    pheno.v <- as.numeric(pheno.v)
    pheno.v <- ifelse(pheno.v == 2, 0, 1)
    
    #Get colnames for pheno.v
    matrixNames <- ls(expressionSet@assayData)
    temp.v <- get(matrixNames[1], envir = expressionSet@assayData)
    names(pheno.v) <- colnames(temp.v)
    
    rm(temp.v)
    return(pheno.v)
}