#'expressionSetDataSet
#'
#'Returns list of data matrices
#'@usage expressionSetDataSet(expressionSet)
#'@param expressionSet ExpressionSet object containing sets of data and 
#'phenotype information
#'@keywords internal
#'@return A list with all the molecular information
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
#'# retrieve dataSet
#'  dataSet <- expressionSetDataSet(expressionSet)
#'  
#'@references Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
#' & Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated 
#' Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
#'  Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
#'@export

expressionSetDataSet <- function(expressionSet) {
    
    if (!(class(expressionSet)[1] == "ExpressionSet")) {
        stop("The object must be an ExpressionSet")
    }
    
    matrixNames <- ls(expressionSet@assayData)
    
    matrixList <- list()
    dataList <- list()
    
    pheno.v <- expressionSetPheno(expressionSet)
    
    for (i in 1:length(matrixNames) ) {
        temp <- get(matrixNames[i], envir = expressionSet@assayData)
        assign(matrixNames[i], temp)
        matrixList[[i]] <- eval(parse(text = matrixNames[i]))
    }
    rm(temp,i)
    
    return(matrixList)
}