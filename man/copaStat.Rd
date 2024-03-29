% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/copaStat.R
\name{copaStat}
\alias{copaStat}
\title{copaStat}
\usage{
copaStat (data, phenotype, tail='right', perms=100, permType='array')
}
\arguments{
\item{data}{A matrix of nGene by nSample}

\item{phenotype}{A vector of 0s and 1s of length nSample, where 1 = case,
0 = control}

\item{tail}{Indicates whether outliers are up (right) or down (left) outliers}

\item{perms}{The number of permutations}

\item{permType}{By all on array or by gene, if by gene increase perms
significantly and plan on lots of time; in theory array should be fine as
genes are rescaled}
}
\value{
A vector with outlier counts by gene
}
\description{
Calculates outlier statistics by the Tibshirani-Hastie method
}
\examples{
data(ExampleData)

#Set up phenotype
phenotype <- pheno
names(phenotype) <- colnames(cnv)

#set up values for expr-meth-cnv in that order
tailLRL <- c('left', 'right', 'left')

#setup dataList
dataSet <- list(expr, meth, cnv)

data <- dataSet[[1]]

tibL <- copaStat(data, phenotype, tail='right', perms=100, permType='array')
}
\references{
Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
& Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated
Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
 Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
}

