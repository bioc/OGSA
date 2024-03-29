% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/outCount.R
\name{outCount}
\alias{outCount}
\title{outCount}
\usage{
outCount (data, phenotype, tail='right', corr=FALSE)
}
\arguments{
\item{data}{A matrix of nGene by nSample}

\item{phenotype}{A vector of 0s and 1s of length nSample, where 1 = case,
0 = control}

\item{tail}{Indicates whether outliers are up (right) or down (left) outliers}

\item{corr}{Whether to correct for normal outliers}
}
\value{
A vector with outlier counts by gene
}
\description{
Counts outliers by the Tibshirani and Hastie method. Adds the ability to
subtract for outliers in the normals using corr = TRUE
}
\examples{
data(ExampleData)
# Set up Phenotype
phenotype <- pheno
names(phenotype) <- colnames(cnv)

#set up datalist
dataSet <- list(expr,meth,cnv)

# set up values for expr-meth-cnv in that order
tailLRL <- c('left', 'right', 'left')

outTibLRL <- outCallTib(dataSet, phenotype=pheno,
                         names=c('Expr', 'Meth', 'CNV'), tail=tailLRL)
}
\references{
Ochs, M. F., Farrar, J. E., Considine, M., Wei, Y., Meshinchi, S.,
& Arceci, R. J. (n.d.). Outlier Analysis and Top Scoring Pair for Integrated
Data Analysis and Biomarker Discovery. IEEE/ACM Transactions on Computational
 Biology and Bioinformatics, 1-1. doi:10.1109/tcbb.2013.153
}

