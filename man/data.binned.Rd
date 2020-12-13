\name{data.binned}
\alias{data.binned}
\docType{data}
\title{
An example of binned P-sites data
}
\description{
A example data containing binned ribosome P-site tracks of 4 replicates on 885 
genes, two biological replicates each for wild type cells and
New1 mutant cells, respectively. It is the output of the data binning function
\code{dataBinning} on P-site coverage data, and input for \code{diffPatternTest}
function for differential pattern analysis.
}
\usage{data("data.binned")}
\format{
    A list of 885 matrices corresponding to 885 genes: in each matrix, rows 
correspond to replicates, columns correspond to bins.
}
\source{
The raw data was  adapted from Kasari et al 2019.
}
\examples{
data(data.binned)
classlabel <- data.frame(condition = c("mutant", "mutant", 
    "wildtype", "wildtype"), comparison = c(2, 2, 1, 1))
rownames(classlabel) <- c("mutant1", "mutant2", "wildtype1", "wildtype2")
result.pst <- diffPatternTest(data = data.binned, 
    classlabel = classlabel, method = c('gtxr', 'qvalue'))
}
\keyword{datasets}
