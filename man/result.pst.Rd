\name{result.pst}
\alias{result.pst}
\docType{data}
\title{
An example of differential pattern analysis result
}
\description{
An example output generated by the differential pattern analysis 
function  \code{diffPatternTest}, including binned data, differential 
pattern results, etc.
}
\usage{data("result.pst")}
\format{
    A list of size 5
\describe{
\item{bin}{A list object of codon/bin-level results. Each element 
of the list is the result from  a gene, containing columns: \code{pvalue},
\code{log2FoldChange}, and the adjusted p-value by method "gtxr"}       
\item{gene}{Gene-level differential pattern results, including T-value, 
p-value, and q-value}
\item{classlabel}{See \code{diffPatternTest}}
\item{data}{The input data for differential pattern analysis in the format of 
a list of named matrices. In each element of the list, rows correspond to 
replicates,columns correspond to bins.}
\item{method}{See \code{diffPatternTest}}
}}
\source{
The data was adapted from Kasari et al 2019.
}
\examples{
data(result.pst)
plotTrack(data = data.psite, genes.list = c("YDR050C", "YDR064W"), 
    replicates = NULL, exons = FALSE)
}
\keyword{datasets}
