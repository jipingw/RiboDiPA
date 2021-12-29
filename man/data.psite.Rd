\name{data.psite}
\alias{data.psite}
\docType{data}
\title{
An example of P-site coverage data
}
\description{
An example data set containing 4 ribo-seq replicates of 885 
genes, two biological replicates each for wild type cells and
New1 mutant cells, respectively. It is the output of P-site mapping funcion
\code{psiteMapping}. It contains the
p-site count at each location of the total transcript within each
replicate.
}
\usage{data("data.psite")}
\format{
    A list of size 4
\describe{
\item{coverage}{ribosome P-site coverage tracks}
\item{counts}{ribosome P-site total count, one count per gene}
\item{psite.mapping}{P-site mapping offset rule}
\item{exons}{relative start and end positions of each exon in the total 
transcript if a given gene, as well as genomic start and end coordinates}
}}
\source{
Raw data was adapted from Kasari et al 2019.
}
\examples{
data(data.psite)
data.binned <- dataBinning(data = data.psite$coverage, bin.width = 0, 
    zero.omit = FALSE, bin.from.5UTR = TRUE, cores = 2)
}
\keyword{datasets}
