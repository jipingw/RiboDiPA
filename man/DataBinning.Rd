\name{DataBinning}
\alias{DataBinning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data binning
}
\description{
This function bins P-sites data matrix into binned matrix with fewer spots for all genes.
}
\usage{
DataBinning(data, bin.width = 0, zero.omit = F, bin.from.5UTR = T, cores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A list of P-site footprint matrices. In each element of the list, rows corrspond to replicates, columns corrspond to nucleotides.
}
  \item{bin.width}{
Binning width per bin. If specified, it is the number of codons merged per bin; if not specified, an adaptive binning width is used.
}  \item{zero.omit}{
If value is \code{TURE}, common zeros are removed.
}  \item{bin.from.5UTR}{
When the coding region length is not any integer multiple of binning width, and if value of \code{bin.from.5UTR} is \code{TRUE}, the uneven width bins will be arranged at right tail.
}  \item{cores}{
The number of cores to use for parallel execution. If not specified, the number of cores is set to the value of \code{detectCores(logical = FALSE)}.
}

}
\details{
We recommend to use an adaptive bin width following Freedman-Diaconis rule. For "zoom-in" purpose, specified \code{bin.width} per bin can be used to check the local differential pattern, though it may lead to low power at small fold change positions and potentially high computational time.
}
\value{
A list of binned P-site footprint matrices: in each matrix, rows corrspond to replicates, columns corrspond to codons/bins.
}
\author{
Keren Li
}
\seealso{
\code{\link[RiboDiPA]{PsiteMapping}}
}
%\examples{
%DATA.BIN <- DataBinning(results$coverage, zero.omit=F)
%DATA.BIN <- DataBinning(results$coverage,bin.width=5,start_codon=F)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data binning}% use one of  RShowDoc("KEYWORDS")
\keyword{bin width}% __ONLY ONE__ keyword per line
