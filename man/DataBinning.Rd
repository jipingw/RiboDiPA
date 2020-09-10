\name{DataBinning}
\alias{DataBinning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data binning
}
\description{
This function bins a mapped P-site data matrix for a given gene into a binned matrix, for statistical testing downstream. Data can be adaptively binned, where each gene has a different number of bins and bin widths, but the bin positions for a given gene are the same across different conditions and replicates. Alternatively, data can also be binned into bins of fixed width, down to the single-codon level.
}
\usage{
DataBinning(data, bin.width = 0, zero.omit = FALSE, bin.from.5UTR = TRUE, cores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
A list of mapped P-site position matrices from the \code{coverage} object of the \code{PsiteMapping} function. In each element of the list, rows correspond to replicates, while columns correspond to nucleotides across the total transcript.
}
\item{bin.width}{
Binning width per bin. If specified, it is the number of codons merged per bin; if not specified, an adaptive binning width method is used.
}
\item{zero.omit}{
If the \code{zero.omit} argument is set to \code{TRUE}, bins with zero mapped P-site counts across all replicates are removed from the differential pattern analysis.
}
\item{bin.from.5UTR}{
When the coding region length is not any integer multiple of binning width, and if value of \code{bin.from.5UTR} is set to \code{TRUE}, the uneven width bins will be arranged at the 3' end of the total transcript. If set to \code{FALSE}, binning will proceed from the 3' end.
}
\item{cores}{
The number of cores to use for parallel execution. If not specified, the number of cores is set to the value of \code{detectCores(logical = FALSE)}.
}

}
\details{
We recommend to use an adaptive bin width \eqn{h} following the Freedman-Diaconis rule, \deqn{h= \frac{2  IQR}{m^{1/3}}}. To see certain regions of transcripts in greater detail (e.g. near the start and stop codons), a specified \code{bin.width} per bin can be used to check the local differential pattern, though it may lead to low power at small fold change positions and potentially high computational time.
}
\value{
A list of binned P-site footprint matrices: in each matrix, rows correspond to replicates, columns correspond to bins.
}
\author{
Keren Li
}
\seealso{
\code{\link[RiboDiPA]{PsiteMapping}}
}
\examples{
load(system.file("extdata","sample.Rdata",package="RiboDiPA"))
data.binned <- DataBinning(data=data.psite$coverage, bin.width=0, zero.omit=FALSE,
    bin.from.5UTR=TRUE, cores=NULL)
data.codon <- DataBinning(data=data.psite$coverage, bin.width=1, zero.omit=FALSE,
    bin.from.5UTR=TRUE, cores=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data binning}% use one of  RShowDoc("KEYWORDS")
\keyword{bin width}% __ONLY ONE__ keyword per line
