\name{psiteCal}
\alias{psiteCal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P-site position of reads
}
\description{
This function outputs the P-site position, provided the CIGAR 
string of the alignment and the start position of read. It is 
implemented as a C++ function using the Rcpp package.
}
\usage{
psiteCal(cigar, start, psitemap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{cigar}{
A vector of CIGAR strings.
}
\item{start}{
A vector of read start positions.
}
\item{psitemap}{
A vector of relative P-site positions, which describe the offset 
from the 5 prime most nucleotide of the read to the P-site.
}
}
\value{
A vector of P-site positions for the reads.
}
\examples{
ex.cigar <- c("21M74731N7M", "2S11M57302N12M3S", "28M", "27M1S")
ex.start <- c(177640, 249163, 249286, 249290)
ex.psitemap <- c(9, 18, 9, 9)
psiteCal(ex.cigar, ex.start, ex.psitemap)
}


