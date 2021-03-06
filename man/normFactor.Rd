\name{normFactor}
\alias{normFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Size factors for normalization by sequencing depth
}
\description{
This function calculate the relative abundance of samples, in essence 
accounting for different sequencing depths across different Ribo-seq 
experiments.
}
\usage{
normFactor(x, condition)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A matrix of mapped P-site positions.
}
\item{condition}{
A vector of indicators. \code{1}'s stand for reference condition, 
\code{2}'s stand for target condtion, \code{0}'s represent replicates
not invloved in the abundance estimation, if present.
}
}
\examples{
data(data.binned)
x <- data.binned$YDR050C
condition <- c(2, 2, 1, 1)
normFactor(x, condition)
}
\value{
A vector of relative abundances.
}

