\name{SizeFactors}
\alias{SizeFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relative abundance
}
\description{
This function calculate the relative abundance of samples.
}
\usage{
SizeFactors(x, condition)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix of footprint.
}
  \item{condition}{
A vector of indicators. \code{1}'s stand for reference condition, \code{2}'s stand for target condtion, \code{0}'s replicates is not invloved in the abundance estimation.
}
}

\value{
A vector of relative abundance.
}

