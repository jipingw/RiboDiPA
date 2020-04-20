\name{psitecal}
\alias{psitecal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P-site position of reads
}
\description{
This function outputs the p-site position providing the cigar string and the start position of read.
}
\usage{
psitecal(cigar, start, psitemap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cigar}{
A vector of cigar strings.}
  \item{start}{
A vector of read start positions.
}
  \item{psitemap}{
A vector of relative P-site position (position of P-site - position of first left position of read).
}
}
\value{
A vector of P-site positions.
}



