\name{plot_test}
\alias{plot_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting ribosome per bin footprints
}
\description{
This function  visualizes the binned tracks of genes with test information annotated
}
\usage{
plot_test(result, genes.list=NULL,threshold=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
Test result from \code{RiboDiPA} function
}  \item{genes.list}{
\code{genes.list} is the subset of investigated genes requested visualization. If \code{genes.list} is not specified, then only genes with significant differential pattern will be output.
} \item{threshold}{
The nomial q-value for significance.
}
}
\value{
Bin-level tracks of genes and test results.
}
