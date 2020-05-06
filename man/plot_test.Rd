\name{plot_test}
\alias{plot_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting ribosome footprint data at bin level.
}
\description{
This function visualizes the ribosome footprint in the form of P-site frequency at per bin level along the total transcript.
}
\usage{
plot_test(result, genes.list=NULL,threshold=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
    Data object resulting from \code{RiboDiPA} function.
}  \item{genes.list}{
\code{genes.list} is the list of genes for visualization. If \code{genes.list} is not specified, then only genes with significant differential patterns specified by q-value \code{threshold} will be plotted. If \code{genes.list} is not NULL, then \code{threshold} argument will be ignored.
} \item{threshold}{
   The q-value threshold for genes whose footprint to be visualized. This argument is ignored if \code{genes.list} is not NULL.
}
}
\value{
Bin-level tracks of genes and test results.
}
