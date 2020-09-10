\name{plot_track}
\alias{plot_track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting ribosome footprint data at the mapped P-site level
}
\description{
This function visualizes the ribosome footprint in the form of P-site frequency at the per nucleotide level along the total transcript. Plotting is implemented with the ggplot2 package.
}
\usage{
plot_track(data, genes.list, replicates=NULL,exons=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
Data object from \code{PsiteMapping} function or wrapper function \code{RiboDiPA}.
}
\item{genes.list}{
A list of genes for visualization.
}
\item{replicates}{
Names of the replicates for which the footprint to visualize. The default is for all.
}\item{exons}{
If value is \code{TRUE}, Ribo-seq footprints per exon of specified genes are also output.
}
}
\value{
Visualizes the Ribo-seq per nucleotide footprint on merged exons of the genes and replicates specified. If \code{exons} is \code{TRUE}, Ribo-seq footprint per exon of specified genes is also output.
}
\examples{
load(system.file("extdata","sample.Rdata",package="RiboDiPA"))
plot_track(data=data.psite,genes.list=c("YDR050C","YDR064W"), replicates=NULL,exons=FALSE)
}
