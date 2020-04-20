\name{plot_track}
\alias{plot_track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting ribosome per nucleotide footprints
}
\description{
This function  visualizes the ribosome per nucleotide footprints with merged exons annotated.
}
\usage{
plot_track(data, genes.list, replicates=NULL, gtf_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Test result from \code{PsiteMapping} function.
}  \item{genes.list}{
Genes required visualization.
} \item{replicates}{
Replicates required visualization.
}\item{gtf_file}{
Annotation file.
}
}
\value{
Nucleotide footprint of genes.
}
