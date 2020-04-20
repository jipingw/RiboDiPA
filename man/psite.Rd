\name{PsiteMapping}
\alias{PsiteMapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{P-site mapping
}
\description{
This function computes P-sites of reads and read counts for all genes and samples.
}
\usage{
PsiteMapping(bam_file_list, gtf_file, cores=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_file_list}{
  A vector of bam file names.
}  \item{gtf_file}{
Annotation file.
}\item{cores}{
The number of cores to use for parallel execution. If not specified, the number of cores is set to the value of \code{detectCores(logical = FALSE)}.
}
}
\details{
Since true RPFs only originate from ribosomes that are actively translating in coding regions, in cases where a gene has multiple splicing isoforms, it is difficult to distinguish which isoform an observed Ribo-seq read comes from. All exons from the same gene are concatenated into a total transcript in order to get a merged picture of translation. Then RPFs are mapped to a given gene into P-site positions.

User can assign the P-site mapping rule or A-site mapping rule by redefine P/A-site rule function \code{.psite}, as needed. Input of \code{.psite} should be length of reads and output corresponding relative P/A-site position (position of P/A-site - position of first left position of read).
}
\value{
\item{coverage}{A List object of matrices. Each element is the P-site footprint of a gene. Rows corrspond to replicates, columns corrspond to nucleotides.
}
\item{counts}{A Matrix object of read counts. Rows corrspond to genes, columns corrspond to replicates.}
}
\author{
Jiping Wang
}
%\examples{
%gtf_file <- "Saccharomyces_cerevisiae.R64-1-1.98.gtf"
%fnames <- list.files(pattern="\\.bam$")
%results <- PsiteMapping(bam_file_list=fnames,gtf_file=gtf_file)
%}

\keyword{P-site}% use one of  RShowDoc("KEYWORDS")
\keyword{A-site}% __ONLY ONE__ keyword per line
