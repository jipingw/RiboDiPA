\name{PsiteMapping}
\alias{PsiteMapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{P-site mapping
}
\description{
This function computes the corresponding P-site locations of all RPFs and the total RPF read counts for all genes and samples.
}
\usage{
PsiteMapping(bam_file_list, gtf_file, psite.mapping="auto",
cores=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_file_list}{
  A vector of bam file names.
}  \item{gtf_file}{
Annotation file.
}\item{psite.mapping}{
For a string input or a user defined matrix assigning the P-site mapping rule or A-site mapping rule. Strings include: "center" for taking center of read as P/A-site, and "auto" for an optimal P-site offset computing.
A user defined matrix should include two columns: "qwidth" and "psite". "qwidth" is the possible range of read length, "psite" is the corresponding P/A-site offset.
}
\item{cores}{
The number of cores to use for parallel execution. If not specified, the number of cores is set to the value of \code{detectCores(logical = FALSE)}.
}
}
\details{
All exons from the same gene are concatenated into a total transcript in order to get a merged picture of translation. Then RPFs are mapped with reference to the total transcript and the P-site positions are inferred accordingly.
If `psite.mapping` is unspecified, a two-step algorithm on start codons of CDS is used to compute optimal P-site offsets based on data following Lauria et al (2018). First, the offset of each read length is defined as the distances between the first nucleotide of start codons and the nucleotide corresponding to the maximum found in the profiles of the start codons. The temporary global offset is defined to be the offset of the read length with maximum count. Next, for each read length, the adjusted offset is defined to be the one corresponding to the local maximum found in the profiles of the start codons closest to the temporary global offset.
}
\value{
\item{coverage}{A list object of matrices. Each element is a matrix representing the P-site footprints of a gene. Rows corrspond to replicates and columns corrspond to nucleotide location with reference to the total transcript.
}
\item{counts}{A matrix object of read counts. Rows corrspond to genes and columns corrspond to replicates.}
\item{exons}{A List object of matrices. Each element contains the relative start and end positions of exons in the gene with reference to the total transcript}\item{psite.mapping}{
The P-site mapping rule or A-site mapping rule used.
}
}
\author{
Keren Li,Ji-Ping Wang
}
\references{
Lauria, F., Tebaldi, T., Bernabò, P., Groen, E., Gillingwater, T. H., & Viero, G. (2018). riboWaltz: Optimization of ribosome P-site positioning in ribosome profiling data. PLoS computational biology, 14(8), e1006169.
}
%\examples{
%gtf_file <- "Saccharomyces_cerevisiae.R64-1-1.98.gtf"
%fnames <- list.files(pattern="\\.bam$")
%results <- PsiteMapping(bam_file_list=fnames,gtf_file=gtf_file)
%}

\keyword{P-site}% use one of  RShowDoc("KEYWORDS")
\keyword{A-site}% __ONLY ONE__ keyword per line
