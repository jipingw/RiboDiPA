\name{psiteMapping}
\alias{psiteMapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{P-site mapping
}
\description{
This function computes the corresponding P-site locations of all 
RPFs and the total RPF read counts for all genes and samples.
}
\usage{
psiteMapping(bam_file_list, gtf_file, psite.mapping = "auto", 
    cores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{bam_file_list}{
A vector of bam file names to be tested. Users should include path 
names if not located in the current working directory. Index files 
(.bai) will be generated if not already present.
}
\item{gtf_file}{
Annotation file used to generate the BAM alignments. Note that a GTF 
file sourced from one organization (e.g. Ensembl) cannot be used with
BAM files aligned with a GTF file sourced from another organization (e.g. 
UCSC).
}
\item{psite.mapping}{
Rules for P-site offsets, to map a given read length of RPF to a P-site.
Input for this parameter is a string input, or a user defined matrix 
assigning the P-site mapping rules. Strings include: "center" for taking 
center of the read as the P-site, and "auto" for an optimal P-site offset,
which is the default. A user defined matrix should include two columns: 
"qwidth" and "psite", where "qwidth" is the range of possible read lengths 
and "psite" is the corresponding offset from the 5' end to map to the P-site.
}
\item{cores}{
The number of cores to use for parallel execution. If not specified, the 
number of cores is set to the value of \code{detectCores(logical = FALSE)}.
}
}
\details{
All exons from the same gene are concatenated into a total transcript in 
order to get a merged picture of translation, using the reduce() function 
from the GRanges package to accomplish the concatenation. Then, RPFs are 
mapped with respect to the total transcript and the P-site positions are 
inferred accordingly.

If `psite.mapping` is unspecified, a two-step algorithm on start codons of
CDS regions is used to compute optimal P-site offsets, following Lauria et
al (2018). First, for a given read length, the offset is calculated by taking
the distance between the first nucleotide of the start codon and the 5' most
nucleotide of the read, and then defining the offset as the 5' position 
with the most reads mapped to it. This process is repeated for all read 
lengths and then the temporary global offset is defined to be the offset 
of the read length with the maximum count. Lastly, for each read length, 
the adjusted offset is defined to be the one corresponding to the local 
maximum found in the profiles of the start codons closest to the temporary 
global offset.

The function will return a list of matrices that can then be used for 
data binning and downstream analysis, among other data objects.
}
\value{
\item{coverage}{A list object of matrices. Each element is a matrix 
representing the P-site footprints of a gene. Rows correspond to replicates
and columns corrspond to nucleotide location with respect to the total 
transcript.
}
\item{counts}{A matrix object of read counts. Rows correspond to genes 
and columns correspond to replicates.}
\item{exons}{A List object of matrices. Each element contains the relative
start and end positions of exons in the gene with respect to the total 
transcript for that given gene}
\item{psite.mapping}{
The P-site (or A-site) mapping rule used to map RPFs to P-site positions.
}
}
\references{
Lauria, F., Tebaldi, T., Bernabò, P., Groen, E., Gillingwater, T. H., 
& Viero, G. (2018). riboWaltz: Optimization of ribosome P-site positioning
in ribosome profiling data. PLoS computational biology, 14(8), e1006169.
}
\examples{
sample_names <- c("WT1.bam", "WT2.bam", "MUT1.bam", "MUT2.bam", 
    "eg.gtf")
url <- "https://github.com/mhope321/RiboDipa/raw/master/"
dest <- paste0(getwd(), "/")
for (sample_file in sample_names){
    download.file(paste0(url, sample_file), paste0(dest, 
        sample_file))
}

bam_file_list <- paste0(dest, c("WT1.bam", "WT2.bam", "MUT1.bam", "MUT2.bam"))
gtf_file <- paste0(dest, "eg.gtf")

data.psite <- psiteMapping(bam_file_list = bam_file_list, 
    gtf_file = gtf_file, psite.mapping = "auto", cores = 2)

}

\keyword{P-site}% use one of  RShowDoc("KEYWORDS")
\keyword{A-site}% __ONLY ONE__ keyword per line
