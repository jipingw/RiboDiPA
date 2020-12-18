\name{RiboDiPA}
\alias{RiboDiPA}
\title{
A wrapper function for the RiboDiPA pipeline}
\description{
A wrapper function for the RiboDiPA pipeline, that will call PsiteMapping, 
DataBinning, and DPTest in order. This function is provided for users' 
convenience and requires BAM files (one per biological replicate), a 
GTF file, and a classlabel object describing what comparisons to make.
The minimal output from the function is a list of genes with signficant 
differential patterns.
}
\usage{
RiboDiPA(bam_file_list, gtf_file, classlabel, psite.mapping = "auto",
    exon.binning = FALSE, bin.width = 0,zero.omit = FALSE, 
    bin.from.5UTR = TRUE, method = c("gtxr", "qvalue"), cores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{bam_file_list}{
A vector of bam file names to be tested. Users should include path names 
if not located in the current working directory. Index files (.bai) will
be generated if not already present.
}
\item{gtf_file}{
Annotation file used to generate the BAM alignments. Note that a GTF 
file sourced from one organization (e.g. Ensembl) cannot be used with 
BAM files aligned with a GTF file sourced from another organization 
(e.g. UCSC).
}
\item{classlabel}{
For matrix input: a DataFrame or data.frame with at least one column
named \code{comparison}. In \code{comparison}, \code{1} stands for the 
reference condition, \code{2} stands for the treatment condtion, and 
\code{0} represents replicates not invloved in the test. Rows of 
\code{classlabel} correspond to the data, which is one row per BAM file.
}
\item{psite.mapping}{
Rules for P-site offsets, to map a given read length of RPF to a P-site.
See \code{psiteMapping} for details.
}
\item{exon.binning}{
Logical indicator. If \code{exon.binning} is \code{TRUE}, use the exon 
boundaries indicated in the GTF file as bins for testing, otherwise, 
adaptive or fixed binning will be performed.
}
\item{bin.width}{
Binning width per bin. \code{0} represents adaptive binning, which is 
the default method. The minimal value for fixed-width binning is \code{1},
which represent single-codon binning. See \code{dataBinning} for details.
}
\item{zero.omit}{
If this parameter is \code{TRUE}, bins with zero reads across all 
replicates for a given gene are removed.
}
\item{bin.from.5UTR}{
When the coding region length is not any integer multiple of binning 
width, and if value of \code{bin.from.5UTR} is \code{TRUE}, the uneven 
width bins will be arranged at the 3' end of the total transcript.
}
\item{method}{
2-component character vector specifies the multiplicity correction 
method for codon/bin-level p-value adjustment. The default See 
\code{diffPatternTest} for details.
}
\item{cores}{
The number of cores to use for parallel execution. If not specified, 
the number of cores is set to the value of \code{detectCores(logical = 
FALSE)}.
}
}
\value{
\item{bin}{A List object of codon/bin-level results. Each element of 
list is of a gene, containing codon/bin results columns: \code{pvalue},
\code{log2FoldChange}, and the adjusted p-value named by the first string 
in \code{method}.
}
\item{gene}{A DataFrame object of gene-level results. It contains columns:
\code{tvalue}, \code{pvalue}, and the adjusted p-value named by the second 
string in \code{method}.
}
\item{small}{Names of genes without sufficient reads}
\item{classlabel}{
The same as  input \code{classlabel}.
}
\item{data}{
Tracks of binned data of all genes reported in \code{bin} and \code{gene}.
}
\item{method}{
The same as input \code{method}.
}
\item{coverage}{A list object of matrices. Each element is a matrix 
representing the P-site footprints of a gene. Rows corrspond to 
replicates and columns corrspond to nucleotide location with reference 
to the total transcript.
}
\item{counts}{A matrix object of read counts. Rows corrspond to genes 
and columns corrspond to replicates.}
\item{exons}{A List object of matrices. Each element contains the 
relative start and end positions of exons in the gene with reference 
to the total transcript}\item{psite.mapping}{
The P-site mapping rule or A-site mapping rule used.
}
}

\examples{
\donttest{
library(BiocFileCache)
file_names <- c("WT1.bam", "WT2.bam", "MUT1.bam", "MUT2.bam", "eg.gtf")
url <- "https://github.com/jipingw/RiboDiPA-data/raw/master/"
bfc <- BiocFileCache()
bam_path <- bfcrpath(bfc,paste0(url,file_names))

classlabel <- data.frame(
    condition = c("mutant", "mutant", "wildtype", "wildtype"),
    comparison = c(2, 2, 1, 1)
)
rownames(classlabel) <- c("mutant1","mutant2","wildtype1","wildtype2") 
result.pip <- RiboDiPA(bam_path[1:4], bam_path[5], classlabel, cores=2)
}
}
\seealso{
\code{\link[RiboDiPA]{psiteMapping}, \link[RiboDiPA]{dataBinning}, 
\link[RiboDiPA]{diffPatternTest}, \link[RiboDiPA]{diffPatternTestExon}}}

