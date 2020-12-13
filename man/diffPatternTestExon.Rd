\name{diffPatternTestExon}
\alias{diffPatternTestExon}
\title{
Main function for differential pattern analysis of exon-binned Ribo-seq data}
\description{
An alternative version of \code{diffPatternTest} for exon level binning. Both
data binning and differential pattern analysis are implemented. Instead of a 
fixed width or adaptive method, the positions of exons in the genome are 
used as bins. Therefore the number of exons per gene and their relative 
sizes determines the bins used for differential pattern testing.}
\usage{
diffPatternTestExon(psitemap, classlabel, 
    method = c("gtxr", "qvalue"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{psitemap}{
A list object from value of \code{psiteMapping} function. In \code{psitemap},
list elements \code{coverage} and \code{exons} are required.
}
\item{classlabel}{
For matrix input: a DataFrame or data.frame with at least a column 
\code{comparison}. In \code{comparison}, \code{1}s stand for the 
reference condition, \code{2}s stand for target condtion, \code{0}s 
represent replicates not invloved in the test, if present. Rows of 
\code{classlabel} correspond to rows of \code{data}.
}
\item{method}{
For a 2-component character vector input: the first argument is the 
multiplicity correction method for exon-level p-value adjustment. The 
second argument is the multiplicity correction method for gene-level 
p-value adjustment. Methods include: "qvalue" for q-value from 
\code{qvalue} pacakge, "gtxr", "holm", "hochberg", "hommel", "bonferroni",
"BH", "BY","fdr", "none" from \code{elitism} package.}
}

\details{
For mammalian species, when the reads are sparse, it’s more meaningful 
to perform a exon level pattern analysis. \code{diffPatternTestExon}() provides
the option of exon level pattern differentiation analysis by treating 
each exon as one bin. But for organisms such as yeast, as most genes 
only contain one exon, the exon-level analysis is not meaningful 
since the analysis will simply result in the RNA-seq type of analysis,
i.e. differential abundance test instead of the pattern analysis. 
Using \code{diffPatternTestExon}() on yeast data is not for organisms with 
minimal alternative splicing or multiple exons. For a given gene, 
the normalizing constant is estimated at codon level.
}
\value{
\item{bin}{A List object of exon-level results. Each element of list
is of a gene, containing exon results columns: \code{pvalue}, 
\code{log2FoldChange}, and the adjusted p-value named by the first 
string in \code{method}.
}
\item{gene}{A DataFrame object of gene-level results. It contains 
columns: \code{tvalue}, \code{pvalue}, and the adjusted p-value 
named by the second string in \code{method}.
}
\item{small}{Names of genes without sufficient reads}
\item{classlabel}{
The same as  input \code{classlabel}.
}
\item{data}{
A list of exon-binned P-site footprint matrices: in each matrix, 
rows corrspond to replicates, columns corrspond to exons. All genes
reported in \code{bin} and \code{gene} are included.
}
\item{method}{
The same as  input \code{method}.
}
}


\examples{
data(data.psite)
classlabel <- data.frame(condition = c("mutant", "mutant", 
    "wildtype", "wildtype"), comparison=c(2, 2, 1, 1))
rownames(classlabel) <- c("mutant1", "mutant2", "wildtype1", "wildtype2")
result.exon <- diffPatternTestExon(psitemap = data.psite, 
    classlabel = classlabel, method = c('gtxr', 'qvalue'))
}
\seealso{
\code{diffPatternTest}
}

