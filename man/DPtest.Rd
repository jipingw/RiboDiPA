\name{DPtest}
\alias{DPtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Differential pattern analysis of Ribo-seq data
}
\description{
The normalized gene data are pooled into a large matrix, where parameter estimations and tests are performed. Within each gene, multiplicity correction are then performed for codon/bin-level p-values. The minimum of adjusted codon/bin-level p-value is defined to be the gene-level p-value.
}
\usage{
DPtest(data, classlabel, method=c('gtxr','qvalue'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
A list of named matrices input from the \code{DataBinning} function. In each element of the list, rows correspond to replicates, columns correspond to bins.
}
\item{classlabel}{
For matrix input: a DataFrame or data.frame with at least a column \code{comparison}. In \code{comparison}, \code{1}s stand for the reference condition, \code{2}s stand for the target condtion, and \code{0}s represent replicates is not invloved in the test, if present. Rows of \code{classlabel} correspond to rows of \code{data}, which are biological replicates.
}
\item{method}{
For a 2-component character vector input: the first argument is the multiplicity correction method for codon/bin-level p-value adjustment. The second argument is the multiplicity correction method for gene-level p-value adjustment. Methods include: "qvalue" for q-value from \code{qvalue} pacakge, "gtxr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none" from the \code{elitism} package.
}
}
\value{
\item{bin}{A List object of codon/bin-level results. Each element of list is of a gene, containing codon/bin results columns: \code{pvalue}, \code{log2FoldChange}, and the adjusted p-value named by the first string in \code{method}.
}
\item{gene}{A DataFrame object of gene-level results. It contains columns: \code{tvalue}, \code{pvalue}, and the adjusted p-value named by the second string in \code{method}.
}
\item{small}{Names of genes without sufficient reads}
\item{classlabel}{
The same as  input \code{classlabel}.
}
\item{data}{
Subset of input \code{data}, including all genes reported in \code{bin} and \code{gene}.
}
\item{method}{
The same as input \code{method}.
}
}
\details{
Using binned data, this function first estimates normalizing constant by exclusing outlier bins which may represent the true differential pattern. An outlier bin is defined as that whose log2-fold change value is more than 1.5 interquartile ranges below the first quartile or above the third quartile. For a given gene, the normalizing constant is defined based on the total read counts from each replicate.

It then performs differential pattern testing on P-site counts bin by bin for each gene. Briefly, counts are modeled by a negative binomial distribution to call bins with statistically significant differences across conditions, bin level p-values are adjusted for multiple hypothesis testing for a given gene, and then the smallest p-value for a gene is adjusted to control for multiple hypothesis testing across all genes.

Additionally, the T-value is a supplementary statistic that quantifies the magnitude of difference between conditions, with larger numbers indicating a greater difference. The $T$-value is defined to be 1-cosine of the angle between the first right singular vectors of the footprint matrices of the two conditions under comparison. It ranges from 0-1, with larger values representing larger differences between conditions, and practically speaking, can be used to identify genes with larger magnitude of pattern difference beyond statistical significance. This might be helpful to investigators to prioritize certain genes for investigation among many that may pass the significance test for differential pattern.
}
\author{
Keren Li
}

\seealso{
\code{\link[elitism]{p.adjust}}
}
\examples{
load(system.file("extdata","sample.Rdata",package="RiboDiPA"))
classlabel <- data.frame(condition = c("mutant","mutant","wildtype","wildtype"),
    comparison=c(2,2,1,1))
rownames(classlabel) <- names.sample
result.pst <- DPtest(data=data.binned, classlabel=classlabel,
    method=c('gtxr','qvalue'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pattern similarity test}% use one of  RShowDoc("KEYWORDS")
\keyword{htest}% __ONLY ONE__ keyword per line
