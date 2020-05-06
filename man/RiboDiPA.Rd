\name{RiboDiPA}
\alias{RiboDiPA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main function for differential pattern analysis of Ribo-seq data
}
\description{
The normalized gene data are pooled into a large matrix, where parameter estimations and tests are performed. Within each gene, multiplicity correction are then performed for codon/bin-level p-values. The minimum of adjusted codon/bin-level p-value is defined to be the gene-level p-value.
}
\usage{
RiboDiPA(data, classlabel, method=c('gtxr','qvalue'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  For list of named matrice input from \code{DataBinning} function. In each element of the list, rows corrspond to replicates, columns corrspond to bins.
}
  \item{classlabel}{
For matrix input: a DataFrame or data.frame with at least a column \code{condition}. In \code{condition}, \code{1}'s stand for reference condition, \code{2}'s stand for target condtion, \code{0}'s replicates is not invloved in the test. Rows of \code{classlabel} correspond to rows of \code{data}.
}
\item{method}{
For a 2-component character vector input: the first argument is the multiplicity correction method for codon/bin-level p-value adjustment. The second argument is the multiplicity correction method for gene-level p-value adjustment. Methods include: "qvalue" for q-value from \code{qvalue} pacakge, "gtxr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none" from \code{elitism} package.
}
}
\value{
 \item{codon}{A List object of codon/bin-level results. Each element of list is of a gene, containing codon/bin results columns: \code{pvalue}, \code{log2FoldChange}, and the adjusted p-value named by the first string in \code{method}.
 }
 \item{gene}{A DataFrame object of gene-level results. It contains columns: \code{tvalue}, \code{pvalue}, and the adjusted p-value named by the second string in \code{method}.
 }
 \item{small}{Names of genes without sufficient reads}
 \item{classlabel}{
The same as  input \code{classlabel}.
}
\item{data}{
subset of input \code{data}, including all genes reported in \code{codon} and \code{gene}.
}
\item{method}{
The same as  input \code{method}.
}
}
\details{
The results when printed will provide information about the ribosome footrint similarity comparison.

The \code{tvalue} is defined to be 1-cos(angle between right singular vectors of comparing conditions).
}
\author{
Keren Li
}

\seealso{
\code{\link[elitism]{p.adjust}}
}
%\examples{
%classlabel <- data.frame("condition"=c(rep("eRF1",2),rep("WT",2)), "type"=names.rep,"comparison"=c(1,1,2,2))
%result <- RiboDiPA(DATA.BIN, classlabel)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pattern similarity test}% use one of  RShowDoc("KEYWORDS")
\keyword{htest}% __ONLY ONE__ keyword per line
