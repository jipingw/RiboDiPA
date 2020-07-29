## ----setup, include = FALSE---------------------------------------------------
library(RiboDiPA)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
## specify bam_files from Ribo-seq, you should replace it by your own bam files
bam_file_list <- list.files(path=system.file("extdata",package="RiboDiPA"),
                        pattern=".bam$", full.names=TRUE)
names.sample <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", bam_file_list)

## -----------------------------------------------------------------------------
## gtf_file you used for Ribo-seq alignment, replace it by corresponding gtf file
gtf_file <- list.files(path=system.file("extdata",package="RiboDiPA"),
                    pattern=".gtf$", full.names=TRUE)

## ---- warning=FALSE, message=FALSE, results='hide', eval=FALSE----------------
#  ## convert RPF reads into the P-site position on merged exons
#  psite.mapping <- cbind(qwidth=20:32,psite=c(12,13,13,13,13,13,13,12,12,13,13,13,13))
#  data.psite <- PsiteMapping(bam_file_list=bam_file_list, gtf_file=gtf_file, psite.mapping=psite.mapping, cores=NULL)

## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  ## merge the P-site data into bins with a fixed or an adaptive width
#  data.binned <- DataBinning(data=data.psite$coverage, bin.width=1, zero.omit=F, bin.from.5UTR=T, cores=NULL)

## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  ## perform differential pattern analysis
#  classlabel <- data.frame("condition"=c(rep("A",2),rep("B",2)),
#                           "type"=names.sample,"comparison"=c(1,1,2,2))
#  rownames(classlabel) <- names.sample
#  result.pst <- RiboDiPA(data=data.binned, classlabel=classlabel,
#                     method=c('gtxr','qvalue'))
#  result.pst$gene

## ---- echo = FALSE------------------------------------------------------------
load(system.file("extdata","1.Rdata",package="RiboDiPA"))
result.pst$gene

## ----fig2, fig.height = 12, fig.width = 9, out.width="90%", out.height="90%", fig.align = "center", warning=FALSE, message=FALSE, results='hide'----
## plot ribosome per nucleotide tracks of specified genes.
plot_track(data=data.psite, genes.list=c("YDR086C","YDR210W"), replicates=NULL,exons=FALSE)

## ----fig3, fig.height = 12, fig.width = 9, out.width="90%", out.height="90%", fig.align = "center", warning=FALSE, message=FALSE, results='hide'----
## plot binned ribosome tracks of siginificant genes: YDR086C and YDR210W.
## you can specify the thrshold to redefine the significant level
plot_test(result=result.pst, genes.list=c("YDR086C","YDR210W"), threshold=0.05) 

