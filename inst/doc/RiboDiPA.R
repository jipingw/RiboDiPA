## ----setup, include = FALSE---------------------------------------------------
library(RiboDiPA)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("RiboDiPA")

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Download sample files from GitHub
sample_names <- c("WT1.bam","WT2.bam","MUT1.bam","MUT2.bam","eg.gtf")
url <- "https://github.com/mhope321/RiboDipa/raw/master/"
dest <- paste0(getwd(),"/")
for (sample_file in sample_names){
    download.file(paste0(url,sample_file),paste0(dest,sample_file))
}

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Get BAM files from the RiboDiPA package
bam_file_list <- list.files(path=dest,
    pattern=".bam$", full.names=TRUE)
names.sample <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", bam_file_list)

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Get the GTF file used for alignment
gtf_file <- list.files(path=dest,pattern=".gtf$", full.names=TRUE)

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Make the class label for the experiment
classlabel <- data.frame(condition = c("mutant","mutant","wildtype","wildtype"),
    comparison=c(2,2,1,1))
rownames(classlabel) <- names.sample

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Run the RiboDiPA pipeline with default parameters
result.pip <- RiboDiPA(bam_file_list,gtf_file,classlabel)

## -----------------------------------------------------------------------------
## View the table of output from RiboDiPA
head(result.pip$gene[order(result.pip$gene$qvalue),])

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Perform individual P-site mapping procedure
data.psite <- PsiteMapping(bam_file_list=bam_file_list, gtf_file=gtf_file, 
    psite.mapping="auto",cores=NULL)

## -----------------------------------------------------------------------------
## P-site mapping offset rule generated
data.psite$psite.mapping

## ---- warning=FALSE, message=FALSE, eval=FALSE--------------------------------
#  ## Use user specified psite mapping offset rule
#  offsets <- cbind(qwidth=c(28,29,30,31,32),psite=c(18,18,18,19,19))
#  data.psite2 <- PsiteMapping(bam_file_list, gtf_file, psite.mapping=offsets,cores=NULL)

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Merge the P-site data into bins with a fixed or an adaptive width
data.binned <- DataBinning(data=data.psite$coverage, bin.width=0, zero.omit=FALSE, 
    bin.from.5UTR=TRUE, cores=NULL)

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Merge the P-site data on each codon
data.codon <- DataBinning(data=data.psite$coverage, bin.width=1, zero.omit=FALSE, 
    bin.from.5UTR=TRUE, cores=NULL)

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Merge the P-site data on each exon and perform differential pattern analysis
result.exon <- DPtest.exon(psitemap=data.psite, classlabel=classlabel, 
    method=c('gtxr','qvalue'))

## ---- warning=FALSE, message=FALSE--------------------------------------------
## Perform differential pattern analysis
result.pst <- DPtest(data=data.binned, classlabel=classlabel, 
    method=c('gtxr','qvalue'))

## ----fig2, fig.height = 12, fig.width = 9, out.width="90%", out.height="90%", fig.align = "center", warning=FALSE, message=FALSE, results='hide'----
## Plot ribosome per nucleotide tracks of specified genes.
plot_track(data=data.psite,genes.list=c("YDR050C","YDR064W"), replicates=NULL,exons=FALSE)

## ----fig3, fig.height = 12, fig.width = 9, out.width="90%", out.height="90%", fig.align = "center", warning=FALSE, message=FALSE, results='hide'----
## Plot binned ribosome tracks of siginificant genes: YDR086C and YDR210W.
## you can specify the thrshold to redefine the significant level
plot_test(result=result.pst, genes.list=NULL, threshold=0.05) 

