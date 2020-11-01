RiboDiPA R package
================

**Maintainer**: Ji-Ping Wang, \<<jzwang@northwestern.edu>\>

**Cite RiboDiPA package**:

L,i K., Hope, C.M., Wang, X.A., Wang, J.-P. (2020) "RiboDiPA: A novel tool for differential pattern analysis in Ribo-seq data." *Nucleic Acid Research* (in revision).


## What is RiboDiPA?

Ribosome profiling (also known as Ribo-seq) is a next-generation sequencing technique to investigate the translational activities of ribosomes across a wide variety of contexts.  Ribo-seq data not only provide the abundance of ribosomes bound to transcripts in the form of counts of ribosome protected fragments (RPFs), but also positional information across transcripts that could be indicative of differences in translational regulation.

**RiboDiPA**, short for **Ribo**some **Di**ferential **P**attern **A**nalysis, is a bioinformatics pipeline developed for analysis of the pattern of Ribo-seq footprint data. RiboDiPA is released as an R package to support statistical inference of translational differences between conditions. Briefly, this involves mapping Ribo-seq data to P-site counts along a total transcript of a gene, followed by binning these counts and performing bin-wise statistical testing.

## Availble format to run RiboDiPA pipeline

RiboDiPA is an R package that utilizes parallel computing functionality with some core functions written in C++, released as part of the Bioconductor suite of tools.

## Installation of RiboDiPA

```r
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RiboDiPA")
```

## Input files required to run RiboDiPA 

1) Ribo-seq alignment files (BAM), one per sample.
2) Gene Transfer Format (GTF) file for the reference genome of interest.

## RiboDiPA main features

The RiboDiPA R package executes four major functions to perform differential pattern analysis of Ribo-seq data, with optional visualization of results. An overview of the process can be seen in Figure 1: 

![A diagram of RiboDiPA.](https://github.com/jipingw/RiboDiPA-data/raw/master/Figure1.png) 

A) **GTF file parsing and exon merging**: For a given gene, all exons annotated in the GTF file are merged into a total transcript. This provides a global picture of changes across conditions for a gene, as the total transcript will capture changes in ribosome occupancy even when transcript isoform usage might change across conditions.

B) **BAM file processing and P-site mapping**: The Ribo-seq alignment files (.bam) are processed to calculate the P-site position for each RPF, with adaptable rules that users’ can specify to call P-sites from the data. The mapped P-site frequency at each nucleotide position along the total transcript is compiled for each gene of each sample.

C) **Data binning**: To overcome the inherent sparseness of Ribo-seq data, P-site positions are merged into bins using one of three methods: 1) an adaptive bin-width method that varies by gene, based on the Freedman-Diaconis rule 2) a fixed bin width method (as small as a single codon) that the user may specify, or 3) binning  by exon, using boundaries specified in the GTF file. 

D) **Differential pattern analysis**: Pattern analysis of genes is performed on binned data for a given gene, comparing bin by bin across conditions to identify regions with statistically significant differences. The results of this testing are output as <em>p</em>-values and <em>q</em>-values for each gene. Additionally, a supplementary statistic, the <em>T</em>-value, is also produced to identify genes with a larger changes across conditions among significant genes, and is calculated based on a singular value decomposition procedure. <em>T</em>-value is intended to account for both the magnitude and number of differential bins, thus providing a way to prioritize subsets of significant genes for further investigation.

Optionally: **Visualization of Ribo-seq footprints**: RiboDiPA also provides functionality for the visualization of mapped P-site frequency data for a given gene, as well as the binned data directly used for testing, with significantly different bins marked.

## Vignette of the RiboDiPA pipeline

The following vignette is intended to provide a walkthrough for running the RiboDiPA R package, pointing out both the main workflow and optional functionality for users. It presumes that you have successfully installed RiboDiPA package from Bioconductor.

The data provided for the purposes of the vignette are adapted from Kasari et al, and represent a high-quality dataset collected in yeast. These data compare wild type cells to cells depleted for New1, which was shown by the authors to be a regulator involved in translation termination. As is often the case for data included in vignettes, the provided files are subsets of the full data set, and are intended to illustrate the functionality of RiboDiPA. We note that a typical full-scale analysis of real data for most users will be computing intensive. The computing time depends upon the number of samples, the sequencing depth of the samples, and the complexity of the organism, in terms of number of genes and exons. For example, the total computing time of the wild type versus New1 comparison (4 samples) on a 20-core node is about 10 minutes. RiboDiPA utilizes the parallel computing functionality of R and automatically detects the number of cores available to run jobs in parallel and improve performance. While a personal computer is more than sufficient for the illustration purposes of this vignette, for optimal performance with real data, we recommend that users run RiboDiPA on a server or computing cluster.

### 0. Ribo-DiPA Wrapper Function

For users' convenience, we have provided a wrapper function to permit execution of the Ribo-DiPA pipeline, which minimally requires a GTF file and BAM files, separated by experiment and replicate.

```r
## Download sample files from GitHub
sample_names <- c("WT1.bam","WT2.bam","MUT1.bam","MUT2.bam","eg.gtf")
url <- "https://github.com/mhope321/RiboDipa/raw/master/"
dest <- paste0(getwd(),"/")
for (sample_file in sample_names){
    download.file(paste0(url,sample_file),paste0(dest,sample_file))
}
```


```r
## Get BAM files from the RiboDiPA package
bam_file_list <- list.files(path=dest,
    pattern=".bam$", full.names=TRUE)
names.sample <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", bam_file_list)
```
This will produce a list of four BAM files: WT1.bam, WT2.bam, MUT1.bam, and MUT2.bam, which represent two biological replicates each of wild type cells and New1 mutant cells, respectively. These BAM files were subset on the interval chrIV:553,166-581,762 using samtools, which is a roughly 30kb region that contains 16 genes. Alternatively, users can declare the names of their BAM files directly in a vector.

```r
## Get the GTF file used for alignment
gtf_file <- list.files(path=dest,pattern=".gtf$", full.names=TRUE)
```
We recommend that users utilize the identical GTF file used to produce the experimental alignments. For example, a GTF file sourced from Ensembl will not work with BAM files aligned with a GTF file sourced from UCSC. The GTF file, "eg.gtf", provided in the package is adapted from Ensembl, Saccharomyces cerevisiae release R64-1-1, and only contains features on chromosome IV. Users may also declare their GTF file directly.

```r
## Make the class label for the experiment
classlabel <- data.frame(condition = c("mutant","mutant","wildtype","wildtype"),
    comparison=c(2,2,1,1))
rownames(classlabel) <- names.sample
```
The class label determines the comparison performed by RiboDiPA, and minimally requires a column named `comparison` which labels the reference condition "1" and the treatment condition "2", with the option of including conditions that should not be compared labeled with "0". In this case "wildtype" represents the reference condition, and "mutant" represents the treatment.

```r
## Run the RiboDiPA pipeline with default parameters
result.pip <- RiboDiPA(bam_file_list,gtf_file,classlabel)
#> Computing P-site offsets
```
The `RiboDiPA()` function is a wrapper function that calls all other necessary functions in the package. The default approaches for this wrapper are to do automatic generation of P-site offsets and adaptive binning of the mapped P-sites, although all parameters available in the other functions of the package are available to be modified in the wrapper.


```r
## View the table of output from RiboDiPA
head(result.pip$gene[order(result.pip$gene$qvalue),])
#>               tvalue       pvalue       qvalue
#> YDR050C   0.07135543 1.039896e-18 8.215177e-17
#> YDR064W   0.06267031 5.501646e-07 2.173150e-05
#> YDR062W   0.04701957 3.319300e-02 8.740823e-01
#> YDL007C-A 0.00000000 9.999990e-01 1.000000e+00
#> YDL013W   0.01329390 8.645236e-01 1.000000e+00
#> YDL017W   0.00000000 9.999989e-01 1.000000e+00
```
The `RiboDiPA()` function outputs a list of genes with their <em>T</em>-value, <em>p</em>-value, and adjusted <em>p</em>-values indicated, stored in the value `gene`, along with other intermediate data objects used in the calculation. In most cases, we expect that users will sort by the adjusted <em>p</em>-value in order to see the most significant genes genome-wide, which we show above. Genes YDR049-YDR065 are located within the interval selected for this vignette, and we can clearly see highly significant gene hits with TPI1 and RPS13 (YDR050C and YDR064W, respectively), with <em>q</em>-values of 8.22e-17 and 2.17e-05. 

This completes the bare minimum requirements for finding genes with significant differential patterns, and subsequent sections will walk through the corresponding functions in more detail.

### 1. P-site mapping
A P-site is the exact position on mRNA that has been translated by the ribosome, where the growing nascent chain of the polypeptide (covalently attached to tRNA) is located. In practice, RPFs that have been aligned to the genome have different lengths, therefore a procedure is required to map a given RPF to a P-site position to get a clear picture of ribosome translational activity. The `PsiteMapping()` function will take the input data, and use user-specified rules to map RPFs to P-sites, or generate those rules automatically using the procedure described in Lauria et al (2018). Additionally, if there are multiple transcript isoforms in a sample that utilize the same exon in the genome, in can be difficult (or impossible) to assign a given RPF to a particular transcript in a Ribo-seq experiment. RiboDiPA circumvents this problem by combining all exons into a "total transcript" and performs testing at the gene level. Therefore, in addition to the P-site offset generation and mapping, the `PsiteMapping()` function also generates total transcript coordinates for each exon.

```r
## Perform individual P-site mapping procedure
data.psite <- PsiteMapping(bam_file_list=bam_file_list, gtf_file=gtf_file, 
    psite.mapping="auto",cores=NULL)
#> Computing P-site offsets
```
P-site mapping outputs a list of values: `coverage`, the coverage across each gene; `counts`, the number of P-sites counts per gene; `exons`, the total transcript coordinates for each exon in the genome; and `psite.mapping`, the P-site mapping offsets. For the `coverage` object, rows correspond to replicates and columns correspond to nucleotide positions with respect to the total transcript. Similarly, for the `counts` object, rows represent genes and columns represent samples. Now, let us examine the offsets generated automatically by the function:

```r
## P-site mapping offset rule generated
data.psite$psite.mapping
```
The read length of the RPF is listed (`qwidth`), followed by the nucleotide offset from the 5' end (`psite`). For instance, reads of length 28 have an offset of 12, meaning that the P-site will be mapped 12 nucleotides from the 5' end of the read. 

As mentioned above, the optimal P-site offsets from the 5' end of reads is calibrated using a two-step algorithm on start codons genome-wide, closely following the procedure of Lauria et al (2018). First, for a given read length, the offset is calculated by taking the distance between the first nucleotide of the start codon and the 5' most nucleotide of the read, and then defining the offset as the 5' position with the most reads mapped to it. This process is repeated for all read lengths and then the temporary global offset is defined to be the offset of the read length with the maximum count. Lastly, for each read length, the adjusted offset is defined to be the one corresponding to the local maximum found in the profiles of the start codons closest to the temporary global offset. 

In case of insufficient reads mapped to the start codons, we recommend users to use the `center` option to take the center of the read as the P-site, or to provide their own offset rules by simply using a matrix with two columns, labeled `qwidth` and `psite`, passed into the `psite.mapping` parameter of the `PsiteMapping()` function. We note that specifying fixed rules for the P-site offsets might be especially useful when comparing across different experiments collected in the same organism, to insure consistency in the comparison.

```r
## Use user specified psite mapping offset rule
offsets <- cbind(qwidth=c(28,29,30,31,32),psite=c(18,18,18,19,19))
data.psite2 <- PsiteMapping(bam_file_list, gtf_file, psite.mapping=offsets,cores=NULL)
```

Lastly, the `PsiteMapping()` function uses the parallel computing package doParallel to speed up the process of mapping P-sites. To utilize this feature, specify the number of cores available for the job using the `cores` parameter. If `cores` is not specified, this function will automatically detect the number of cores on your computer to run jobs in parallel.

### 2. Data binning
Once reads have been mapped to P-sites in the various experiments, the next step is to bin mapped P-sites together to permit statistical testing. The smallest bin one could imagine is a single-codon (three nucleotides) which would provide the highest resolution of binning, but entails some practical problems. For instance, very long genes will have more codons, therefore after correction for multiple hypothesis testing, only the most pronounced perturbations would show statistical significance at large genes. Alternatively, the largest bin imaginable is to use an entire gene as one bin, although positional information across the gene will be lost. Therefore, a robust method to choose the right bin size per gene to permit discovery is needed, which is the essence of the RiboDiPA adaptive binning method. 

The adaptive method uses a procedure based on the Freedman-Diaconisis rule to pick an optimal number of bins of equal width for each gene, where different genes will have different bin widths, but the positions and number of bins for a gene will be the same across replicates and conditions to permit testing.    

```r
## Merge the P-site data into bins with a fixed or an adaptive width
data.binned <- DataBinning(data=data.psite$coverage, bin.width=0, zero.omit=FALSE, 
    bin.from.5UTR=TRUE, cores=NULL)
```
The function `DataBinning()` returns a list of binned P-site footprint matrices. In each matrix, rows correspond to replicates, and columns correspond to bins. If the parameter `bin.width` is not specified or set to zero, this indicates that the function will run in the adaptive binning mode (as opposed to fixed-width mode, see below). In general, we recommend to use adaptive binning, due to the fact that most Ribo-seq experiments are sparse and have few numbers of reads, on a per codon basis. 

If the `zero.omit` argument is set to `TRUE`, bins with all zeros across all replicates are removed from the differential pattern analysis. When the length of total transcript is not an integer multiple of the bin width, binning will start from the 5' end if `bin.from.5UTR` argument is `TRUE`, or from the 3' end otherwise. In general, bin width is equal for every bin across the total transcript, except for the last two bins, which are adjusted to prevent the last bin from being very small in the case where the bin width does not divide the total transcript length evenly.

```r
## Merge the P-site data on each codon
data.codon <- DataBinning(data=data.psite$coverage, bin.width=1, zero.omit=FALSE, 
    bin.from.5UTR=TRUE, cores=NULL)
```
In cases where coverage permits, users can also specify a fixed width of bin, all the way down to 1, which represents single-codon resolution. This can be useful for examining details at regions that are very likely to have changes in translational regulation, namely near the start and stop codons. For instance, we examined 50 codons upstream and downstream of the stop and start codons respectively to identify a patterns of stacked ribosomes near the stop codon in the case of New1 deletion (see Li et al, 2020).

```r
## Merge the P-site data on each exon and perform differential pattern analysis
result.exon <- DPtest.exon(psitemap=data.psite, classlabel=classlabel, 
    method=c('gtxr','qvalue'))
```
In cases where users would prefer to use exons as the bins for statistical testing, we have provided a function called `DPtest.exon()`. This function rolls data binning and differential pattern testing into one function and outputs the same structure qw `DPtest()` function. For organisms like yeast where alternative splicing is minimal, this option may not be very useful, but for higher organisms with many exons and much more alternative splicing, it may provide useful insight.

### 3. Differential pattern analysis
Once appropriate bins have been generated, the RiboDiPA package performs differential pattern testing on P-site counts bin by bin for each gene. Briefly, counts are modeled by a negative binomial distribution to call bins with statistically significant differences across conditions, and then the smallest p-value for a given gene is adjusted to control for multiple hypothesis testing across all genes. 

```r
## Perform differential pattern analysis
result.pst <- DPtest(data=data.binned, classlabel=classlabel, 
    method=c('gtxr','qvalue'))
```
The `DPtest()` function takes the output from data binning as input, and also requires a class label object, describing the comparison to be made. The class label object is minimally a data frame with two columns, `condition` and `comparison`, where `condition` labels the conditions tested, and `comparison` labels the experimental layout, where "1" indicates the control condition, "2" indicates the treatment condition, and "0" indicates replicates that should not be compared, if present.
The output of this function is a list called result, which contains a data frame object `gene` as well as other objects that contain intermediate calculations. `gene` contains <em>T</em>-value, <em>p</em>-value, and adjusted <em>p</em>-value (in this case measured by the qvalue) of all genes. The adjusted <em>p</em>-value is the main result of the testing, and therefore the main output of the package. In this case, we subset the list on genes with a <em>q</em>-value less than or equal to 0.05, which are DP positive at this threshold. The <em>p</em>-value of a gene is the minimal bin-level adjusted <em>p</em>-values for that gene. 
Additionally, the <em>T</em>-value is a supplementary statistic that quantifies the magnitude of difference between conditions, with larger numbers indicating a greater difference. The <em>T</em>-value is defined to be 1-cosine of the angle between the first right singular vectors of the footprint matrices of the two conditions under comparison. It ranges from 0-1, with larger values representing larger differences between conditions, and practically speaking, can be used to identify genes with larger magnitude of pattern difference beyond statistical significance. This might be helpful to investigators to prioritize certain genes for investigation among many that may pass the significance test for differential pattern.
Optionally, users may specify which method to use for correction for multiple hypothesis testing. The <em>q</em>-value method from `qvalue` package is the default method of multiplicity correction at gene-level, and the hybrid Hochberg-Hommel method `gtxr` from `elitism` pacakge is the default method of multiplicity correction at bin-level. Other options defined by the package `elitism` is invoked by the option to the parameter method.


### 4. Plotting
Lastly, we anticipate that users will want to plot track level data of mapped P-sites and bins for the genes that test positive for differential patterns. Therefore, we have included two plotting functions, one for mapped P-sites and one for binned data that are generated from the mapped P-sites.  The plotting is implemented with the package `ggplot2`.


```r
## Plot ribosome per nucleotide tracks of specified genes.
plot_track(data=data.psite,genes.list=c("YDR050C","YDR064W"), replicates=NULL,exons=FALSE)
```

<img src="https://github.com/jipingw/RiboDiPA-data/raw/master/fig2-1.png" title="plot of chunk fig2" alt="plot of chunk fig2" width="50%" height="50%" style="display: block; margin: auto;" /><img src="https://github.com/jipingw/RiboDiPA-data/raw/master/fig2-2.png" title="plot of chunk fig2" alt="plot of chunk fig2" width="50%" height="50%" style="display: block; margin: auto;" />
The `plot_track()` function visualizes reads mapped to P-site positions on a per gene basis. The input argument `data` is the output object of `PsiteMapping()` or `RiboDiPA()` function. The counts of RPFs mapped to P-sites is shown on the y-axis, while the total transcript in nucleotides is shown on the x-axis. Regardless of which strand the gene is located on in the genome, `plot_track()` always shows the total transcript with the 5' end on the left and the 3' end on the right. The mapped P-sites passed in the data argument, while the genes to be plotted are passed to the genes.list argument. A single gene can be plotted in this manner, as well as multiple genes with one gene per plot shown. If the exons argument is set to `TRUE`, RPFs per exon of the specified genes are also output.


```r
## Plot binned ribosome tracks of siginificant genes: YDR086C and YDR210W.
## you can specify the thrshold to redefine the significant level
plot_test(result=result.pst, genes.list=NULL, threshold=0.05) 
```
<img src="https://github.com/jipingw/RiboDiPA-data/raw/master/fig3-1.png" title="plot of chunk fig3" alt="plot of chunk fig3" width="50%" height="50%" style="display: block; margin: auto;" /><img src="https://github.com/jipingw/RiboDiPA-data/raw/master/fig3-2.png" title="plot of chunk fig3" alt="plot of chunk fig3" width="50%" height="50%" style="display: block; margin: auto;" />

The `plot_test()` function similarly visualizes RPF counts mapped to P-sites, but with data binned into the bins used for differential pattern testing. The data is passed in to the result argument, using the output object of `DPtest()` or `DPtest.exon()` or `RiboDiPA()` function. For replicates marked as "1" in `classlabel` (see `DPtest()` function), the tracks are colored blue and replicates marked as "2"  are colored red. Differential bins are colored black, with bin-level adjusted <em>p</em>-value annotated underneath the the track of the last replicate. If `genes.list` is not specified, all genes with significant differential pattern will be output. Lastly, the threshold parameter allows users to specify a threshold of signifance for choosing which genes to plots when the gene list is large. A threshold value of 0.05 will only plot genes with an adjusted <em>p</em>-value of less than or equal to 0.05.

### Conclusion
Thus concludes our vignette. For additional information, please consult the reference manual for each individual function, as well as the associated paper for this package for methodological details (Li et al, 2020)

\
\
<small>

**Other references of illustration data sets and p-site mapping**:

Kasari, V., Pochopien, A.A., Margus, T., Murina, V., Turnbull, K., Zhou, Y., Nissan, T., Graf, M., Novacek, J., Atkinson, G.C., Johansson, M.J.O,. Wilson, D.N., Hauryliuk, V. (09, 2019) "A role for the Saccharomyces cerevisiae ABCF protein New1 in translation termination/recycling." *Nucleic Acids Research*, 47(16), 8807–8820.

Lauria, F., Tebaldi, T., Bernabò, P., Groen, E.J.N., Gillingwater, T.H., Viero, G. (2018) "riboWaltz: Optimization of ribosome P-site positioning in ribosome profiling data." *PLoS Computational Biology*, 14(8):e1006169.

Freedman, D. and Diaconis, P. (1981) "On the histogram as a density estimator: L2 theory." *Zeitschrift fu ̈r Wahrscheinlichkeitstheorie und Verwandte Gebiet*e, 57(4), 453–476.

Gou, J., Tamhane, A. C., Xi, D., and Rom, D. (2014) "A class of improved hybrid Hochberg-hommel type step-up multiple test procedures." *Biometrika*, 101(4), 899–911.

Li K., Hope C.M., Wang X.A., Wang J.-P. (2020) "RiboDiPA: A novel tool for differential pattern analysis in Ribo-seq data." *Nucleic Acid Research* (in revision).



</div>
