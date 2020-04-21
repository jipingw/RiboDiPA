RiboDiPA R package 
===========================================================================================================

#### Keren Li, Matthew Hope, Xiaozhong A. Wang, Ji-Ping Wang 

#### 04/20/2020 

#### Maintainer: Ji-Ping Wang, \<[jzwang@northwestern.edu](mailto:jzwang@northwestern.edu)\>

**Reference**: Li, K., Hope, M., Wang, X., Wang, J.-P., RiboDiPA: A
novel tool for differential pattern analysis in Ribo-seq data, 2020

What is RiboDiPA?
-----------------

Ribosome profiling (also known as Ribo-seq) is an cutting-edge technique
to investigate translational activities across a wide variety of
contexts. Ribo-seq data not only provides the abundance of ribosomes
bound to transcripts for quantification of translational efficiency, but
also positional information across transcripts that could be indicative
of differences in translation dynamics between conditions. We are
particularly interested in statistical inference of translational
differences between conditions using Ribo-seq data. In this regard, most
existing tools focuse on differentiation of translational efficiency
between conditions based on the gene-wise Ribo-seq counts. The
translational efficiency for each gene is typically quantified by
Ribo-seq read count normalized by the RNA-seq gene expression. While
these approaches can provide insight into translational efficiency, the
ribosome binding pattern along the gene body, which we believe can
provide a more complete picture of transcriptional and translational
variations across conditions, has been mostly ignored.

**RiboDiPA**, short for **Ribo**some **Di**ferential **P**attern
**A**nalysis, is a bioinformatics pipeline developed for differential
pattern analysis of Ribo-seq footprint data that allows quick
identification of genes with statistically significant differences in
ribosome footprint occupancy patterns between conditions.

RiboDiPA pipline available format
---------------------------------

RiboDiPA is an R package that utilizes parallel computing functionality
with some core functions written in C++.

RiboDiPA input
--------------

1.  Ribo-seq alignment (.bam) files, one per sample.
2.  Genome Transfer File (GTF) for the reference genome.

RiboDiPA main features
----------------------

RiboDiPA R package contains three major functions.

1.  **GTF file parsing and exon merging**: A total transcript is compile
    for each gene by concatenating all exons from the same gene.

2.  **BAM file processeing and P-site mapping**: process the Ribo-seq
    alignment file (.bam) to calculate the P-site position for each RPF
    with reference to the total transcript of each gene. The Ribo-seq
    footprint, represented by the mapped P-site frequency at each
    nucleotide position along the total transcript is compiled for each
    gene of each sample.

3.  **Data binning**: merge P-site postion data into bins with a
    user-customizable fixed bin width (as small as per codon) or an
    adapte width using a built-in binning method.

4.  **Differential pattern analysis**: perform the differential pattern
    analysis on the data, output statistical significance including
    \\(p\\)-value and \\(q\\)-value together with a supplementary
    statistics \\(T\\)-value. The \\(T\\)-value can be used to identify
    genes that contain larger regions of differential codons/bins among
    the significant genes.

5.  **Visualization of Ribo-seq footprint**: RiboDiPA also provides a
    function for visualization the Ribo-seq footprint of genes.

The following diagram illustrates the flow of RiboDiPA pipeline.


<img src="http://bioinfo.stats.northwestern.edu/~jzwang/RiboDiPA/Figure2.png" width="800" height="350">

----------------------

The following vignette is intended to provide example codes for running
RiboDiPA R package. It presumes that you have successfully installed
RiboDiPA package. We illustrate below how to: 1) process bam file and
map the P-site positions with reference to total transcript for all
genes within all samples; 2) bin the mapped P-site frequency footprint;
and 3) perform the differential pattern analysis.

All three steps are computing intensive. The computing time depends upon
the number of samples, the sequencing depth and the complexity of the
organism. For example the total computing time of the WT vs. erf1d
comparison from yeast data from Wu data set (4 samples) on a 20-core
node is about 10 minutes. In the P-site mapping and data binning steps,
RiboDiPA utilizes the parallel computing functionality of R and
automatically detects the number of cores on your computer to run jobs
in parallel. Due to the large size of bam file and limited computing
power of personal computer, we recommend users to run it in servers or
computing clusters.

### 1. P-site mapping

##### Set up inputs: alignment files (.bam) and GTF files (.gtf).

``` {.sourceCode .r}
## specify bam_files from Ribo-seq, you should replace it by your own bam files
bam_file_list <- list.files(path=system.file("extdata",package="RiboDiPA"),
                        pattern=".bam$", full.names=TRUE)
names.sample <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", bam_file_list)
```

The four bam files were subsetted from a specific region of chromosome
21 from the origianl data due to package size limit. The reference for
original data files can be found in the manuscript above.

``` {.sourceCode .r}
## gtf_file you used for Ribo-seq alignment, replace it by corresponding gtf file
gtf_file <- list.files(path=system.file("extdata",package="RiboDiPA"),
                    pattern=".gtf$", full.names=TRUE)
```

##### P-site mapping

``` {.sourceCode .r}
## convert RPF reads into the P-site position on merged exons
data.psite <- PsiteMapping(bam_file_list=bam_file_list, gtf_file=gtf_file, 
                           cores=NULL)
```

The function `PsiteMapping` returns a list of two elements: `coverage`
and `counts`. `coverage`is a list of matrices, with each element
represents the Ribo-seq (P-site) footprint of a gene. Rows corrspond to
replicates and columns corrspond to nucleotide positions with reference
to the total transcript. `counts` is a mtrix of read counts where each
row stand for a gene and each column for a sample.

If `cores` is not specified, this function will automatically detect the
number of cores on your computer to run jobs in parallel.

### 2. Data binning

``` {.sourceCode .r}
## merge the P-site data into bins with a fixed or an adaptive width
data.binned <- DataBinning(data=data.psite$coverage, bin.width=1, zero.omit=F, 
                           bin.from.5UTR=T, cores=NULL)
```

This function `DataBinning` returns a list of binned P-site footprint
matrices. In each matrix, rows corrspond to replicates, columns
corrspond to codons/bins. If value of ‘bin.width’ is specified, the
footprint data is binned with ‘bin.width’ codons per bin. If ‘bin.width’
is not specified or equal `0`, an adaptive bin width is calcualted using
the Freedman-Diaconisis rule. In general we recommend to use adpative
binning due to the low read count per codon typically observed for
Ribo-seq data. Alterntively uses can specify a smaller value for
‘bin.width’ (a minimum of 1 for codon level analysis) for fine-scale
pattern analysis

When the length of total transcript is not an integer multiple of the
binning width, binning will start from the 5’ end if `bin.from.5UTR`
argument is `TRUE`, or from the 3’ end otherwise. If the `zero.omit`
argument is `TRUE`, bins with all zeros across replicates are removed
from the differential pattern analysis.

### 3. Differential pattern analysis algorithm

``` {.sourceCode .r}
## perform differential pattern analysis
coldata <- data.frame("condition"=c(rep("A",2),rep("B",2)),
                         "type"=names.sample,"comparison"=c(1,1,2,2))
rownames(coldata) <- names.sample
result.pst <- RiboDiPA(data=data.binned, genes.list=NULL, coldata=coldata, 
                   method=c('gtxr','qvalue'))
result.pst$gene
```

    #>             tvalue       pvalue         padj
    #> YCR076C 0.08024215 9.993758e-01 9.998607e-01
    #> YDL160C 0.06009355 9.242224e-01 9.998607e-01
    #> YDR086C 0.34976814 4.102696e-13 4.102696e-12
    #> YDR148C 0.08099411 3.333678e-01 9.998607e-01
    #> YDR210W 0.31982445 7.658630e-10 3.829315e-09
    #> YDR400W 0.07960790 9.991307e-01 9.998607e-01
    #> YFR047C 0.08576016 9.943428e-01 9.998607e-01
    #> YGL097W 0.08407233 9.987185e-01 9.998607e-01
    #> YGR193C 0.08001507 9.998607e-01 9.998607e-01
    #> YOL146W 0.06674621 9.976214e-01 9.998607e-01

The function `RiboDiPA` performs the differential pattern analysis. It
first normalizes the Ribo-seq footprint data withint each gene, then
pools the normalized data from all genesfor parameter estimations and
differential abundance test. For each gene, RiboDiPA outputs a
gene-level \\(p\\)-value, and further, an adjusted \\(p\\)-value with
multiple testing correction (method can be speficied) and a
\\(q\\)-value (from `qvalue` package) for false discovery rate control.

RiboDiPA also ouputs a supplementary measure called \\(T\\)-value, which
is defined to be 1-cosine of the angle between the first right singular
vectors of the footprint matrices of the two conditions under
comparison. It can be used to identify genes with larger magnitude of
pattern difference beyond statistical significance.

### 4. Plotting

``` {.sourceCode .r}
## plot binned ribosome tracks of siginificant genes: YDR086C and YDR210W.
## you can specify the thrshold to redefine the significant level
plot_track(result=result.pst, genes.list=c("YDR086C","YDR210W"), threshold=0.05) 
```

This function `plot_track` visualizes the Ribo-seq footprint of the
genes specified in the `genes.list`. For replicates marked as `1` in
`coldata` (see `RiboDiPA` function), the tracks are colored blue and
replicates marked as `2` are colored red. Differential bins are colored
black, with bin-level adjusted \\(p\\)-value annotated underneath the
the track of the last replicate. If `genes.list` is not specified, all
genes with significant differential pattern will be output.
