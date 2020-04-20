
#### p-site mapping on the merged exons ####

PsiteMapping <- function(bam_file_list, gtf_file, cores=NULL){
  options(warn=-1)

  # parse gtf file to create GRangesList object
  txdb <- makeTxDbFromGFF(gtf_file)
  all_genes <- exonsBy(txdb, by="gene")

  # exclude genes that are on two strands or genes exist on multiple chromosomes
  all_genes <- all_genes[sapply(unique(runValue(strand(all_genes))),length)==1]
  all_genes <- all_genes[sapply(unique(runValue(seqnames(all_genes))),length)==1]

  # merge overlapping exons
  all_genes <- IRanges::reduce(all_genes)

  for (k in 1:length(bam_file_list)){
    #if bam index not exist, generate index
    if(!file.exists(paste(bam_file_list[k],".bai",sep=""))){
      cat("Index", bam_file_list[k],"\n")
      indexBam(bam_file_list[k])
    }
    cat ("Computing coverage for", bam_file_list[k],"\n")
    chrom <- seqnames(seqinfo(all_genes))

    if(is.null(cores)){
      cores <- detectCores(logical = FALSE)
    }
    cl <- parallel::makeCluster(cores-1)
    registerDoParallel(cl)
    # registerDoSEQ()
    # read coverage
    i <- NULL
    results <- foreach(i=1:length(chrom),
                       .export=c(".mapping_by_ch",".IntersectionStrict2","psitecal",".psite"),
                       .packages = c("Rsamtools", "GenomicAlignments","GenomicFeatures","data.table","Rcpp"),
                       .multicombine=TRUE) %dopar%
      {
        res <- .mapping_by_ch(bam_file_list[k],all_genes,chrom[i])
        return(res)
      }

    stopImplicitCluster()
    stopCluster(cl)

    ##separate coverage and read_counts
    coverage_RLE <- RleList()    ## coverage in Rle format
    read_counts <- data.frame()  ## data frame for read-counts

    read_counts <- do.call(rbind,lapply(results,`[[`,2))
    coverage_RLE <- lapply(results,`[[`,1)

    ## convert Rle object into vector coverage
    coverage_RLE <- coverage_RLE[lengths(coverage_RLE)>0]
    coverage_RLE <- do.call('c',coverage_RLE)
    coverage_vector <- lapply(coverage_RLE,as.vector)

    cat (bam_file_list[k], "is done!", "\n")

    if (k==1){
      coverage_matrix <- coverage_vector
      read_counts_matrix <- read_counts
    }else {
      coverage_matrix <- mapply(rbind,coverage_matrix,coverage_vector)
      read_counts_matrix <- cbind(read_counts_matrix,read_counts)
    }
  }
  colnames(read_counts_matrix) <- bam_file_list
  return(list(coverage=coverage_matrix,counts=read_counts_matrix))
}


.mapping_by_ch <- function(bam,all_genes,ch){
  ##function to calculate the coverage score and return in Rle objects

  bf = BamFile(bam)
  tx_compat_cvg=RleList()
  read_counts=data.frame()

  ## slow
  genes=all_genes[sapply(runValue(seqnames(all_genes)),unique)==ch]

  if(length(genes)){
    gr = as(seqinfo(bf), "GRanges")
    if(ch %in% runValue(seqnames(gr))){
      param = ScanBamParam(which = gr[ch],flag=scanBamFlag(isSecondaryAlignment=FALSE))
      gal2=readGAlignments(bf, param = param)
      gal2=gal2[strand(gal2)!="*"]
      gal2=gal2[is.na(seqnames(gal2))==F]

      if(length(gal2)){
        gal3 <- data.frame(gal2)
        gal3.psite <- .psite(gal3$qwidth)
        gal3.psite[gal3$strand=='-'] <- gal3$qwidth[gal3$strand=='-']-gal3.psite[gal3$strand=='-']-1

        center <- psitecal(gal3$cigar,gal3$start,gal3.psite)
        gal3$qwidth <- 1
        gal3$cigar <- "1M"
        gal3$start <- center
        gal3$end <- gal3$start
        gal3$width <- 1
        gal3$njunc <- 0

        gal3 <- makeGRangesFromDataFrame(gal3)
        gal3 <- as(gal3,"GAlignments")
      }else{
        gal3 <- gal2
      }
      #only keep genes on the current chromosome
      results <- .IntersectionStrict2(genes,gal3)
      gal3 <- gal3[queryHits(results)]
      results <- .IntersectionStrict2(genes,gal3)

      tx2reads <- setNames(as(t(results), "List"), names(genes))
      compat_reads_by_tx <- extractList(gal3, tx2reads)
      read_counts=data.frame(lengths(tx2reads))
      tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                             genes,
                                             ignore.strand=FALSE)
    }else{
      tx2reads1 <- IntegerList(vector("list", length(genes)))
      names(tx2reads1) <- names(genes)
      compat_reads_by_tx <- extractList(gr, tx2reads1)
      read_counts=data.frame(lengths(tx2reads1))
      tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                             genes,
                                             ignore.strand=FALSE)
    }
  }
  return(list(tx_compat_cvg,read_counts))
}


.IntersectionStrict2 <-  function (features, reads, ignore.strand = FALSE, inter.feature = TRUE)
{
  ##function adapted to compile compatible reads (single or paired-end) as HTseq
  ov <- findOverlaps(reads, features, type = "within", ignore.strand = ignore.strand)
  if (inter.feature) {
    reads_to_keep <- which(countQueryHits(ov) == 1L)
    ov <- ov[queryHits(ov) %in% reads_to_keep]
  }
  return(ov)
}

.psite <- function(x){
  (x-1)*(x<=13)+12*(x==20|x==27|x==28)+13*(x>13&x!=20&x!=27&x!=28)
}

#### data binning ###

.app2exact.self <- function(n, p) {
  # p=p[1:m] is an allocation, p_i >=0, sum_i p_i = 1
  ans=floor(n*p);
  k=n-sum(ans);
  if(k>0) {
    stemp=n*p-ans;
    otemp=order(-stemp)[1:k];
    ans[otemp]=ans[otemp]+1;
  };
  return(ans)
}

DataBinning <- function(data, bin.width=0, zero.omit=F, bin.from.5UTR=T, cores=NULL){
  if(is.null(cores)){
    cores <- detectCores(logical = FALSE)
  }
  cl <- makeCluster(cores-1)
  registerDoParallel(cl)

  genes.list <- names(data)
  gene <- NULL
  DATA.BIN <- foreach(gene = genes.list,.final = function(x) setNames(x, genes.list),
                      .packages ="reldist", .export=".app2exact.self")%dopar%
    {
      data1 <- data[[gene]]
      n <- nrow(data1)
      p <- ncol(data1)
      # merge every 3 nt into a codon
      p.codon <- ceiling(p/3)
      codon.size <- .app2exact.self(p,rep(1/p.codon,p.codon))
      if(!bin.from.5UTR){
        codon.size <- rev(codon.size)
      }
      codon.cum <- c(0,cumsum(codon.size))
      data.codon <- matrix(0,nrow=n,ncol=p.codon)
      for(j in 1:p.codon){
        data.codon[,j] <- rowSums(matrix(data1[,(codon.cum[j]+1):(codon.cum[j+1])],nrow=n))
      }

      # adpative or fixed bin width
      if(bin.width==1){
        data.bin <- data.codon
      }else if(sum(data1)){
        if(bin.width){
          # fixed
          p.bin <- ceiling(p.codon/bin.width)
        }else{
          # adaptive
          n.data <- mean(rowSums(data.codon))
          iqr.data <- abs(wtd.iqr(x=1:p.codon,weight=colMeans(data.codon)))
          n.bin <- ceiling(p.codon/(2*iqr.data/(n.data^(1/3))))
          p.bin <- min(p.codon, n.bin)
        }

        bin.size <- .app2exact.self(ncol(data.codon),rep(1/p.bin,p.bin))
        if(!bin.from.5UTR){
          bin.size <- rev(bin.size)
        }
        bin.cum <- c(0,cumsum(bin.size))
        data.bin <- matrix(0,nrow=n,ncol=p.bin)
        for(j in 1:p.bin){
          data.bin[,j] <- rowSums(matrix(data.codon[,(bin.cum[j]+1):(bin.cum[j+1])],nrow=n))
        }
      }else{
        data.bin <- matrix(0,nrow=n,ncol=1)
      }

      # remove common zero bins
      if(zero.omit){
        data.bin <- data.bin[,colSums(data.bin)>0]
      }
      if(is.vector(data.bin)){
        data.bin <- matrix(data.bin,ncol=1)
      }
      data.bin;
    }
  stopCluster(cl)
  return(DATA.BIN)
}



#### svd function ####
.svdv1 <- function(x){
  ifelse(is.vector(x),
         result <- x/sqrt(sum(x^2)),
         result <- svd(x)$v[,1]
  )
  return(result)
}

#### abundance estimation function ####
SizeFactors <- function(x, condition){
  s.tmp <- rowSums(x)/median(rowSums(x))
  if(all(s.tmp>0)){
    y1 <- colSums2(x,rows=(condition==1))
    y2 <- colSums2(x,rows=(condition==2))
    ind1 <- (y1>0)&(y2>0)
    if(sum(ind1)>1){
      x <- x[,ind1]
      y1 <- colSums2(x,rows=(condition==1))
      y2 <- colSums2(x,rows=(condition==2))
      x.total <- sum(x)
      tmp <- log(y1)-log(y2)
      tmp2 <- as.numeric(quantile(tmp,probs=c(.25,.75), na.rm = T)+c(-1.5,1.5)*IQR(tmp, na.rm = T))
      ind2 <- tmp>=tmp2[1] &tmp<=tmp2[2]
      if(sum(ind2)>1){
        x <- x[,ind2]
        if(sum(x)/x.total>0.1&all(rowSums(x)>0)){
          s.tmp <- rowSums(x)/median(rowSums(x))
        }
      }
    }
  }
  s.tmp[s.tmp==0] <- 1e-8
  s.tmp;
}

#### main function ####
RiboDiPA <- function(data, coldata, genes.list=NULL,  method=c('gtxr','qvalue')){
  options(warn=-1)

  ## prepare data
  if(is.null(genes.list)){
    genes.list <- names(data)
  }else{
    data <- data[genes.list]
  }

  data <- lapply(data,'[', which(coldata$comparison!=0),)
  noread <- which(do.call('c',lapply(data,is.vector)))
  genes.list <- genes.list[which(do.call('c',lapply(data,is.matrix)))]
  data <- data[genes.list]

  condition <- coldata$comparison[coldata$comparison!=0]

  ## t-value
  tvalue <- 1-do.call('c',lapply(data, function(x) abs(sum(.svdv1(x[which(condition==1),])*.svdv1(x[which(condition==2),])))))

  ## size factor
  S.hat <- lapply(data,SizeFactors, condition=condition)
  normFactors <- t(do.call('cbind',mapply(function(x,y) replicate(ncol(y),x),S.hat,data)))

  # link all codons
  DATA.com <- t(do.call('cbind',data))

  # codon level test
  coldata <- coldata[coldata$comparison!=0,]
  coldata$comparison <- as.factor(coldata$comparison)
  colnames(DATA.com) <- coldata$type
  dds <- DESeqDataSetFromMatrix(countData = DATA.com,
                                colData = coldata,
                                design = ~ comparison)
  # assign size factors
  normalizationFactors(dds) <- normFactors
  dds <- DESeq(dds)
  res <- as.matrix(results(dds))[,c("pvalue","log2FoldChange")]

  # split test results by genes
  LL <- unlist(lapply(data,ncol))
  cLL <- c(0,cumsum(LL))
  res.codon <- NULL
  pvalue.gene <- NULL
  for(l in 1:length(LL)){
    tmp <- res[(cLL[l]+1):(cLL[l]+LL[l]),]
    tmp <- cbind(tmp,padj =NA)
    tmp[which(!is.na(tmp[,"pvalue"])),"padj"] <- elitism::p.adjust(na.omit(tmp[,"pvalue"]), method=method[1])
    res.codon[[names(LL)[l]]] <- tmp
    pvalue.gene <- c(pvalue.gene,min(tmp[,"padj"],na.rm = T))
  }
  names(pvalue.gene) <- names(LL)

  # genome wide family control
  if(method[2]=="qvalue"){
    padj.gene <- qvalue(pvalue.gene)$qvalues
  }else{
    padj.gene <- elitism::p.adjust(pvalue.gene, method=method[2])
  }

  return(list(codon=res.codon[genes.list],
       gene=data.frame(cbind(tvalue=tvalue,pvalue=pvalue.gene,padj=padj.gene))[genes.list,],
       small=rownames(noread),coldata=coldata,data=data))
}

#### plotting ####
.reduceself <- function(ints){
  for(x in 1:nrow(ints)){
    for(y in 1:nrow(ints)){
      if(data.table::between(ints[x,1], ints[y,1], ints[y,2])){
        ints[x,1] <- ints[y,1]
        if(ints[y,2] > ints[x,2]){
          ints[x,2] <- ints[y,2]
        } else {
          ints[y,2] <- ints[x,2]
        }
      }
    }
  }
  tmp <- unique(ints, margin = 1)
  res <- tmp[order(tmp[,1]),]
  if(is.vector(res)){
    res <- matrix(res,nrow=1)
  }
  return(res)
}

.segcal <- function(x,ints){
  l <- sum(x>ints[,2])
  if(!l){
    result <- x-ints[1,1]
  }else{
    result <- diff(colSums(ints[1:l,]))+(x-ints[l+1,1])
  }
  if(!result) result <- 1
  result;
}

plot_track <- function(data, genes.list,replicates=NULL,gtf_file){
  options(warn=-1)
  txdb <- makeTxDbFromGFF(gtf_file)
  all_genes <- exonsBy(txdb, by="gene")
  if(is.null(replicates)){
    replicates <- 1:ncol(data$counts)
  }
  par(oma = c(1,4,3,0) + 0.5,
      mar = c(4,1,1,1) + 0.5)
  layout(matrix(1:length(replicates), ncol=1))
  for(gene in genes.list){
    ranges(all_genes[[gene]])

    data1 <- data$coverage[[gene]][replicates,]
    n <- nrow(data1)
    for(i in replicates){
        barplot(data1[i,], main=colnames(data$counts)[i],col = "black",cex.axis=2,
                border=NA,ylim=c(-max(data1[i,])/20,max(data1[i,])))
    }

    u <- par("usr")
    code.tmp <- cbind(start(all_genes[[gene]]),end(all_genes[[gene]]))
    code.reduced2 <- .reduceself(code.tmp)
    code.reduced2 <- .reduceself(code.reduced2)

    tmtmp <- NULL
    for(i in 1:nrow(code.tmp)){
      tmtmp <- rbind(tmtmp,c(start=.segcal(code.tmp[i,1],code.reduced2),
        end=.segcal(code.tmp[i,2],code.reduced2)))
    }
    tmtmp2 <- cbind(unlist(tmtmp[,1])*u[2]*0.96/max(as.numeric(tmtmp)),
                    unlist(tmtmp[,2])*u[2]*0.96/max(as.numeric(tmtmp)))
    for(i in 1:nrow(tmtmp2)){
      rect(tmtmp2[i,1],-max(data1[n,])/40, tmtmp2[i,2],0,col='green')
    }

    title(gene,cex.main=2, outer=TRUE)
  }
}

plot_test <- function(result, genes.list=NULL,threshold=0.05){
  condition <- result$coldata$comparison
  if(is.null(genes.list)){
    genes.list <- rownames(result$gene[which(result$gene$padj<=threshold),])
  }
  par(oma = c(1,4,3,0) + 0.5,
      mar = c(4,1,1,1) + 0.5)
  layout(matrix(1:sum(condition!=0), ncol=1))
  for(gene in genes.list){

    data1 <- result$data[[gene]]
    p <- ncol(data1)
    tmp <- result$codon[[gene]][,'padj']
    col.red <- rep("red",p)
    col.red[which(tmp<=threshold)] <- "black"
    col.blue <- rep("blue",p)
    col.blue[which(tmp<=threshold)] <- "black"
    tmp[which(tmp>threshold)] <- NA
    tmp[which(tmp<=threshold)] <- formatC(tmp[which(tmp<=threshold)], format = "e", digits = 2)
    for(i in 1:length(condition)){
      if(condition[i]==1){
        mids <- barplot(data1[i,], main=rownames(result$coldata)[i],col = col.blue,cex.axis=2,border=NA)
      }
    }
    for(i in 1:length(condition)){
      if(condition[i]==2){
        mids <- barplot(data1[i,], main=rownames(result$coldata)[i],col = col.red,cex.axis=2,border=NA)
      }
    }
    axis(1, at=mids,labels=tmp,las=3)
    title(paste(gene,", p-value =",formatC(result$gene[gene,"padj"], format = "e", digits = 2), " T-value =",
                formatC(result$gene[gene,"tvalue"], format = "e", digits = 2), sep=""),cex.main=2, outer=TRUE)
  }
}
