
#### p-site mapping on the merged exons ####

PsiteMapping <- function(bam_file_list, gtf_file, psite.mapping="auto",cores=NULL){
  options(warn=-1)
  names.sample <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", bam_file_list)

  # parse gtf file to create GRangesList object
  txdb <- makeTxDbFromGFF(gtf_file)
  all_genes <- exonsBy(txdb, by="gene")

  # exclude genes that are on two strands or genes exist on multiple chromosomes
  all_genes <- all_genes[sapply(unique(runValue(strand(all_genes))),length)==1]
  all_genes <- all_genes[sapply(unique(runValue(seqnames(all_genes))),length)==1]

  exons.merge <- NULL
  for(gene in names(all_genes)){
    single.gene <- all_genes[[gene]]
    range.gene <- cbind(start(single.gene),end(single.gene))
    range.gene2 <- .reduceself(range.gene)
    range.gene2 <- .reduceself(range.gene2)

    range.relative <- NULL
    for(i in 1:nrow(range.gene)){
      range.relative <- rbind(range.relative,c(start=.segcal(range.gene[i,1],range.gene2),
                             end=.segcal(range.gene[i,2],range.gene2)))
    }
    rownames(range.relative) <- single.gene$exon_name
    exons.merge[[gene]] <- range.relative
  }

  # merge overlapping exons
  all_genes <- IRanges::reduce(all_genes)

  for (k in 1:length(bam_file_list)){
    #if bam index not exist, generate index
    if(!file.exists(paste(bam_file_list[k],".bai",sep=""))){
      cat("Index", bam_file_list[k],"\n")
      indexBam(bam_file_list[k])
    }
    chrom <- seqnames(seqinfo(all_genes))
  }

  if(psite.mapping=="auto"){
    cat ("Computing P-site offsets \n")
    psite.mapping <- .psite_ribowaltz(bam_file_list, txdb)
  }

  for (k in 1:length(bam_file_list)){
    cat ("Computing coverage for", bam_file_list[k],"\n")
    if(is.null(cores)){
      cores <- detectCores(logical = FALSE)
    }
    cl <- parallel::makeCluster(cores-1)
    registerDoParallel(cl)
    # read coverage
    i <- NULL
    results <- foreach(i=1:length(chrom),
                       .export=c(".mapping_by_ch",".IntersectionStrict2","psitecal",".create_anno",".psite_ribowaltz"),
                       .packages = c("Rsamtools", "GenomicAlignments","GenomicFeatures","data.table","Rcpp"),
                       .multicombine=TRUE) %dopar%
      {
        res <- .mapping_by_ch(bam_file_list[k], all_genes, chrom[i], psite.mapping)
        return(res)
      }
    stopImplicitCluster()
    stopCluster(cl)

    ##separate coverage and read_counts
    coverage_RLE <- RleList()    ## coverage in Rle format
    read_counts <- data.frame()  ## data.frame for read-counts

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
  colnames(read_counts_matrix) <- names.sample
  return(list(coverage=coverage_matrix,counts=read_counts_matrix,
              exons=exons.merge,psite.mapping=psite.mapping))
}


.create_anno <-  function(txdb) {
  group_name <- NULL
  exon <- suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "tx",use.names=T))
  utr5<- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdb,use.names=T))
  cds <- suppressWarnings(GenomicFeatures::cdsBy(txdb, by = "tx", use.names=T))
  utr3<- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdb,use.names=T))

  exon <- as.data.table(exon[unique(names(exon))])
  utr5 <- as.data.table(utr5[unique(names(utr5))])
  cds <- as.data.table(cds[unique(names(cds))])
  utr3 <-as.data.table(utr3[unique(names(utr3))])

  anno_df <- exon[, list(l_tr = sum(width)), by = list(transcript = group_name)]
  l_utr5 <- utr5[, list(l_utr5 = sum(width)), by = list(transcript = group_name)]
  l_cds <- cds[, list(l_cds = sum(width)), by = list(transcript = group_name)]
  l_utr3 <- utr3[, list(l_utr3 = sum(width)), by = list(transcript = group_name)]

  merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
  anno_df  <-  Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
  anno_df[is.na(anno_df)] <- 0

  return(anno_df)
}

.psite_ribowaltz <- function(bam_file_list,txdb){
  cds_start <- cds_stop <- count <- end3 <- end5 <- i.l_cds <- .SD <- NULL
  i.l_utr5 <- offset <- onsite <- site_dist_end3 <- site_dist_end5 <- transcript <- NULL
  # prepare annotation
  annotation <- .create_anno(txdb)
  all_tx <- exonsBy(txdb, by="tx",use.names=T)
  name.tx <- names(unlist(all_tx))
  all_tmp <- as.data.table(all_tx)
  all_tmp[,start:=start-30][,end:=end+30][,width:=width+60]
  all_tmp <- GRanges(all_tmp)
  names(all_tmp) <- name.tx

  site_sub <- NULL
  for(bam in bam_file_list){
    data <- GenomicAlignments::readGAlignments(bam)
    qwidth.data <- qwidth(data)
    data <- GenomicRanges::granges(data)
    data <- as.data.table(GenomicFeatures::mapToTranscripts(data, all_tmp))
    data[,start:=start-30][,end:=end-30]
    data$width <- qwidth.data[data$xHits]

    data <- data[, c("seqnames", "start", "end", "width"), with=FALSE]
    setnames(data, c("transcript", "end5", "end3","length"))

    nreads <- nrow(data)
    data <- data[as.character(transcript) %in% as.character(annotation$transcript)]

    data[annotation, on = "transcript", `:=`(c("cds_start",
                                               "cds_stop"), list(i.l_utr5 + 1, i.l_utr5 + i.l_cds))]
    data[cds_start == 1 & cds_stop == 0, `:=`(cds_start, 0)]

    lev <- sort(unique(data$length))
    data[, site_dist_end5 := end5 - cds_start]
    data[, site_dist_end3 := end3 - cds_start]
    site_sub <- rbind(site_sub, data[site_dist_end5 <= -6 & site_dist_end3 >= 5])
  }
  site_sub <- site_sub[length==end3-end5+1]
  minlen <- min(site_sub$length)
  maxlen <- max(site_sub$length)
  t <- table(factor(site_sub$length, levels = lev))

  # offset
  offsets <- data.table(length = as.numeric(as.character(names(t))), count = as.vector(t))
  offsets[, onsite :=T][count == 0, onsite := F]
  tempoff <- function(v_dist) {
    ttable <- sort(table(v_dist), decreasing = T)
    ttable_sr <- ttable[as.character(as.numeric(names(ttable)) + 1)]
    ttable_sl <- ttable[as.character(as.numeric(names(ttable)) - 1)]
    tsel <- rowSums(cbind(ttable > ttable_sr, ttable > ttable_sl), na.rm = T)
    return(as.numeric(names(tsel[tsel == 2][1])))
  }

  offset5 <- site_sub[, list(offset = tempoff(.SD$site_dist_end5)), by = length]
  offsets.tmp  <-  merge(offsets, offset5, all.x = TRUE, by = "length")
  best.offset <- as.numeric(offsets.tmp[!is.na(offset), list(count = sum(count)), by = offset
                                        ][count == max(count)][, offset])
  # adjusted offset
  adj.off <- function(dtsite, bestoff){
    t <- table(factor(dtsite, levels = seq(min(dtsite) - 2, max(dtsite) + 1)))
    t[1:2] <- t[3] + 1
    locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == -2)]))) + 1
    adjoff <- locmax[which.min(abs(locmax - bestoff))]
    ifelse(length(adjoff) != 0, adjoff, bestoff)
  }
  adj.tab <- site_sub[, list(offset = adj.off(site_dist_end5, best.offset)), by = length]
  offsets <- merge(offsets, adj.tab, all.x = TRUE, by = "length")
  offsets[is.na(offset), offset := best.offset
          ][, offset := abs(offset)][, count := NULL][,onsite:=NULL]

  setnames(offsets,c("qwidth","psite"))
  return(offsets)
}


.mapping_by_ch <- function(bam, all_genes, ch, psite.mapping){
  ##function to calculate the coverage score and return in Rle objects
  psite <- center <- NULL
  bf <-  BamFile(bam)
  tx_compat_cvg <- RleList()
  read_counts <- data.frame()

  genes <- all_genes[sapply(runValue(seqnames(all_genes)),unique)==ch]

  if(length(genes)){
    gr  <-  as(seqinfo(bf), "GRanges")
    if(ch %in% runValue(seqnames(gr))){
      param <-  ScanBamParam(which = gr[ch],flag=scanBamFlag(isSecondaryAlignment=FALSE))
      gal2 <- readGAlignments(bf, param = param)
      gal2 <- gal2[strand(gal2)!="*"]
      gal2 <- gal2[is.na(seqnames(gal2))==F]

      if(length(gal2)){
        gal3 <- as.data.table(gal2)
        if(all(psite.mapping=="center")){
          gal3 <- gal3[,psite:=floor(qwidth/2)]
        }else{
          gal3 <- merge.data.table(gal3,psite.mapping,by="qwidth")
        }
        gal3[strand=='-',psite:=qwidth-psite-1]
        gal3[,center:=psitecal(cigar,start,psite)]
        gal3[,qwidth := 1][,width := 1][,njunc := 0][,cigar := "1M"][,start:=center
                         ][,end:=center][,psite:=NULL][,center:=NULL]
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
      read_counts <- data.frame(lengths(tx2reads))
      tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                             genes,
                                             ignore.strand=FALSE)
    }else{
      tx2reads1 <- IntegerList(vector("list", length(genes)))
      names(tx2reads1) <- names(genes)
      compat_reads_by_tx <- extractList(gr, tx2reads1)
      read_counts <- data.frame(lengths(tx2reads1))
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
    result <- x-ints[1,1]+1
  }else{
    result <- sum(ints[1:l,2])-sum(ints[1:l,1])+(x-ints[l+1,1])+1
  }
  return(result)
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
RiboDiPA <- function(data, classlabel, method=c('gtxr','qvalue')){
  options(warn=-1)

  ## prepare data
  genes.list <- names(data)
  data <- lapply(data,'[', which(classlabel$comparison!=0),)
  noread <- which(do.call('c',lapply(data,is.vector)))
  genes.list <- genes.list[which(do.call('c',lapply(data,is.matrix)))]
  data <- data[genes.list]

  condition <- classlabel$comparison[classlabel$comparison!=0]

  ## t-value
  tvalue <- 1-do.call('c',lapply(data, function(x) abs(sum(.svdv1(x[which(condition==1),])*.svdv1(x[which(condition==2),])))))

  ## size factor
  S.hat <- lapply(data,SizeFactors, condition=condition)
  normFactors <- t(do.call('cbind',mapply(function(x,y) replicate(ncol(y),x),S.hat,data)))

  # link all codons
  DATA.com <- t(do.call('cbind',data))

  # codon level test
  classlabel <- classlabel[classlabel$comparison!=0,]
  classlabel$comparison <- as.factor(classlabel$comparison)
  colnames(DATA.com) <- classlabel$type
  dds <- DESeqDataSetFromMatrix(countData = DATA.com,
                                colData = classlabel,
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
    colnames(tmp)[3] <- method[1]
    res.codon[[names(LL)[l]]] <- tmp
    pvalue.gene <- c(pvalue.gene,min(tmp[,method[1]],na.rm = T))
  }
  names(pvalue.gene) <- names(LL)

  # genome wide family control
  if(method[2]=="qvalue"){
    padj.gene <- qvalue(pvalue.gene)$qvalues
  }else{
    padj.gene <- elitism::p.adjust(pvalue.gene, method=method[2])
  }
  gene.result <- data.frame(cbind(tvalue=tvalue,pvalue=pvalue.gene,padj=padj.gene))[genes.list,]
  colnames(gene.result)[3] <- method[2]
  return(list(codon=res.codon[genes.list],
       gene=gene.result,small=rownames(noread),classlabel=classlabel,data=data,method=method))
}

#### plotting ####

plot_track <- function(data, genes.list,replicates=NULL,exons=FALSE){
  options(warn=-1)
  if(is.null(replicates)){
    replicates <- 1:ncol(data$counts)
  }
  par(oma = c(1,4,3,0) + 0.5,
      mar = c(4,1,1,1) + 0.5)
  layout(matrix(1:length(replicates), ncol=1))
  for(gene in genes.list){
    data1 <- data$coverage[[gene]]
    for(i in replicates){
        barplot(data1[i,], main=colnames(data$counts)[i],col = "black",cex.axis=2,
                border=NA,ylim=c(-max(data1[i,])/20,max(data1[i,])),space=0)
    }

    range.relative <- data$exons[[gene]]
    for(i in 1:nrow(range.relative)){
      rect(range.relative[i,1]-1,-max(data1[replicates,])/20, range.relative[i,2],0,col='green')
    }
    title(paste("gene: ",gene,sep=""),cex.main=2, outer=TRUE)

    if(exons==TRUE){
      exons.gene <- data$exons[[gene]]
      for(j in 1:nrow(exons.gene)){
        data2 <- data1[,(exons.gene[j,1]):(exons.gene[j,2])]
        for(i in replicates){
          barplot(data2[i,], main=colnames(data$counts)[i],col = "black",cex.axis=2,
                  border=NA,ylim=c(-max(data2[i,])/20,max(data2[i,])),space=0)
        }
        title(paste("gene: ",gene,', exon: ', rownames(exons.gene)[j],sep=""),cex.main=2, outer=TRUE)
      }
    }
  }
}

plot_test <- function(result, genes.list=NULL,threshold=0.05){
  method <- result$method
  condition <- result$classlabel$comparison
  if(is.null(genes.list)){
    genes.list <- rownames(result$gene[which(result$gene[,method[2]]<=threshold),])
  }
  par(oma = c(1,4,3,0) + 0.5,
      mar = c(4,1,1,1) + 0.5)
  layout(matrix(1:sum(condition!=0), ncol=1))
  for(gene in genes.list){

    data1 <- result$data[[gene]]
    p <- ncol(data1)
    tmp <- result$codon[[gene]][,method[1]]
    col.red <- rep("red",p)
    col.red[which(tmp<=threshold)] <- "black"
    col.blue <- rep("blue",p)
    col.blue[which(tmp<=threshold)] <- "black"

    for(i in 1:length(condition)){
      if(condition[i]==1){
        mids <- barplot(data1[i,], main=rownames(result$classlabel)[i],col = col.blue,cex.axis=2,border=NA,space=0)
      }
    }
    for(i in 1:length(condition)){
      if(condition[i]==2){
        mids <- barplot(data1[i,], main=rownames(result$classlabel)[i],col = col.red,cex.axis=2,border=NA,space=0)
      }
    }
    at.plot <- mids[which(tmp<=threshold)]
    label.plot <- formatC(tmp[which(tmp<=threshold)], format = "e", digits = 2)
    if(tmp[1]>threshold){
      at.plot <- c(mids[1], at.plot)
      label.plot <- c(NA, label.plot)
    }
    if(tmp[p]>threshold){
      at.plot <- c(at.plot,mids[p])
      label.plot <- c(label.plot,NA)
    }

    axis(1, at=at.plot,labels=label.plot,las=3, cex.axis=0.8)
    title(paste(gene,", ",method[2]," =",formatC(result$gene[gene,method[2]], format = "e", digits = 2), ", T-value =",
                formatC(result$gene[gene,"tvalue"], format = "e", digits = 2), sep=""),cex.main=2, outer=TRUE)
  }
}
