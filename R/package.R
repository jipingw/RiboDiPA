#### wrap function ####
RiboDiPA <- function(bam_file_list, gtf_file, classlabel, psite.mapping="auto",
    exon.binning=FALSE, bin.width=0, zero.omit=FALSE, bin.from.5UTR=TRUE,
    method=c("gtxr", "qvalue"), cores=NULL) {
    if (is.null(cores)) {
        cores <- detectCores(logical=FALSE)
    }
    data.psite <- psiteMapping(bam_file_list=bam_file_list, gtf_file=gtf_file,
        psite.mapping=psite.mapping, cores=cores)
    ifelse(exon.binning, {
        result <- diffPatternTestExon(psitemap=data.psite, 
            classlabel=classlabel, method=method)
        }, {
        data.binned <- dataBinning(data=data.psite$coverage,
            bin.width=bin.width, zero.omit=zero.omit,
            bin.from.5UTR=bin.from.5UTR, cores=cores)
        result <- diffPatternTest(data=data.binned, classlabel=classlabel,
            method=method)
    })
    return(c(result, data.psite))
}

#### p-site mapping on the merged exons ####

psiteMapping <- function(bam_file_list, gtf_file, psite.mapping="auto",
    cores=NULL) {
    options(warn=-1)
    # Get names of samples from BAM file list, removing file paths and suffixes
    names.sample <- sub("(.*\\/)([^.]+)(\\.[[:alnum:]]+$)", "\\2", 
        bam_file_list)
    names.sample <- sub(".bam", "", names.sample)
    # Parse gtf file to create GRangesList object
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)
    all_exons <- GenomicFeatures::exons(txdb, 
        columns=c("EXONNAME","TXNAME","GENEID"), use.names = TRUE)
    all_genes <- split(all_exons, unlist(values(all_exons)$GENEID))
    # exclude genes on two strands
    all_genes <- all_genes[vapply(runValue(strand(all_genes)), length,
        integer(1)) == 1]
    # exclude genes on multiple chromosomes
    all_genes <- all_genes[vapply(runValue(seqnames(all_genes)), length,
        integer(1)) == 1]
    message("Computing total transcript coordinate for exons ...")
    exons.coord <- lapply(all_genes, .totalTransciptCoordinate)
    message("done\n")
    # merge overlapping exons in total transcript
    all_genes <- IRanges::reduce(all_genes)
    # p-site mapping
    results_all <- .pMappingAll(bam_file_list=bam_file_list,
        psite.mapping=psite.mapping, txdb=txdb, all_genes=all_genes,
        cores=cores)
    colnames(results_all$counts) <- names.sample
    # assign genomic coordinates as names of sites
    tmp <- mapply(.coordGene, results_all$coverage, exons.coord)
    results_all$coverage <- tmp
    results_all$exons <- exons.coord
    return(results_all)
}

.coordGene <- function(coverage.gene,coord.gene){
    strands <- coord.gene[,"strand"][1]
    if(strands == "+"){
        coord.tmp <- IRanges(coord.gene[,"start_genome"],
            coord.gene[,"end_genome"])
        coord.tmp <- BiocGenerics::union(coord.tmp, coord.tmp)
        result <- NULL
        for(i in seq_len(length(coord.tmp))){
            result <- c(result,seq(from=start(coord.tmp)[i], 
                to=end(coord.tmp)[i]))
        }
        colnames(coverage.gene) <- result
    }else{
        coord.tmp <- IRanges(coord.gene[,"end_genome"],
            coord.gene[,"start_genome"])
        coord.tmp <- BiocGenerics::union(coord.tmp, coord.tmp)
        result <- NULL
        for(i in seq_len(length(coord.tmp))){
            result <- c(result,seq(from=start(coord.tmp)[i], 
                to=end(coord.tmp)[i]))
        }
        colnames(coverage.gene) <- rev(result)
    }
    return(coverage.gene)
}

.createAnno <- function(txdb) {
    group_name <- NULL
    exon <- suppressWarnings(GenomicFeatures::exonsBy(txdb,
        by="tx", use.names=TRUE))
    utr5 <- suppressWarnings(GenomicFeatures::fiveUTRsByTranscript(txdb,
        use.names=TRUE))
    cds <- suppressWarnings(GenomicFeatures::cdsBy(txdb, by="tx",
        use.names=TRUE))
    utr3 <- suppressWarnings(GenomicFeatures::threeUTRsByTranscript(txdb,
        use.names=TRUE))

    exon <- data.table::as.data.table(exon[unique(names(exon))])
    utr5 <- data.table::as.data.table(utr5[unique(names(utr5))])
    cds <- data.table::as.data.table(cds[unique(names(cds))])
    utr3 <- data.table::as.data.table(utr3[unique(names(utr3))])

    anno_df <- exon[, list(l_tr=sum(width)), by=list(transcript=group_name)]
    l_utr5 <- utr5[, list(l_utr5=sum(width)), by=list(transcript=group_name)]
    l_cds <- cds[, list(l_cds=sum(width)), by=list(transcript=group_name)]
    l_utr3 <- utr3[, list(l_utr3=sum(width)), by=list(transcript=group_name)]

    merge_allx <- function(x, y) merge(x, y, all.x=TRUE)
    anno_df <- Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
    # standardGeneric for 'Reduce' defined from package 'BiocGenerics'
    anno_df[is.na(anno_df)] <- 0

    return(anno_df)
}

.psiteOffset <- function(bam_file_list, txdb) {
    cds_start <- cds_stop <- end3 <- end5 <- i.l_cds <- offset <- NULL
    i.l_utr5 <-  site_dist_end3 <- site_dist_end5 <- transcript <- NULL
    # prepare annotation
    annotation <- .createAnno(txdb)
    all_tx <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)
    name.tx <- names(unlist(all_tx))
    all_tmp <- data.table::as.data.table(all_tx)
    all_tmp[, `:=`(start, start - 30)][, `:=`(end, end + 30)][, `:=`(width,
        width + 60)]
    all_tmp <- GRanges(all_tmp)
    names(all_tmp) <- name.tx

    site_sub <- NULL
    for (bam in bam_file_list) {
        data <- GenomicAlignments::readGAlignments(bam)
        qwidth.data <- qwidth(data)
        data <- GenomicRanges::granges(data)
        data <- as.data.table(GenomicFeatures::mapToTranscripts(data, all_tmp))
        data[, `:=`(start, start - 30)][, `:=`(end, end - 30)]
        data$width <- qwidth.data[data$xHits]

        data <- data[, c("seqnames", "start", "end", "width"), with=FALSE]
        setnames(data, c("transcript", "end5", "end3", "length"))

        nreads <- nrow(data)
        data <- data[as.character(transcript) %in%
            as.character(annotation$transcript)]
        data[annotation, on="transcript", `:=`(c("cds_start", "cds_stop"),
            list(i.l_utr5 + 1, i.l_utr5 + i.l_cds))]
        data[cds_start == 1 & cds_stop == 0, `:=`(cds_start, 0)]

        lev <- sort(unique(data$length))
        data[, `:=`(site_dist_end5, end5 - cds_start)]
        data[, `:=`(site_dist_end3, end3 - cds_start)]
        site_sub <- rbind(site_sub, data[site_dist_end5 <= -6 &
            site_dist_end3 >= 5])
    }
    site_sub <- site_sub[length == end3 - end5 + 1]
    offsets <- .offsetFun(site_sub=site_sub, lev=lev)
    return(offsets)
}

.tempOff <- function(v_dist) {
    ttable <- sort(table(v_dist), decreasing=TRUE)
    ttable_sr <- ttable[as.character(as.numeric(names(ttable)) + 1)]
    ttable_sl <- ttable[as.character(as.numeric(names(ttable)) - 1)]
    tsel <- rowSums(cbind(ttable > ttable_sr, ttable > ttable_sl), na.rm=TRUE)
    return(as.numeric(names(tsel[tsel == 2][1])))
}

.adjOff <- function(dtsite, bestoff) {
    t <- table(factor(dtsite, levels=seq(min(dtsite) - 2, max(dtsite) + 1)))
    t[c(1, 2)] <- t[3] + 1
    locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) ==
        -2)]))) + 1
    adjoff <- locmax[which.min(abs(locmax - bestoff))]
    ifelse(length(adjoff) != 0, adjoff, bestoff)
}

.offsetFun <- function(site_sub, lev) {
    onsite <- count <- .SD <- site_dist_end5 <- NULL
    minlen <- min(site_sub$length)
    maxlen <- max(site_sub$length)
    t <- table(factor(site_sub$length, levels=lev))
    # offset
    offsets <- data.table(length=as.numeric(as.character(names(t))),
        count=as.vector(t))
    offsets[, `:=`(onsite, TRUE)][count == 0, `:=`(onsite, FALSE)]
    offset5 <- site_sub[, list(offset=.tempOff(.SD$site_dist_end5)), by=length]
    offsets.tmp <- merge(offsets, offset5, all.x=TRUE, by="length")
    best.offset <- as.numeric(offsets.tmp[!is.na(offset),
        list(count=sum(count)), by=offset][count == max(count)][, offset])
    # adjusted offset
    adj.tab <- site_sub[, list(offset=.adjOff(site_dist_end5, best.offset)),
        by=length]
    offsets <- merge(offsets, adj.tab, all.x=TRUE, by="length")
    offsets[is.na(offset), `:=`(offset, best.offset)][,
        `:=`(offset, abs(offset))][, `:=`(count, NULL)][, `:=`(onsite, NULL)]
    setnames(offsets, c("qwidth", "psite"))
}

.pMappingAll <- function(bam_file_list, psite.mapping, txdb, all_genes, cores){
    chrom <- seqnames(seqinfo(all_genes))  # extract all chromosome names
    for (k in seq_along(bam_file_list)) {
        if (!file.exists(gsub(".bam", ".bai", bam_file_list[k]))) {
            message("Indexing", bam_file_list[k], "\n")
            Rsamtools::indexBam(bam_file_list[k])
        }
    }
    if (psite.mapping[1] == "auto") {
        message("Computing P-site offsets \n")
        psite.mapping <- .psiteOffset(bam_file_list=bam_file_list, txdb=txdb)
    }
    if (is.null(cores)) {
        cores <- detectCores(logical=FALSE)
    }
    for (k in seq_along(bam_file_list)) {
        message("Computing coverage for", bam_file_list[k], "...")
        cl <- parallel::makeCluster(cores - 1)
        if (cores > 1) {
            registerDoParallel(cl)
        }
        i <- NULL
        results <- foreach(i=seq_along(chrom), .export=c(".pMappingChr",
        ".intersectionStrict", "psiteCal", ".createAnno", ".psiteOffset"),
        .packages=c("Rsamtools", "GenomicAlignments", "GenomicFeatures",
        "data.table", "Rcpp"), .multicombine=TRUE) %dopar%
            {
                return(.pMappingChr(bam=bam_file_list[k], all_genes=all_genes,
                    ch=chrom[i], psite.mapping=psite.mapping))
            }
        stopImplicitCluster()
        stopCluster(cl)
        read_counts <- do.call(rbind, lapply(results, `[[`, 2)) # separate read
        coverage_RLE <- lapply(results, `[[`, 1)  # separate coverage
        coverage_RLE <- coverage_RLE[lengths(coverage_RLE) > 0] # Rle to vector
        coverage_RLE <- do.call("c", coverage_RLE)
        coverage_vector <- lapply(coverage_RLE, as.vector)
        ifelse(k == 1, {
            coverage_matrix <- coverage_vector
            read_counts_matrix <- read_counts
        }, {
            coverage_matrix <- mapply(rbind, coverage_matrix, coverage_vector)
            read_counts_matrix <- cbind(read_counts_matrix, read_counts)
        })
        message("done! \n")
    }
    return(list(coverage=coverage_matrix, counts=read_counts_matrix,
        psite.mapping=psite.mapping))
}

.pMappingChr <- function(bam, all_genes, ch, psite.mapping) {
    ## function to calculate the coverage score and return in Rle objects
    psite <- center <- NULL
    bf <- BamFile(bam)
    tx_compat_cvg <- RleList()
    read_counts <- data.frame()
    genes <- all_genes[vapply(runValue(seqnames(all_genes)),
        function(x) as.character(unique(x)), character(1)) == ch]
    if (!length(genes)) return(list(tx_compat_cvg, read_counts))
    gr <- as(seqinfo(bf), "GRanges")
    if (ch %in% runValue(seqnames(gr))) {
        param <- ScanBamParam(which=gr[ch],
            flag=scanBamFlag(isSecondaryAlignment=FALSE))
        gal2 <- readGAlignments(bf, param=param)
        gal2 <- gal2[strand(gal2) != "*"]
        gal2 <- gal2[is.na(seqnames(gal2)) == FALSE]
        if (length(gal2)) {
            gal3 <- as.data.table(gal2)
            ifelse(all(psite.mapping[1] == "center"),
                gal3[, `:=`(psite, floor(qwidth/2))],
                gal3 <- merge.data.table(gal3, psite.mapping, by="qwidth"))
            gal3[strand == "-", `:=`(psite, qwidth - psite - 1)][, `:=`(center,
                psiteCal(cigar, start, psite))]
            gal3 <- gal3[center > 0, ]
            gal3[, `:=`(qwidth, 1)][, `:=`(width, 1)][, `:=`(njunc, 0)][,
                `:=`(cigar, "1M")]
            gal3[, `:=`(start, center)][, `:=`(end, center)][,
                `:=`(psite, NULL)][, `:=`(center, NULL)]
            gal3 <- makeGRangesFromDataFrame(gal3)
            gal3 <- as(gal3, "GAlignments")
        } else {
            gal3 <- gal2
        }
        results <- .intersectionStrict(genes, gal3)
        gal3 <- gal3[queryHits(results)]
        results <- .intersectionStrict(genes, gal3)
        tx2reads <- setNames(as(t(results), "List"), names(genes))
        compat_reads_by_tx <- extractList(gal3, tx2reads)
        read_counts <- data.frame(lengths(tx2reads))
    } else {
        tx2reads1 <- IntegerList(vector("list", length(genes)))
        names(tx2reads1) <- names(genes)
        compat_reads_by_tx <- extractList(gr, tx2reads1)
        read_counts <- data.frame(lengths(tx2reads1))
    }
    tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx, genes,
        ignore.strand=FALSE)
    return(list(tx_compat_cvg, read_counts))
}

.intersectionStrict <- function(features, reads, ignore.strand=FALSE,
    inter.feature=TRUE){
    ## function adapted to compile compatible reads
    ov <- findOverlaps(reads, features, type="within",
        ignore.strand=ignore.strand)
    if (inter.feature) {
        reads_to_keep <- which(countQueryHits(ov) == 1L)
        ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    return(ov)
}

.segCal <- function(x, ints) {
    l <- sum(x >= ints[, 1])
    cumints <- c(0, cumsum(ints[, 2]))
    return(cumints[l] + x - ints[l, 1] + 1)
}

.totalTransciptCoordinate <- function(x) {
    range.gene <- cbind(start(x), end(x))
    strands <- as.character(strand(x))[1]
    starts <- start(x)
    ends <- end(x)
    if(strands=="+"){
        range.gene2 <- IRanges(starts, ends)
        range.gene2 <- as.matrix(BiocGenerics::union(range.gene2, range.gene2))
        range.relative <- matrix(vapply(range.gene,
            function(x) .segCal(x, range.gene2), numeric(1)), ncol=2)
        range.relative <- cbind(range.relative, range.gene)
        rownames(range.relative) <- unlist(x$EXONNAME) 
    }else{
        m <- max(start(x),end(x))+1
        range.gene1 <- cbind(rev(m-ends),rev(m-starts))
        range.gene2 <- IRanges(range.gene1[,1],range.gene1[,2])
        range.gene2 <- as.matrix(BiocGenerics::union(range.gene2, range.gene2))
        range.relative <- matrix(vapply(range.gene1,
            function(x) .segCal(x, range.gene2), numeric(1)), ncol=2)
        range.relative <- cbind(range.relative,
            rev(range.gene[,2]), rev(range.gene[,1]))
        rownames(range.relative) <- rev(unlist(x$EXONNAME))
    }
    colnames(range.relative) <- c("start_tranx", "end_tranx","start_genome",
        "end_genome")
    range.relative <- as.data.frame(range.relative)
    range.relative$seqnames <- as.character(seqnames(x))
    range.relative$strand <- strands
    range.relative$txname <- unlist(x$TXNAME)
    return(range.relative)
}

#### data binning ####

dataBinning <- function(data, bin.width=0, zero.omit=FALSE, bin.from.5UTR=TRUE,
    cores=NULL) {
    options(warn=-1)
    if (is.null(cores)) cores <- detectCores(logical=FALSE)
    cl <- makeCluster(cores - 1)
    if (cores > 1) registerDoParallel(cl)
    genes.list <- names(data)
    gene <- NULL
    message("Data binning ...")
    DATA.BIN <- foreach(gene=genes.list, .final=function(x) setNames(x,
        genes.list), .packages="reldist", 
        .export=c(".app2Exact",".binning"))  %dopar% {
            data1 <- data[[gene]]
            .binning(data1, bin.width, zero.omit, bin.from.5UTR)
    }
    stopImplicitCluster()
    stopCluster(cl)
    message("done! \n")
    return(DATA.BIN)
}

.binning <- function(data1, bin.width, zero.omit, bin.from.5UTR){
    n <- nrow(data1)
    p <- ncol(data1)
    p.codon <- ceiling(p/3)  # merge every 3 nt into a codon
    ifelse(bin.from.5UTR, bin.codon <- (seq_len(p) - 1)%/%3,
        bin.codon <- rev((rev(seq_len(p)) - 1)%/%3))
    data.codon <- t(apply(data1, 1, function(x) 
        unname(tapply(x, bin.codon, sum))))
    names.codon <- unname(tapply(colnames(data1), bin.codon, function(x)
        paste(x[1],tail(x,n=1), sep="-")))
    colnames(data.codon) <- names.codon
    if (bin.width == 1) return(data.codon) # return codon level
    if (sum(data1)) {
        if (bin.width) {
            p.bin <- ceiling(p.codon/bin.width)  # fixed width
        } else {
            n.data <- mean(rowSums(data.codon))
            iqr.data <- abs(wtd.iqr(x=seq_len(p.codon),
                weight=colMeans(data.codon)))
            n.bin <- ceiling(p.codon/(2 * iqr.data/(n.data^(1/3))))
            p.bin <- min(p.codon, n.bin)  # adaptive width
        }
        bin.size <- .app2Exact(p.codon, rep(1/p.bin, p.bin))
        if (!bin.from.5UTR) bin.size <- rev(bin.size)
        Bin.size <- c(0, cumsum(bin.size))
        data.bin <- do.call("cbind",lapply(split(seq_len(p.codon),
            cut(seq_len(p.codon), Bin.size)), function(x) {
            rowSums(data.codon[, x, drop=FALSE])}))
        bin.size2 <- do.call(c,lapply(split(unname(table((
            seq_len(p) - 1)%/%3)),cut(seq_len(p.codon), Bin.size)),sum))
        Bin.size2 <- c(0, cumsum(unname(bin.size2)))
        names.bin <- lapply(split(seq_len(p),cut(seq_len(p), Bin.size2)),
            function(x) paste(colnames(data1)[x[1]],
            colnames(data1)[tail(x,n=1)], sep="-"))
        colnames(data.bin) <- unname(names.bin)
    } else {
        data.bin <- matrix(0, nrow=n, ncol=1) # no reads
        colnames(data.bin) <- paste(colnames(data1)[1],
            tail(colnames(data1),n=1), sep="-")
    }
    if (zero.omit) data.bin <- data.bin[,colSums(data.bin) > 0, drop=FALSE]
    if (is.vector(data.bin)){
        data.bin <- matrix(data.bin, ncol=1) 
        colnames(data.bin) <- paste(colnames(data1)[1],
            tail(colnames(data1),n=1), sep="-")
    } 
    return(data.bin)
}
    
.app2Exact <- function(n, p) {
    ## p=p[1:m] is an allocation, p[i] >=0, sum(p)=1
    ans <- floor(n * p)
    k <- n - sum(ans)
    if (k > 0) {
        stemp <- n * p - ans
        otemp <- order(-stemp)[seq_len(k)]
        ans[otemp] <- ans[otemp] + 1
    }
    return(ans)
}

#### main function ####
diffPatternTest <- function(data, classlabel, method=c("gtxr", "qvalue")) {
    options(warn=-1)

    ## prepare data
    genes.list <- names(data)
    ## remove other replicates not involved in test
    data <- lapply(data, "[", which(classlabel$comparison != 0), )
    ## remove genes with only a vector
    noread <- which(do.call("c", lapply(data, is.vector)))
    genes.list <- genes.list[which(do.call("c", lapply(data, is.matrix)))]
    data <- data[genes.list]
    genes.list <- genes.list[which(apply(do.call("rbind",
        lapply(data, rowSums)), 1, min) > 0)]
    ## remove genes with no read
    noread2 <- which(apply(do.call("rbind", lapply(data, rowSums)),1,min) == 0)
    data <- data[genes.list]

    result.dp <- .diffPatternTest2(data=data, classlabel=classlabel[classlabel$
        comparison != 0, ], method=method)
    result.dp$small <- names(c(noread, noread2))
    result.dp$data <- data
    result.dp$classlabel <- classlabel
    return(result.dp)
}

diffPatternTestExon <- function(psitemap, classlabel, method=c("gtxr", 
    "qvalue")) {
    options(warn=-1)

    ## prepare data
    data <- psitemap$coverage
    genes.list <- names(data)
    ## remove other replicates not involved in test
    data <- lapply(data, "[", which(classlabel$comparison != 0), )
    genes.list <- genes.list[which(apply(do.call("rbind",
        lapply(data, rowSums)), 1, min) > 0)]
    ## remove genes with no read
    noread <- which(apply(do.call("rbind", lapply(data, rowSums)),1,min) == 0)
    data <- data[genes.list]
    sgtf <- psitemap$exons[genes.list]

    result.dp <- .diffPatternTestExon2(sgtf=sgtf, data=data, classlabel=
        classlabel[classlabel$comparison != 0, ], method=method)
    result.dp$small <- names(noread)
    result.dp$data <- data
    result.dp$classlabel <- classlabel
    tmp <- mapply(base::cbind,result.dp$bin,sgtf, SIMPLIFY = FALSE)
    result.dp$bin <- tmp
    return(result.dp)
}

.diffPatternTest2 <- function(data, classlabel, method) {
    genes.list <- names(data)
    condition <- classlabel$comparison

    tvalue <- 1 - do.call("c", lapply(data, function(x) {
        abs(sum(.svdSelf(x[which(condition == 1), ]) *
            .svdSelf(x[which(condition == 2), ])))
    }))

    message("Bin level testing ...\n")
    DATA.com <- t(do.call("cbind", data))  # link all bins
    classlabel$comparison <- as.factor(classlabel$comparison)
    colnames(DATA.com) <- classlabel$type
    dds <- DESeqDataSetFromMatrix(countData=DATA.com, colData=classlabel,
        design=~comparison)
    S.hat <- lapply(data, normFactor, condition=condition)  # size factor
    normFactors <- t(do.call("cbind", mapply(function(x, y)
        replicate(ncol(y), x), S.hat, data)))
    normalizationFactors(dds) <- normFactors  # assign size factors
    dds <- DESeq(dds)
    res <- as.matrix(results(dds))[, c("pvalue", "log2FoldChange")]

    message("done!\nSplitting test results by genes ...")
    LL <- unlist(lapply(data, ncol))
    res.bin <- apply(res, 2, function(x) split(x, rep(names(LL), LL)))
    res.bin <- mapply(cbind, res.bin$pvalue, res.bin$log2FoldChange)
    res.bin <- lapply(res.bin, function(x) {
        x <- cbind(x, padj=NA)
        colnames(x) <- c("pvalue", "log2FoldChange", method[1])
        x[which(!is.na(x[, "pvalue"])), method[1]] <-
            elitism::p.adjust(na.omit(x[, "pvalue"]), method=method[1])
        return(x)
    })
    pvalue.gene <- do.call("c", lapply(res.bin, function(x) min(x[, method[1]],
        na.rm=TRUE)))
    pvalue.gene <- sort(pvalue.gene)
    genes.list <- names(pvalue.gene)
    message("done!\nGenome wide family control ...")
    ifelse(method[2] == "qvalue", padj.gene <- qvalue(pvalue.gene)$qvalues,
        padj.gene <- elitism::p.adjust(pvalue.gene, method=method[2]))
    gene.result <- data.frame(tvalue=tvalue[genes.list],
        pvalue=pvalue.gene[genes.list], padj=padj.gene[genes.list])
    colnames(gene.result)[3] <- method[2]
    message("done! \n")
    return(list(bin=res.bin[genes.list], gene=gene.result, method=method))
}

.diffPatternTestExon2 <- function(sgtf, data, classlabel, method) {
    genes.list <- names(data)
    condition <- classlabel$comparison
    message("Exon binning ...\n")
    data.exon <- mapply(function(x, y) {
    apply(y, 1, function(z) rowSums2(x, cols=seq(z[1], z[2])))}, data, sgtf)
    tvalue <- 1 - do.call("c", lapply(data.exon, function(x) {
        ifelse(ncol(x) == 1, 1, abs(sum(.svdSelf(x[which(condition == 1), ]) *
            .svdSelf(x[which(condition == 2), ]))))
    }))
    tvalue <- pmax(0, tvalue)
    message("done!\nExon level testing ...\n")
    factor.exon <- lapply(data, function(x) {
        x.codon <- t(apply(x, 1, function(v) {
            unname(tapply(v, (seq_along(v) - 1)%/%3, sum))}))
        normFactor(x.codon, condition)
    })
    normFactors <- t(do.call("cbind", mapply(function(x, y)
        replicate(ncol(y), x), factor.exon, data.exon)))
    DATA.com <- t(do.call("cbind", data.exon))
    classlabel <- classlabel[classlabel$comparison != 0, ]
    classlabel$comparison <- as.factor(classlabel$comparison)
    colnames(DATA.com) <- classlabel$type
    dds <- DESeqDataSetFromMatrix(countData=DATA.com, colData=classlabel,
        design=~comparison)
    normalizationFactors(dds) <- normFactors  # assign size factors
    dds <- DESeq(dds)
    res <- as.matrix(results(dds))[, c("pvalue", "log2FoldChange")]
    message("Exon level testing done!\nSplitting test results by genes ...")
    LL <- unlist(lapply(data.exon, ncol))
    res.bin <- apply(res, 2, function(x) split(x, rep(names(LL), LL)))
    res.bin <- mapply(cbind, res.bin$pvalue, res.bin$log2FoldChange)
    res.bin <- lapply(res.bin, function(x) {
        x <- cbind(x, padj=NA)
        colnames(x) <- c("pvalue", "log2FoldChange", method[1])
        x[which(!is.na(x[, "pvalue"])), method[1]] <-
            elitism::p.adjust(na.omit(x[, "pvalue"]), method=method[1])
        return(x)
    })
    pvalue.gene <- do.call("c", lapply(res.bin, function(x)
        min(x[, method[1]], na.rm=TRUE)))
    ifelse(method[2] == "qvalue", padj.gene <- qvalue(pvalue.gene)$qvalues,
        padj.gene <- elitism::p.adjust(pvalue.gene,method=method[2]))
    gene.result <- data.frame(tvalue=tvalue[genes.list],
        pvalue=pvalue.gene[genes.list],padj=padj.gene[genes.list])
    colnames(gene.result)[3] <- method[2]
    message("done! \n")
    return(list(bin=res.bin[genes.list], gene=gene.result, method=method,
        exon=data.exon))
}

## abundance estimation function
normFactor <- function(x, condition) {
    s.tmp <- rowSums(x)/median(rowSums(x))
    if (all(s.tmp > 0)) {
        y1 <- colSums2(x, rows=(condition == 1))
        y2 <- colSums2(x, rows=(condition == 2))
        ind1 <- (y1 > 0) & (y2 > 0)
        if (sum(ind1) > 1) {
            x <- x[, ind1]
            y1 <- colSums2(x, rows=(condition == 1))
            y2 <- colSums2(x, rows=(condition == 2))
            x.total <- sum(x)
            tmp <- log(y1) - log(y2)
            tmp2 <- as.numeric(quantile(tmp, probs=c(0.25, 0.75), na.rm=TRUE) +
                c(-1.5, 1.5) * IQR(tmp, na.rm=TRUE))
            ind2 <- tmp >= tmp2[1] & tmp <= tmp2[2]
            if (sum(ind2) > 1) {
                x <- x[, ind2]
                if (sum(x)/x.total > 0.1 & all(rowSums(x) > 0)) {
                    s.tmp <- rowSums(x)/median(rowSums(x))
                }
            }
        }
    }
    s.tmp[s.tmp == 0] <- 1e-08
    s.tmp
}

.svdSelf <- function(x) {
    ifelse(is.vector(x), result <- x/sqrt(sum(x^2)), result <- svd(x)$v[, 1])
    return(result)
}


#### plotting ####

plotTrack <- function(data, genes.list, replicates=NULL, exons=FALSE) {
    nt <- RPF <- NULL
    options(warn=-1)
    if (is.null(replicates)) {
        replicates <- seq_len(ncol(data$counts))
    }
    gplot <- NULL
    for (gene in genes.list) {
        gplot[[gene]] <- NULL
        data0 <- data$coverage[[gene]]
        nn <- nrow(data$exons[[gene]])
        x.tick <- c(data$exons[[gene]][,1],data$exons[[gene]][nn,2])
        x.label <- c(data$exons[[gene]][,3],data$exons[[gene]][nn,4])
        track0 <- NULL
        for (i in replicates) {
            rep.i <- data.frame(nt=seq_len(ncol(data0)), RPF=data0[i, ],
                                replicate=colnames(data$counts)[i])
            track0 <- rbind(track0, rep.i)
        }
        gplot[[gene]][["nt"]] <- ggplot(track0, aes(x=nt, y=RPF)) +
            geom_bar(stat="identity") + theme_minimal() +
            facet_wrap(~replicate, ncol=1) +
            scale_x_continuous(breaks=x.tick, labels =x.label) +
            geom_vline(xintercept=x.tick, alpha=0.1)+
            ggtitle(paste("gene: ", gene, ", ", ncol(data0), " nt", sep="")) +
            theme(plot.title=element_text(hjust=0.5), 
                    axis.text.x = element_text(angle = 45, vjust = 1, 
                                            size = 8, hjust = 1))
        if (exons) {
            gplot[[gene]][["exon"]] <- NULL
            exons.gene <- data$exons[[gene]]
            for (j in seq_len(nrow(exons.gene))) {
                exon <- rownames(exons.gene)[j]
                track2 <- track0[track0$nt >= exons.gene[j, 1] & track0$nt <=
                                    exons.gene[j, 2], ]
                gplot[[gene]][["exon"]][[exon]] <- ggplot(track2,
                            aes(x=nt,y=RPF))+ geom_bar(stat="identity") +
                theme_minimal() + facet_wrap(~replicate,ncol=1) +
                scale_x_continuous(breaks = unname(unlist(exons.gene[j, 1:2])),
                            labels = unname(unlist(exons.gene[j, 3:4]))) +
                ggtitle(paste("gene: ", gene, ", exon: ", exon, ", ",
                        diff(unlist(exons.gene[j,1:2])) + 1, " nt", sep="")) +
                    theme(plot.title=element_text(hjust=0.5))
            }
        }
    }
    return(gplot)
}

plotTest <- function(result, genes.list=NULL, threshold=0.05) {
    bin <- RPF <- condition <- NULL
    method <- result$method
    classlabel <- result$classlabel[result$classlabel$comparison != 0, ]
    if (is.null(genes.list)) {
        genes.list <- rownames(result$gene[which(result$gene[, method[2]] <=
                                                threshold),])
    }
    gplot <- NULL
    for (gene in genes.list) {
        data1 <- result$data[[gene]]
        tmp <- result$bin[[gene]][, method[1]]
        diff.spot <- which(tmp <= threshold)
        track1 <- data.frame()
        for (i in seq_len(nrow(classlabel))) {
            rep.i <- data.frame(bin=seq_len(ncol(data1)), RPF=data1[i, ],
                                replicate=rownames(classlabel)[i],
                                condition=classlabel$condition[i])
            rep.i$condition <- as.character(rep.i$condition)
            rep.i[diff.spot, "condition"] <- "DP spot"
            track1 <- rbind(track1, rep.i)
        }
        rownames(track1) <- NULL
        # name.track <- rownames(result$bin[[gene]])
        track1$replicate <- as.factor(track1$replicate)
        track1$replicate <- factor(track1$replicate,
                    levels(track1$replicate)[order(classlabel$comparison)])
        gplot[[gene]] <- ggplot(track1, aes(x=bin, y=RPF, fill=condition)) +
            geom_bar(stat="identity") + theme_minimal() +facet_wrap(~replicate,
                            ncol=1) + scale_x_continuous(breaks = diff.spot,
                labels = formatC(tmp[diff.spot], format = "e", digits = 2)) +
            ggtitle(paste(gene,", ", method[2], " =", formatC(result$gene[gene,
                method[2]], format="e", digits=2), ", T-value =",
                formatC(result$gene[gene, "tvalue"], format="e", digits=2),
                        ", ", ncol(data1), " bins", sep="")) +
            theme(plot.title=element_text(hjust=0.5), 
                    axis.text.x = element_text(size = 8, 
                                        angle = 90, vjust = 0.5, hjust=1))
    }
    return(gplot)
}

#### track output for genome browser ####
binTrack <- function(data, exon.anno){
    classlabel <- data$classlabel
    track.bin <- vector(mode="list", length=nrow(classlabel))
    names(track.bin) <- rownames(classlabel)
    for(n.rep in seq_len(nrow(classlabel))){
        track <- rownames(classlabel)[n.rep]
        result.track <- GRanges(seqnames = NULL, ranges = NULL,
            score = NULL, genes = NULL, tracks = NULL, bin = NULL)
        for(gene in names(data$data)){
            data.bin <- data$data[[gene]]
            bin.range <- colnames(data.bin)
            range.bin <- do.call(rbind, do.call(c, lapply(bin.range,
                strsplit, split = '-')))
            chr <- exon.anno[[gene]]$seqnames[1]
            strand <- exon.anno[[gene]]$strand[1]
            ifelse(strand=="+", {
                starts <- as.integer(range.bin[,1])
                ends <- as.integer(range.bin[,2])
                exon.tmp <- data.table(start = exon.anno[[gene]][,3],
                    end = exon.anno[[gene]][,4])
            },{
                ends <- as.integer(range.bin[,1])
                starts <- as.integer(range.bin[,2])
                exon.tmp <- data.table(start = exon.anno[[gene]][,4],
                    end = exon.anno[[gene]][,3])
            })
            track.tmp <- data.table(start = starts, end = ends,
                score = data$data[[gene]][n.rep,],
                bin = seq_len(nrow(range.bin)),
                padj = data$bin[[gene]][,3])
            setkey(exon.tmp, start, end)
            overlap.bin <- as.matrix(data.table::foverlaps(track.tmp,
                exon.tmp, nomatch = 0L))
            
            track.gene <- GRanges(seqnames = chr, 
                ranges=IRanges(start = overlap.bin[,3] - 1, 
                end = overlap.bin[,4]), strand = strand,
                score = overlap.bin[,5], genes = gene, tracks = track, 
                bin = as.character(overlap.bin[,6]), padj = overlap.bin[,7])
            exon.int <- GRanges(seqnames = chr, ranges = 
                IRanges(start = overlap.bin[,1] - 1, end = overlap.bin[,2]))
            track.gene <- pintersect(track.gene, exon.int)[,1:5]
            result.track <- c(result.track, track.gene)
        }
        track.bin[[track]] <- result.track
    }
    return(track.bin)
}

bpTrack <- function(data, names.rep = NULL, genes.list = NULL){
    if(is.null(names.rep)){
        names.rep <- colnames(data$counts)
    }
    if(is.null(genes.list)){
        genes.list <- names(data$coverage)
    }
    track.bp <- vector(mode = "list", length = length(names.rep))
    names(track.bp) <- names.rep
    for(n.rep in seq_len(length(names.rep))){
        track <- names.rep[n.rep]
        result.track <- GRanges(seqnames = NULL, ranges = NULL,
            coverahe = NULL, genes = NULL)
        for(gene in genes.list){
            chr <- data$exons[[gene]]$seqnames[1]
            strand <- data$exons[[gene]]$strand[1]
            sites <- as.numeric(colnames(data$coverage[[gene]]))
            coverages <- data$coverage[[gene]][n.rep,]
            sites <- sites[which(coverages>0)]
            coverages <- coverages[which(coverages>0)]
            if(length(sites)){
                track.gene <- GRanges(seqnames = chr, 
                    ranges = IRanges(start = sites - 1, end = sites),
                    strand = strand, coverage = coverages)
                result.track <- c(result.track, track.gene)
            }
        }
        track.bp[[track]] <- result.track
    }
    return(track.bp)
}

exonTrack <- function(data, gene){
    classlabel <- data$classlabel
    track.bin <- vector(mode = "list", length = nrow(classlabel))
    names(track.bin) <- rownames(classlabel)
    data.gene <- data$bin[[gene]]
    txnames <- unique(data.gene[,"txname"])
    strands <- data.gene[1, "strand"]
    chr <- data.gene[1, "seqnames"]
    for(n.rep in seq_len(nrow(classlabel))){
        name.track <- rownames(classlabel)[n.rep]
        track.bin[[name.track]] <- vector(mode = "list",
            length = length(txnames))
        names(track.bin[[name.track]]) <- txnames
        for(tx in txnames){
            data.tx <- data.gene[data.gene[,"txname"] == tx,]
            ifelse(strands == "+", {
                starts <- data.tx[,"start_genome"]
                ends <- data.tx[,"end_genome"]
            },{
                ends <- data.tx[,"start_genome"]
                starts <- data.tx[,"end_genome"]
            })
            track.tx <- GRanges(seqnames = chr, 
                ranges = IRanges(start = starts - 1, end = ends - 1),
                strand = strands,
                score = data$exon[[gene]][n.rep, rownames(data.tx)], 
                genes = gene, tracks = name.track, 
                exon = rownames(data.tx))
            track.bin[[name.track]][[tx]] <- track.tx
        }
    }
    return(track.bin)
}
