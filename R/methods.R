## show a list of available organisms from UCSC database
ShowAvailableOrg <- function() {
  data(Organisms)
}

GetPromoterAnno <- function(organism, save = FALSE, Dir = NULL) {
  ShowAvailableOrg()
  if (!organism %in% Organisms[,2])
    stop(
      'Organism not found. Please call funcition ShowAvailabelOrg() to retrive a list of available organisms from UCSC. Remember to use accession code rather than the organism name directly.'
    )
  ucsc_db <-
    src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
  Ref_Flat <-
    tbl(ucsc_db, sql("SELECT geneName, chrom, strand, txStart, txEnd FROM refFlat"))
  Promoter_Anno <- GetPromoters(as.data.frame(Ref_Flat))
  if (save) {
    if (is.null(Dir))
      stop('please provide a path where annotation file will be saved')
    write.table(
      Promoter_Anno, paste(Dir, '/Promoter_Anno.bed', sep = ''), quote = F, sep = '\t', row.names = F
    )
  }
  return(Promoter_Anno)
}

bedtoolsr <- function(bam = NULL, bamdir = NULL, bed) {
  if (is.null(bamdir)) {
    if (is.null(bam) | !file.exists(bam))
      stop('bam file not provided or can not be found.')
    else
      bamfiles <- bam
  }
  else
    bamfiles <- grep('.bam', dir(bamdir), value = TRUE)
  if (length(bamfiles))
    stop('There are no bam files in the directory')
  else{
    bamfiles <- paste(bamdir, '/', bamfiles, sep = '')
    temp1 <- c()
    for (i in length(bamfiles)) {
      command <-
        paste('coverageBed -abam', bamfiles[i], '-b', bed, sep = ' ')
      cat('processing: ', bamfiles[i], sep = '\t')
      temp <- fread(input = command, data.table = FALSE)
      temp <- arrange(temp, V1, V2, V3)
      temp1 <- bind_cols(temp1, temp$V6)
    }
    rownames(temp1) <- temp$V4
    return(list(count = temp1,
                gaplength = temp[,(ncol(temp) - 1)]))
  }
}


GetPromoters <-
  function(uctable, upstream = 2000, downstream = 2000) {
    uctable <- uctable %>%
      filter(sapply(strsplit(chrom, '_'), length) == 1)
    uniqid <- uctable %>%
      select(geneName) %>%
      distinct(geneName)
    ### A single gene might encode transcript on different chromsome,
    ### strand or transcriptional start site
    uctable1 <- data.frame()
    for (i in 1:length(uniqid$geneName)) {
      cat(i, '\n')
      index <-
        grep(paste('^', uniqid$geneName[i], '$', sep = ''), uctable$geneName)
      temp <- unique(uctable$chrom[index])
      temp1 <- unique(uctable$strand[index])
      for (j in 1:length(temp)) {
        for (k in 1:length(temp1)) {
          index1 <-
            which(uctable[index,2] == temp[j] & uctable[index, 3] == temp1[k])
          if (length(index1)) {
            if (temp1[k] == '+')
              temp2 <- unique(uctable[index[index1],4])
            else
              temp2 <- unique(uctable[index[index1],5])
            temp3 <-
              data.frame(
                Symbol = paste(uniqid$geneName[i], '#', 1:length(temp2), sep = ''),
                Chrom = rep(temp[j], length(temp2)),
                Strand = rep(temp1[k], length(temp2)),
                TSS = temp2, stringsAsFactors = FALSE
              )
            uctable1 <- rbind(uctable1, temp3)
          }
        }
      }
    }
    Promoter4k <- data.frame(
      chrom = uctable1$Chrom,
      chromStart = ifelse(
        uctable1$Strand == '+', uctable1$TSS - upstream, uctable1$TSS - downstream
      ),
      chromEnd = ifelse(
        uctable1$Strand == '+', uctable1$TSS + downstream, uctable1$TSS + upstream
      ),
      name = uctable1$Symbol,
      Strand = uctable1$Strand
    )
    return(Promoter4k)
  }

## retrive chromsome length information from UCSC database.
GetChromLength <- function(organism) {
  cat('Retriving chromsome length information from UCSC database', '\n')
  ShowAvailableOrg()
  if (!organism %in% Organisms[,2])
    stop(
      'Organism not found. Please call funcition ShowAvailabelOrg() to retrive a list of available organisms from UCSC. Remember to use accession code rather than the organism name directly.'
    )
  ucsc_db <-
    src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
  Ref_Flat <- tbl(ucsc_db, sql("SELECT chrom, size FROM chromInfo"))
}

#### create window function to create 100bp window across whole genome
CreateWindows <- function(chromlen, binsize) {
  cat(paste("Generate ", binsize, " bp windows across the whole genome", sep = ''), "\n")
  Wholegenome_bin <-
    foreach (
      i = 1:length(chromlen$size), .options.multicore = list(preschedule = FALSE),
      .combine = bind_rows, .inorder = TRUE, .verbose =
        TRUE,
      .errorhandling = 'stop', .multicombine =
        TRUE
    ) %dopar% {
      cat(i, '\t')
      bin <-
        seq(0, chromlen$size[i], by = binsize)
      temp <-
        data.frame(
          chrom = rep(chromlen$chrom[i], length(bin)),
          chromStart = bin,
          chromEnd = c(bin[2:length(bin)], chromlen$size[i])
        )
      temp
    }
  return(Wholegenome_bin)
}

########### exclude regions: 1. in promoter regions, 2. in CpG island regions 3. ATCG content is not 100%
Getfasta <- function(fa, bed) {
  cat("Building sequence file for specified binsize", "\n")
  if (is.null(fa) | !file.exists(fa))
    stop('Sequence file not provided or cannot be found')
  if (is.null(bed) | !file.exists(bed))
    stop('Bed file not provided or cannot be found')
  command <-
    paste('bedtools getfasta -fi', fa, '-bed', bed, '-fo', gsub('bed$', 'fa', bed), sep = ' ')
  system(command)
}

ConstructExRegions <- function(organism, promoter_anno, CpG = TRUE) {
  cat("Construct regions to exclude when building background noise", "\n")
  promoter_anno <-
    select(promoter_anno, chrom, chromStart, chromEnd)
  if (CpG) {
    ShowAvailableOrg()
    if (!organism %in% Organisms[,2])
      stop(
        'Organism not found. Please call funcition ShowAvailabelOrg() to retrive a list of available organisms from UCSC. Remember to use accession code rather than the organism name directly.'
      )
    ucsc_db <-
      src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
    CpGisland <-
      as.data.frame(tbl(
        ucsc_db, sql(
          "SELECT chrom, chromStart, chromEnd FROM cpgIslandExtUnmasked"
        )
      ))
    colnames(promoter_anno) <-
      colnames(CpGisland) <- c('chrom', 'start', 'end')
    bind_rows(promoter_anno, as.data.frame(CpGisland))
  }
  else
    promoter_anno
}

FilterRegions <- function(bed, fasta, Exbed) {
  cat("Build preliminary background annotation", "\n")
  Exclude_index <- ExcludeIntersection(bed, Exbed)
  rm(Exbed)
  letter_freq <- CountFreqency(fasta)
  Exclude_index1 <- which(letter_freq$ATCG != 1)
  Exclude_index_all <-
    Reduce(dplyr::union, list(Exclude_index, Exclude_index1))
  rm(list = c("Exclude_index", "Exclude_index1"))
  Filter_regions <-
    bind_cols(bed[-Exclude_index_all,], letter_freq[-Exclude_index_all,][1])
  rm(list = c("Exclude_index_all", "bed", "letter_freq"))
  colnames(Filter_regions) <-
    c('chrom', 'chromStart', 'chromEnd', 'GC')
  return(Filter_regions)
}

CountFreqency <- function(fasta, CG = TRUE, ATCG = TRUE) {
  fasta_file <- readDNAStringSet(fasta,use.names = TRUE)
  Freqmat <- data.frame(matrix(NA, nrow = length(fasta_file), ncol = 2))
  colnames(Freqmat) <- c('CG', 'ATCG')
  if (CG)
    Freqmat$CG <-
    letterFrequency(fasta_file, letters = 'CG', as.prob = TRUE)[,1]
  if (ATCG)
    Freqmat$ATCG <-
    letterFrequency(fasta_file, letters = 'ATCG', as.prob = TRUE)[,1]
  return(Freqmat)
}

########## Exclude intersection between regions
ExcludeIntersection <-
  function(data1, data2) {
    if (ncol(data1) != ncol(data2))
      stop('The two datasets must have equal column size.')
    binsize <- data1[1,3] - data1[1,2]
    foreach (
      i = 1:nrow(data2), .options.multicore = list(preschedule = FALSE),
      .combine = c, .inorder = TRUE, .verbose = TRUE,
      .errorhandling = 'stop', .multicombine = TRUE
    ) %dopar% {
      cat(i, '\t')
      index <- which(data1[,1] == data2[i,1])
      index1 <- which(data1[index,2] + binsize > data2[i,2] & data1[index,3] + binsize < data2[i,3])
      if (length(index1))
        index[index1]
    }
  }

######### identify background region for each promoter and calculate background counts for each promoter
### for each gene, get 40 100bp windows that have lowest average rpkm, then add the count
### together as the background estimation for this gene
### all the genes have at least 4k up or down stream
### cannot find enough windows in -4k + 4k windows
### will try to find closest 80 windows (GC < 0.4), and choose 40 among them based on methylation levels

IdentifyBackground <-
  function(organism, bed_path, binsize = 100, promo_bed, fa) {
    chromlen <- as.data.frame(GetChromLength(organism))
    chromlen <- chromlen[which(sapply(strsplit(chromlen$chrom, "_"), length) ==1),]
    bed_100bp <- as.data.frame(CreateWindows(chromlen, binsize))
    rm(chromlen)
    bed_path <- paste(bed_path, 'bed_100bp.bed', sep = '')
    storage.mode(bed_100bp[,2]) <-
      storage.mode(bed_100bp[,3]) <- 'integer'
    cat("Write bed file to disk", "\n")
    write.table(
      bed_100bp, bed_path, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE
    )
    Exclude_regions <-
      as.data.frame(ConstructExRegions(organism, promo_bed))
    Getfasta(fa, bed_path)
    bed_100bp_filtered <-
      FilterRegions(bed_100bp, gsub('bed$', 'fa', bed_path), Exclude_regions)
    rm(Exclude_regions)
    Background_region <-
      foreach (
        i = 1:nrow(promo_bed), .options.multicore = list(preschedule = FALSE),
        .combine = bind_rows, .inorder = TRUE, .verbose =
          TRUE,
        .errorhandling = 'stop', .multicombine =
          TRUE
      ) %dopar% {
        cat('promoter', i, '\t')
        preset <-
          filter(bed_100bp_filtered, chrom == promo_bed$chrom[i])
        preset <-
          mutate(
            preset, proximity = abs(chromStart - promo_bed$chromStart[i]),
            geneName = promo_bed$name[i]
          )
        preset <-
          arrange(preset, proximity)
        k <- 0
        bgregion <- c()
        for (j in 1:nrow(preset)) {
          if (preset[j,4] < 0.4) {
            k <- k + 1
            bgregion <-
              bind_rows(bgregion, preset[j,])
          }
          if (k == 80)
            break
        }
        bgregion
      }
    return(Background_region)
  }

CalculateTPM <- function(dataMat, gaplength) {
  if (nrow(dataMat) != length(gaplength))
    stop(
      'The row dimension of the dataset has to equal to the length of gap length. Remember to use accession code rather than the organism name directly.'
    )
  libsize <- apply(dataMat, 2, sum)
  gaplength <- t(pracma::repmat(gaplength, ncol(dataMat), 1))
  dataMat_TPM <- (dataMat / gaplength) * 10 ^ 6
  return(apply(dataMat_TPM, 1, mean))
}

Summerizebg <- function(dataMat_bg, gaplength) {
  dataMat_bg$TPM <- CalculateTPM(dataMat_bg, gaplength)
  dataMat_bg <- dataMat_bg %>%
    group_by(name) %>%
    arrange(TPM) %>%
    dplyr::slice(1:40) %>%
    summarise_each(funs(sum)) %>%
    arrange(name) %>%
    select(-1)
}

MBDDiff <-
  function(promoter, background, conditions, method = "pooled",
           sharingMode = "maximum", fitType = "local", pvals_only = FALSE, paraMethod =
             'NP') {
    MBD <- XBSeqDataSet(promoter, background, conditions)
    MBD <- estimateRealCount(MBD)
    MBD <- estimateSizeFactors(MBD)
    MBD <-
      estimateSCV(MBD, method = method, sharingMode = sharingMode, fitType =
                    fitType)
    Teststas <-
      XBSeqTest(MBD, levels(conditions)[1L], levels(conditions)[2L], method =
                  paraMethod)
    return(list(MBD, Teststas))
  }


######### GC enrichment test by yidong's matlab program 