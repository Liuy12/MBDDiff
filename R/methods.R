# library(dplyr)
# library(RMySQL)
# library(pryr)
# library(doParallel)
# library(data.table)
# library(pracma)
# library(rgl)
# library(rCharts)
# library(plotly)
# library(slidify)
# library(slidifyLibraries)

## show a list of available organisms from UCSC database
ShowAvailableOrg <- function(){
  return(data(Organisms))
}

GetPromoterAnno <- function(organism, save = T, Dir = NULL){
  if(!organism %in% ShowAvailableOrg()[,2])
    stop('Organism not found. Please call funcition ShowAvailabelOrg() to retrive a list of available organisms from UCSC. Remember to use accession code rather than the organism name directly.')
  ucsc_db <- src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
  Ref_Flat <- tbl(ucsc_db, sql("SELECT geneName, chrom, strand, txStart, txEnd FROM refFlat"))
  Promoter_Anno <- GetPromoters(as.data.frame(Ref_Flat))
  if(save){
    if(is.null(Dir))
      stop('please provide a path where annotation file will be saved')
    write.table(Promoter_Anno, paste(Dir, '/Promoter_Anno.bed', sep = ''), quote = F, sep = '\t', row.names = F)
  }
  return(Promoter_Anno)
}

bedtoolsr <- function(bam = NULL, bamdir = NULL, bed){
  if(is.null(bamdir)){
    if(is.null(bam) | !file.exists(bam))
      stop('bam file not provided or can not be found.')
    else
      bamfiles <- bam
  }
  else
    bamfiles <- grep('.bam', dir(bamdir), value = T)
  if(length(bamfiles))
    stop('There are no bam files in the directory')
  else{
    bamfiles <- paste(bamdir, '/', bamfiles, sep ='')
    temp1 <- c()
    for(i in length(bamfiles)){
      command <- paste('coverageBed -abam', bamfiles[i], '-b', bed, sep = ' ')
      cat('processing: ', bamfiles[i], sep = '\t')
      temp <- fread(input = command, data.table = F)
      temp <- arrange(temp, V1, V2, V3)
      temp1 <- bind_cols(temp1, temp$V6)
    }
    rownames(temp1) <- temp$V4
    return(list(
      count = temp1,
      gaplength = temp[,(ncol(temp)-1)]
      ))
  }
}


GetPromoters <- function(uctable, upstream = 2000, downstream = 2000){
  uctable <- uctable %>% 
    filter(sapply(strsplit(chrom, '_'), length) == 1)
  uniqid <- uctable %>%
    select(geneName) %>%
    distinct(geneName)
  ### A single gene might encode transcript on different chromsome, 
  ### strand or transcriptional start site
  uctable1 <- data.frame()
  for(i in 1:length(uniqid$geneName)){
    cat(i, '\n')
    index <- grep(paste('^', uniqid$geneName[i], '$', sep=''), uctable$geneName)
    temp <- unique(uctable$chrom[index])
    temp1 <- unique(uctable$strand[index])
    for(j in 1:length(temp)){
      for(k in 1:length(temp1)){
        index1 <- which(uctable[index,2] == temp[j] & uctable[index, 3] == temp1[k])
        if(length(index1)){
          if(temp1[k] == '+')
            temp2 <- unique(uctable[index[index1],4])
          else
            temp2 <- unique(uctable[index[index1],5])
          temp3 <- data.frame(Symbol = paste(uniqid$geneName[i], '#', 1:length(temp2), sep=''),
                              Chrom = rep(temp[j], length(temp2)),
                              Strand = rep(temp1[k], length(temp2)),
                              TSS = temp2, stringsAsFactors = F)
          uctable1 <- rbind(uctable1, temp3)
        }
      }
    }
  }
  Promoter4k <- data.frame(chrom = uctable1$Chrom,
                           chromStart = ifelse(uctable1$Strand == '+', uctable1$TSS - upstream, uctable1$TSS - downstream),
                           chromEnd = ifelse(uctable1$Strand == '+', uctable1$TSS + downstream, uctable1$TSS + upstream),
                           name = uctable1$Symbol,
                           Strand = uctable1$Strand
  )
  return(Promoter4k)
}

## retrive chromsome length information from UCSC database.
GetChromLength <- function(organism){
  if(!organism %in% ShowAvailableOrg()[,2])
    stop('Organism not found. Please call funcition ShowAvailabelOrg() to retrive a list of available organisms from UCSC. Remember to use accession code rather than the organism name directly.')
  ucsc_db <- src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
  Ref_Flat <- tbl(ucsc_db, sql("SELECT chrom, size FROM chromInfo"))
}

#### create window function to create 100bp window across whole genome
CreateWindows <- function(chromlen, binsize, cores = detectCores()){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  Wholegenome_bin <- foreach (i=1:length(chromlen), .options.multicore=list(preschedule=FALSE), 
                                    .combine=bind_rows, .inorder=TRUE, .verbose=TRUE, 
                                    .errorhandling='stop', .multicombine=TRUE) %dopar% {
                                      cat(i, '\t')
                                      bin <- seq(0, chromlen[i], by = binsize)
                                      temp <- data.frame(chrom = rep(names(chromlen)[i], length(bin)),
                                                         chromStart = bin,
                                                         chromEnd = c(bin[2:length(bin)], chromlen[i]),
                                                         name = paste(names(chromlen)[i], '_', 1:length(bin), sep ='')
                                      )
                                      temp
                                    }
  return(Wholegenome_bin)
}

########### exclude regions: 1. in promoter regions, 2. in CpG island regions 3. ATCG content is not 100%
Getfasta <- function(fa, bed){
  if(is.null(fa) | !file.exists(fa))
    stop('Sequence file not provided or cannot be found')
  if(is.null(bed) | !file.exists(bed))
    stop('Bed file not provided or cannot be found')
  command <- paste('bedtools getfasta -fi', fa, '-bed', bed, '-fo', gsub('.bed', '.fa', bed), sep = ' ')
  system(command)
}

ConstructExRegions <- function(organism, promoter_anno, CpG = T){
  if(CpG){
    if(!organism %in% ShowAvailableOrg()[,2])
      stop('Organism not found. Please call funcition ShowAvailabelOrg() to retrive a list of available organisms from UCSC. Remember to use accession code rather than the organism name directly.')
    ucsc_db <- src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
    CpGisland <- as.data.frame(tbl(ucsc_db, sql("SELECT chrom, chromStart, chromEnd FROM cpgIslandExtUnmasked")))
    colnames(promoter_anno) <- colnames(CpGisland) <- c('chrom', 'start', 'end')
    bind_rows(promoter_anno, as.data.frame(CpGisland))
  }
  else
    promoter_anno
}

FilterRegions <- function(bed, fasta, Exbed, cores = detectCores()){
  Exclude_index <- ExcludeIntersection(bed, Exbed, cores = cores)
  letter_freq <- CountFreqency(fasta)
  Exclude_index1 <- which(letter_freq$ATCG != 1)
  Exclude_index_all <- Reduce(dplyr::union, list(Exclude_index, Exclude_index2, Exclude_index3))
  return(bind_cols(bed[-Exclude_index_all,], letter_freq$CG[-Exclude_index_all]))
}

CountFreqency <- function(fasta, CG = T, ATCG = T){
  fasta_file <- readDNAStringSet(fasta,use.names = T)
  Freqmat <- data.frame()
  if(CG)
    Freqmat$CG <- letterFrequency(fasta_file, letters = 'CG', as.prob = T)
  if(ATCG)
    Freqmat$ATCG <- letterFrequency(fasta_file, letters = 'ATCG', as.prob = T)
  return(Freqmat)
}

########## Exclude intersection between regions
ExcludeIntersection <- function(data1, data2, cores = detectCores()){
  if(ncol(data1) != ncol(data2))
    stop('The two datasets must have equal column size.')
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  foreach (i=1:nrow(data2), .options.multicore=list(preschedule=FALSE), 
           .combine=c, .inorder=TRUE, .verbose=TRUE, 
           .errorhandling='stop', .multicombine=TRUE) %dopar% {
             cat(i, '\t')
             index <- which(data1[,1] == data2[i,1])
             index1 <- which(data1[index,2] > data2[i,2] & data1[index,3] < data2[i,3])
             index2 <- which(data1[index,2] < data2[i,2] & data1[index,3] > data2[i,2])
             index3 <- which(data1[index,2] < data2[i,3] & data1[index,3] > data2[i,3])
             index4 <- index[c(index1, index2, index3)]
             if(length(index4))
               index4
           }
}

######### identify background region for each promoter and calculate background counts for each promoter 
### for each gene, get 40 100bp windows that have lowest average rpkm, then add the count
### together as the background estimation for this gene 
### all the genes have at least 4k up or down stream
### cannot find enough windows in -4k + 4k windows
### will try to find closest 80 windows (GC < 0.4), and choose 40 among them based on methylation levels 

IdentifyBackground <- function(organism, bed_path, binsize, promo_bed, cores = detectCores()){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  chromlen <- GetChromLength(organism)
  bed_100bp <- CreateWindows(chromlen, binsize)
  write.table(bed_100bp, paste(bed_path, '/bed_100bp.bed', sep = ''), quote = F, sep = '\t', row.names = F)
  Exclude_regions <- ConstructExRegions(organism, promo_bed)
  Getfasta(fa, bed_path)
  bed_100bp_filtered <- FilterRegions(bed_100bp, gsub('.bed', '.fa', bed_path), Exclude_regions) 
  Background_region <- foreach (i=1:nrow(promo_bed), .options.multicore=list(preschedule=FALSE), 
                                .combine=bind_rows, .inorder=TRUE, .verbose=TRUE, 
                                .errorhandling='stop', .multicombine=TRUE) %dopar% {
                                  cat('promoter', i, '\t')
                                  preset <- filter(bed, Chrom == promo_bed[i,1])
                                  preset <- mutate(preset, proximity = abs(Start - promo_bed[i,7]))
                                  preset <- arrange(preset, proximity)
                                  k <- 0 
                                  bgregion <- c()
                                  for(j in 1:nrow(preset)){
                                    if(preset[j,5] < 0.4){
                                      k <- k + 1
                                      bgregion <- bind_rows(bgregion, preset[j,])
                                    }
                                    if(k == 80)
                                      break
                                  }
                                  bgregion
                                }
}

CalculateTPM <- function(dataMat, gaplength){
  if(nrow(dataMat) != length(gaplength))
    stop('The row dimension of the dataset has to equal to the length of gap length. Remember to use accession code rather than the organism name directly.')
  libsize <- apply(dataMat, 2, sum)
  gaplength <- t(pracma::repmat(gaplength, ncol(dataMat), 1))
  dataMat_TPM <- (dataMat/gaplength)*10^6
  return(apply(dataMat_TPM, 1, mean))
}

Summerizebg <- function(dataMat_bg, gaplength){
  dataMat_bg$TPM <- CalculateTPM(dataMat_bg, gaplength)
  dataMat_bg <- dataMat_bg %>% 
    group_by(name) %>% 
    arrange(TPM) %>% 
    dplyr::slice(1:40) %>% 
    summarise_each(funs(sum)) %>% 
    arrange(name) %>% 
    select(-1)
}

MBDDiff <- function(promoter, background, conditions, method = "pooled", 
                    sharingMode = "maximum", fitType = "local", pvals_only = FALSE, paraMethod='NP'){
  MBD <- XBSeqDataSet(promoter, background, conditions)
  MBD <- estimateRealCount(MBD)
  MBD <- estimateSizeFactors(MBD)
  MBD <- estimateSCV(MBD, method=method, sharingMode=sharingMode, fitType=fitType)
  Teststas <- XBSeqTest(MBD, levels(conditions)[1L], levels(conditions)[2L], method =paraMethod)
  return(list(MBD, Teststas))
}


######### GC enrichment test by yidong's matlab program 