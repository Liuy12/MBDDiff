library(dplyr)
library(RMySQL)
library(pryr)

GetPromoterAnno <- function(organism, save = F, Dir = NULL){
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

bedtoolsr <- function(bamdir, bed){
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


#### create window function to create 100bp window across whole genome
registerDoMC(80)
t1 <- Sys.time()
CreateWindows <- function(chromlen, binsize){
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

t2 <- Sys.time()
t2 - t1

########### exclude regions: 1. in promoter regions, 2. in CpG island regions 3. ATCG content is not 100%
# linux code # get sequence file for each bin
# bedtools getfasta -fi ../../genomics/genome.fa -bed Wholegenome_100bp.bed -fo Wholegenome_100bp.fa -name

Getfasta <- function(fa, bed, path = NULL){
  if(is.null(fa))
    stop('please provide sequence file')
  if(is.null(path))
    stop('please provide path where the sequence file will be saved')
  command <- paste('bedtools getfasta -fi', 'fa', '-bed', bed, '-fo', gsub('.bed', '.fa', bed), sep = ' ')
  try(system(command))
}

ConstructExRegions <- function(organism, promoter_anno, CpG = T){
  if(CpG){
    ucsc_db <- src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
    CpGisland <- as.data.frame(tbl(ucsc_db, sql("SELECT chrom, chromStart, chromEnd FROM cpgIslandExtUnmasked")))
    colnames(promoter_anno) <- colnames(CpGisland) <- c('chrom', 'start', 'end')
    bind_rows(promoter_anno, as.data.frame(CpGisland))
  }
  else
    promoter_anno
}


FilterRegions <- function(bed, fasta, Exbed){
  registerDoMC(40)
  t1 <- Sys.time()
  Exclude_index <- ExcludeIntersection(bed, Exbed)
  t2 <- Sys.time()
  t2 - t1
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
ExcludeIntersection <- function(data1, data2){
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




######### use bedtools to count reads aligned to each 100bp window

## linux code
# time for i in $(ls | grep .bam)
# do   
# echo $i
# coverageBed -abam $i -b Wholegenome_100bp.bed > $i-wholegenome-100bp.txt
# done


######### identify background region for each promoter and calculate background counts for each promoter 


### for each gene, get 40 100bp windows that have lowest average rpkm, then add the count
### together as the background estimation for this gene 
### all the genes have at least 4k up or down stream
### cannot find enough windows in -4k + 4k windows
### will try to find closest 80 windows (GC < 0.4), and choose 40 among them. 
IdentifyBackground <- function(bed, promo_bed){
  registerDoMC(80)
  t1 <- Sys.time()
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
#                                  bgregion <- arrange(bgregion, RPKM_average)
#                                  bgregion <- bgregion %>% dplyr::slice(1:40) %>% mutate(name = hg4kpromoter[i,4])
                                  bgregion
                                }
  t2 <- Sys.time()
  t2 - t1
}


bgregion <- arrange(bgregion, RPKM_average)
bgregion <- bgregion %>% dplyr::slice(1:40) %>% mutate(name = hg4kpromoter[i,4])


Background_region_aggre <- Background_region %>% group_by(name) %>% summarise(sum(RPKM_average)/40) %>%
  mutate(group = 'Background')
colnames(Background_region_aggre) <- c('name', 'RPKM_average', 'group')





######### GC enrichment test by yidong's matlab program 




#Usage:


Ref_Flat <- GetPromoterAnno('hg19')
Ref_Flat <- as.data.frame(Ref_Flat)

t1 <- Sys.time()

temp <- GetPromoters(Ref_Flat[1:100,])

t2 <- Sys.time()
t2 - t1
