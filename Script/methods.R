library(dplyr)
library(RMySQL)
library(pryr)

GetPromoterAnno <- function(organism){
  ucsc_db <- src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
  Ref_Flat <- tbl(ucsc_db, sql("SELECT geneName, chrom, strand, txStart, txEnd FROM refFlat"))
  Promoter_Anno <- GetPromoters(as.data.frame(Ref_Flat))
  return(Promoter_Anno)
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
                                                         chromEnd = c(bin[2:length(bin)], chromlen[i])
                                      )
                                      temp
                                    }
  return(Wholegenome_bin)
}

t2 <- Sys.time()
t2 - t1

########### exclude regions: 1. in promoter regions, 2. in CpG island regions 3. 
ATCGcont <- letterFrequency(Wholegenome_100bp_fa, letters = 'ATCG', as.prob = T)

########## using CpG island program from UCSC
CpGisland_new <- fread('bam file/CpGisland_new.bed')

registerDoMC(40)
t1 <- Sys.time()
Exclude_index <- ExcludeIntersection(Wholegenome_100bp_bed, hg4kpromoter)
t2 <- Sys.time()
t2 - t1

registerDoMC(80)
t1 <- Sys.time()
Exclude_index1 <- ExcludeIntersection(Wholegenome_100bp_bed, CpGisland_new)
t2 <- Sys.time()
t2 - t1

registerDoMC(80)
t1 <- Sys.time()
Exclude_index2 <- ExcludeIntersection(Wholegenome_100bp_bed, CpGisland_unmask)
t2 <- Sys.time()
t2 - t1

hg4kpromoter <- mutate(hg4kpromoter, V7 = V2+2000)
# hg4kpromoter_unique <- hg4kpromoter %>% 
#   group_by(V1) %>%
#   distinct(V7)
Exclude_index3 <- which(Wholegenome_100bp_bed$ATCG != 1)
#Exclude_index_100bp <- dplyr::union(Exclude_index, Exclude_index3)
Exclude_index_100bp <- Reduce(dplyr::union, list(Exclude_index, Exclude_index2, Exclude_index3))

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



######### calculate GC content for each 100bp window

# linux code
#bedtools getfasta -fi ../../genomics/genome.fa -bed Wholegenome_100bp.bed -fo Wholegenome_100bp.fa -name

Wholegenome_100bp_fa <- readDNAStringSet('bam file/Wholegenome_100bp.fa',use.names = T)
GCcont <- letterFrequency(Wholegenome_100bp_fa, letters = 'CG', as.prob = T)


######### identify background region for each promoter and calculate background counts for each promoter 


### for each gene, get 40 100bp windows that have lowest average rpkm, then add the count
### together as the background estimation for this gene 
### all the genes have at least 4k up or down stream
### cannot find enough windows in -4k + 4k windows
### will try to find closest 80 windows (GC < 0.4), and choose 40 among them. 
temp <- bind_cols(as.data.frame(RPKM_combined_100bp[-Exclude_index_100bp, c(1:3, 10)]), as.data.frame(Wholegenome_100bp_bed[-Exclude_index_100bp, 4]))
colnames(temp) <- c(colnames(temp)[1:4], 'GC_content')

registerDoMC(80)
t1 <- Sys.time()
Background_region <- foreach (i=1:nrow(hg4kpromoter), .options.multicore=list(preschedule=FALSE), 
                              .combine=bind_rows, .inorder=TRUE, .verbose=TRUE, 
                              .errorhandling='stop', .multicombine=TRUE) %dopar% {
                                cat('promoter', i, '\t')
                                preset <- filter(temp, Chrom == hg4kpromoter[i,1])
                                preset <- mutate(preset, proximity = abs(Start - hg4kpromoter[i,7]))
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
                                bgregion <- arrange(bgregion, RPKM_average)
                                bgregion <- bgregion %>% dplyr::slice(1:40) %>% mutate(name = hg4kpromoter[i,4])
                                bgregion
                              }
t2 <- Sys.time()
t2 - t1

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
