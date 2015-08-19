GetPromoterAnno <- function(organism){
  ucsc_db <- src_mysql(organism, 'genome-mysql.cse.ucsc.edu', user = 'genome')
  Ref_Flat <- tbl(ucsc_db, sql("SELECT geneName, chrom, strand, txStart, txEnd FROM refFlat"))
  Promoter_Anno <- GetPromoters(Ref_Flat)
}

###### parallel computing
library(doParallel)
library(doMC)

cl <- makeCluster(3)
registerDoParallel(cl)

benchmark(a <- foreach(i = 1:3, .combine = c) %dopar% sqrt(4),
          b <- parSapply(cl, list(1,2,3), sqrt),
          c <- bplapply(1:3, sqrt, BPPARAM = DoparParam()),
          d <- foreach(i = 1:3, .combine = c) %do% sqrt(i),
          replications = 1000,
          columns=c('test', 'elapsed', 'replications'))

GetPromoters <- function(uctable){
  uctable <- uctable %>% 
    filter(sapply(strsplit(chrom, '_'), length) == 1)
  uniqid <- uctable %>%
    select(geneName) %>%
    distinct(geneName)
  ### A single gene might encode transcript on different chromsome, 
  ### strand or transcriptional start site
  uctable1 <- data.frame()
  for(i in 1:length(uniqid)){
    cat(i, '\n')
    index <- grep(paste('^', uniqid[i], '$', sep=''), uctable$V1)
    temp <- unique(uctable$V3[index])
    temp1 <- unique(uctable$V4[index])
    for(j in 1:length(temp)){
      for(k in 1:length(temp1)){
        index1 <- which(uctable[index,3] == temp[j] & uctable[index, 4] == temp1[k])
        if(length(index1)){
          if(temp1[k] == '+')
            temp2 <- unique(uctable[index[index1],5])
          else
            temp2 <- unique(uctable[index[index1],6])
          temp3 <- data.frame(Symbol = rep(uniqid[i], length(temp2)),
                              Chrom = rep(temp[j], length(temp2)),
                              Strand = rep(temp1[k], length(temp2)),
                              TSS = temp2, stringsAsFactors = F)
          uctable1 <- rbind(uctable1, temp3)
        }
      }
    }
  }
  
  for(i in 1:nrow(uctable1)){
    cat(i, '\n')
    if(uctable1$Strand[i] == '+'){
      index <- with(uctable1, which(uctable$V1 == Symbol[i] &
                                      uctable$V3 == Chrom[i] &
                                      uctable$V4 == Strand[i] &
                                      uctable$V5 == TSS[i]))
      uctable1$TES[i] <- max(uctable$V6[index])
    }
    else{
      index <- with(uctable1, which(uctable$V1 == Symbol[i] &
                                      uctable$V3 == Chrom[i] &
                                      uctable$V4 == Strand[i] &
                                      uctable$V6 == TSS[i]))
      uctable1$TES[i] <- min(uctable$V5[index])
    }
  }
  
  for(i in 1:length(uniqid)){
    cat(i, '\n')
    index <- grep(paste('^', uniqid[i], '$', sep=''), uctable1$Symbol)
    if(length(index) > 1)
      uctable1$Symbol[index] <- paste(uniqid[i], '#', 1:length(index), sep='')
  }
  
  Promoter4k <- data.frame(chrom = uctable1$Chrom,
                           chromStart = uctable1$TSS - 2000,
                           chromEnd = uctable1$TSS + 2000,
                           name = uctable1$Symbol,
                           Strand = uctable1$Strand
  )
  
}


Wholegenome <- readDNAStringSet('../genomics/genome.fa', use.names = T)
chromlen <- width(Wholegenome)
chromname <- names(Wholegenome)


#### create window function
CreateWindows <- function(data, binsize){
  bedwindow <- foreach (i=1:nrow(data), .options.multicore=list(preschedule=FALSE), 
                        .combine=bind_rows, .inorder=TRUE, .verbose=TRUE, 
                        .errorhandling='stop', .multicombine=TRUE) %dopar% {
                          cat(i, '\n')
                          if(data$Strand[i] == '+')
                            sec <- seq(data$TSS[i], data$TES[i], by = binsize)
                          else
                            sec <- seq(data$TES[i], data$TSS[i], by = binsize)
                          temp <- data.frame(chrom = rep(data$Chrom[i], length(sec)-1),
                                             chromStart = sec[1:(length(sec)-1)],
                                             chromEnd = sec[2:length(sec)],
                                             name = paste(data$Symbol[i], '-', 1:(length(sec)-1),sep=''),
                                             Strand = rep(data$Strand[i], length(sec) -1))
                          temp
                        }
  return(bedwindow)
}

################
### construct promoter regions from ucsc refflat table 
setwd('Dropbox/research/DNAmethylation/')
uctable <- read.table(file= 'hg19refFlat', header = FALSE, sep= '\t', stringsAsFactors = F)
########## Remove chromsomes like chr6_amm_... (Hypolotipe or not able to accurately
########## identity their coordinates)
uctable <- uctable[sapply(strsplit(uctable$V3,'_'), length) == 1,]
CpGisland <- read.table('bam file/CpGisland.bed', sep='\t',stringsAsFactors = F)
CpGisland <- CpGisland[sapply(strsplit(CpGisland$V1,'_'), length) == 1,]
write.table(CpGisland, 'bam\ file/CpGisland.bed', sep='\t', row.names = F, col.names = F, quote=F)

uniqid <- unique(uctable$V1)

uctable1 <- data.frame()
for(i in 1:length(uniqid)){
  cat(i, '\n')
  index <- grep(paste('^', uniqid[i], '$', sep=''), uctable$V1)
  temp <- unique(uctable$V3[index])
  temp1 <- unique(uctable$V4[index])
  for(j in 1:length(temp)){
    for(k in 1:length(temp1)){
      index1 <- which(uctable[index,3] == temp[j] & uctable[index, 4] == temp1[k])
      if(length(index1)){
        if(temp1[k] == '+')
          temp2 <- unique(uctable[index[index1],5])
        else
          temp2 <- unique(uctable[index[index1],6])
        temp3 <- data.frame(Symbol = rep(uniqid[i], length(temp2)),
                            Chrom = rep(temp[j], length(temp2)),
                            Strand = rep(temp1[k], length(temp2)),
                            TSS = temp2, stringsAsFactors = F)
        uctable1 <- rbind(uctable1, temp3)
      }
    }
  }
}


#### get the TES
for(i in 1:nrow(uctable1)){
  cat(i, '\n')
  if(uctable1$Strand[i] == '+'){
    index <- with(uctable1, which(uctable$V1 == Symbol[i] &
                                    uctable$V3 == Chrom[i] &
                                    uctable$V4 == Strand[i] &
                                    uctable$V5 == TSS[i]))
    uctable1$TES[i] <- max(uctable$V6[index])
  }
  else{
    index <- with(uctable1, which(uctable$V1 == Symbol[i] &
                                    uctable$V3 == Chrom[i] &
                                    uctable$V4 == Strand[i] &
                                    uctable$V6 == TSS[i]))
    uctable1$TES[i] <- min(uctable$V5[index])
  }
}

############################################
######### assign unique id to entires with the same gene symbol
for(i in 1:length(uniqid)){
  cat(i, '\n')
  index <- grep(paste('^', uniqid[i], '$', sep=''), uctable1$Symbol)
  if(length(index) > 1)
    uctable1$Symbol[index] <- paste(uniqid[i], '#', 1:length(index), sep='')
}

Promoter4k <- data.frame(chrom = uctable1$Chrom,
                         chromStart = uctable1$TSS - 2000,
                         chromEnd = uctable1$TSS + 2000,
                         name = uctable1$Symbol,
                         Strand = uctable1$Strand
)


