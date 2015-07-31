library('ggplot2')
setwd('../Dropbox/research/BackgroundCorrection/')
####################### Assume that background noise and true signal have equal p
simulateNB(100000, 10, 100, 0.5)
simulateNB(100000, 10, 100, 0.1)
simulateNB(100000, 10, 100, 0.2)
simulateNB(100000, 10, 100, 0.3)
simulateNB(100000, 10, 100, 0.4)
simulateNB(100000, 10, 100, 0.6)
simulateNB(100000, 10, 100, 0.7)
simulateNB(100000, 10, 100, 0.8)
simulateNB(100000, 10, 100, 0.9)




######################################################################
########################## 4K window promoter region reads 
######################################################################
# HG_LCH_49T <- read.table('Promoter_count/HG_LCH_49T_unique_sorted-Promoter4k.txt', sep = '\t', stringsAsFactors=F)
# HG_LCH_53T <- read.table('Promoter_count/HG_LCH_53T_unique_sorted-Promoter4k.txt', sep = '\t', stringsAsFactors=F)
# HG_LCH_55T <- read.table('Promoter_count/HG_LCH_55T_unique_sorted-Promoter4k.txt', sep = '\t', stringsAsFactors=F)
# HG_LCH_56T <- read.table('Promoter_count/HG_LCH_56T_unique_sorted-Promoter4k.txt', sep = '\t', stringsAsFactors=F)
# T_31A1 <- read.table('Promoter_count/T-31A1_unique_sorted-Promoter4k.txt',sep='\t', stringsAsFactors=F)
# T47K <- read.table('Promoter_count/T47K_unique_sorted-Promoter4k.txt', sep='\t', stringsAsFactors=F)

############ change using fread
HG_LCH_49T <- as.data.frame(fread('bam file/HG_LCH_49T_unique_sorted-Promoter4k.txt'))
HG_LCH_53T <- as.data.frame(fread('bam file/HG_LCH_53T_unique_sorted-Promoter4k.txt'))
HG_LCH_55T <- as.data.frame(fread('bam file/HG_LCH_55T_unique_sorted-Promoter4k.txt'))
HG_LCH_56T <- as.data.frame(fread('bam file/HG_LCH_56T_unique_sorted-Promoter4k.txt'))
T_31A1 <- as.data.frame(fread('bam file/T-31A1_unique_sorted-Promoter4k.txt'))
T47K <- as.data.frame(fread('bam file/T47K_unique_sorted-Promoter4k.txt'))

HG_LCH_49T <- mutate(HG_LCH_49T, V10 = V6*50/V8)
HG_LCH_53T <- mutate(HG_LCH_53T, V10 = V6*50/V8)
HG_LCH_55T <- mutate(HG_LCH_55T, V10 = V6*50/V8)
HG_LCH_56T <- mutate(HG_LCH_56T, V10 = V6*50/V8)
T_31A1 <- mutate(T_31A1, V10 = V6*50/V8)
T47K <- mutate(T47K, V10 = V6*50/V8)

HG_LCH_49T <- mutate(HG_LCH_49T, V11 = V6*10^9/(sum(HG_LCH_49T_100bp[,4])*4000))
HG_LCH_53T <- mutate(HG_LCH_53T, V11 = V6*10^9/(sum(HG_LCH_53T_100bp[,4])*4000))
HG_LCH_55T <- mutate(HG_LCH_55T, V11 = V6*10^9/(sum(HG_LCH_55T_100bp[,4])*4000))
HG_LCH_56T <- mutate(HG_LCH_56T, V11 = V6*10^9/(sum(HG_LCH_56T_100bp[,4])*4000))
T_31A1 <- mutate(T_31A1, V11 = V6*10^9/(sum(T_31A1_100bp[,4])*4000))
T47K <- mutate(T47K, V11 = V6*10^9/(sum(T47K_100bp[,4])*4000))

HG_LCH_49T <- arrange(HG_LCH_49T, V4)
HG_LCH_53T <- arrange(HG_LCH_53T, V4)
HG_LCH_55T <- arrange(HG_LCH_55T, V4)
HG_LCH_56T <- arrange(HG_LCH_56T, V4)
T_31A1 <- arrange(T_31A1, V4)
T47K <- arrange(T47K, V4)

#######################################################################
###################### linux code 
###################### coverage bed to calculate CpG island reads for 
###################### each sample 
for i in $(cat samplelist.txt)
do echo $i.bam 
coverageBed -abam $i.bam -b CpGisland.bed > $i-CpGisland.txt
done
exit 0

for i in $(cat samplelist.txt)
do echo $i.bam 
coverageBed -abam $i.bam -b hg19Promoter4k.bed > $i-Promoter4k.txt
done
exit 0

####################################################################
HG_LCH_49T_CpG <- read.table('bam file/HG_LCH_49T_unique_sorted-CpGisland.txt', sep= '\t', stringsAsFactors=F)
HG_LCH_53T_CpG <- read.table('bam file/HG_LCH_53T_unique_sorted-CpGisland.txt', sep= '\t', stringsAsFactors=F)
HG_LCH_55T_CpG <- read.table('bam file/HG_LCH_55T_unique_sorted-CpGisland.txt', sep= '\t', stringsAsFactors=F)
HG_LCH_56T_CpG <- read.table('bam file/HG_LCH_56T_unique_sorted-CpGisland.txt', sep= '\t', stringsAsFactors=F)
T_31A1_CpG <- read.table('bam file/T-31A1_unique_sorted-CpGisland.txt', sep= '\t', stringsAsFactors=F)
T47K_CpG <- read.table('bam file/T47K_unique_sorted-CpGisland.txt', sep= '\t', stringsAsFactors=F)

HG_LCH_49T_CpG <- mutate(HG_LCH_49T_CpG, V9 = V5*50/V7)
HG_LCH_53T_CpG <- mutate(HG_LCH_53T_CpG, V9 = V5*50/V7)
HG_LCH_55T_CpG <- mutate(HG_LCH_55T_CpG, V9 = V5*50/V7)
HG_LCH_56T_CpG <- mutate(HG_LCH_56T_CpG, V9 = V5*50/V7)
T_31A1_CpG <- mutate(T_31A1_CpG, V9 = V5*50/V7)
T47K_CpG <- mutate(T47K_CpG, V9 = V5*50/V7)

library('rowr')
combined_coverage <- cbind.fill(HG_LCH_49T[10], HG_LCH_53T[10], HG_LCH_55T[10], 
                           HG_LCH_56T[10], T_31A1[10], T47K[10],
                           HG_LCH_49T_CpG[9], HG_LCH_53T_CpG[9], HG_LCH_55T_CpG[9], 
                           HG_LCH_56T_CpG[9], T_31A1_CpG[9], T47K_CpG[9], fill = NA)

colnames(combined_coverage) <- c('HG_LCH_49T', 'HG_LCH_53T', 'HG_LCH_55T',
                                 'HG_LCH_56T', 'T_31A1', 'T47K', 
                                 'HG_LCH_49T_CpG', 'HG_LCH_53T_CpG', 'HG_LCH_55T_CpG',
                                 'HG_LCH_56T_CpG', 'T_31A1_CpG', 'T47K_CpG')

temp <- stack(as.data.frame(combined_coverage))
pdf('promoter_coverage.pdf')
ggplot() + geom_density(aes(x = values, color = ind, linetype = ind), data = temp) +
  scale_colour_manual(values = rep(c("black", "blue", "purple", "gray", "tan3", "red"), each = 2)) +
  scale_linetype_manual(values = rep(c(1,2), times = 6)) +
  guides(colour=guide_legend(''), linetype=guide_legend('')) +
  labs(title='Density plot', x='Coverage', y='Density') + xlim(c(0,5))
dev.off()

temp <- stack(as.data.frame(log2(combined_coverage[,1:6])))
pdf('promoter_coverage_log2.pdf')
ggplot() + geom_density(aes(x = values, color = ind), data = temp) +
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  guides(colour=guide_legend('')) +
  labs(title='Density plot', x='log2 Coverage', y='Density')
dev.off()


#######################################################################################
################# get sequence file of promoter region and calculate cg content
#### linux code 
bedtools getfasta -fi ../../genomics/genome.fa -bed hg19Promoter4k.bed -fo hg19reflat4kpromoter.fa -name
bedtools getfasta -fi ../../genomics/genome.fa -bed CpGisland.bed -fo CpGisland.fa -name
#######################################################################################

#################### using biostring package to calcualte CG concentration of promoter region
Promoter_fa <- readDNAStringSet('bam file/hg19reflat4kpromoter.fa',use.names = T)
CGconc_promoter <- letterFrequency(Promoter_fa, letters = 'CG', as.prob = T)

hg4kpromoter <- fread('bam file/hg19Promoter4k.bed')
hg4kpromoter$V6 <- CGconc_promoter
hg4kpromoter <- arrange(hg4kpromoter, V1, V2, V3)
hg4kpromoter <- as.data.frame(hg4kpromoter)

CGconc_promoter1 <- c()
for(i in 1:nrow(T_31A1)){
  index <- with(hg4kpromoter, which(V1 == T_31A1$V1[i] &
                                      V2 == T_31A1$V2[i] &
                                      V4 == T_31A1$V4[i]))
  CGconc_promoter1 <- c(CGconc_promoter1, hg4kpromoter$V6[index])
  cat(i, '\n')
}

CpGisland_fa <- readDNAStringSet('bam file/CpGisland.fa', use.names = T)
CGconc_CpGisland <- letterFrequency(CpGisland_fa, letters = 'CG', as.prob = T)

CpGisland <- fread('bam file/CpGisland.bed')
CpGisland_unmask <- fread('bam file/CpGisland_unmask.bed')
CGconc_CpGisland1 <- c()
CpGisland$V5 <- CGconc_CpGisland

CGconc_CpGisland <- c()
for(i in 1:nrow(T_31A1_CpG)){
  index <- with(CpGisland, which(V1 == T_31A1_CpG$V1[i] &
                                      V2 == T_31A1_CpG$V2[i] &
                                      V4 == T_31A1_CpG$V4[i]))
  CGconc_CpGisland <- c(CGconc_CpGisland, CpGisland$V5[index])
  cat(i, '\n')
}

Combined_GC <- stack(as.data.frame(combined_coverage))
temp <- cbind.fill(matrix(rep(CGconc_promoter1, times = 6), ncol = 6),
              matrix(rep(CGconc_CpGisland, times = 6), ncol = 6))
colnames(temp) <- c('HG_LCH_49T', 'HG_LCH_53T', 'HG_LCH_55T',
                                 'HG_LCH_56T', 'T_31A1', 'T47K', 
                                 'HG_LCH_49T_CpG', 'HG_LCH_53T_CpG', 'HG_LCH_55T_CpG',
                                 'HG_LCH_56T_CpG', 'T_31A1_CpG', 'T47K_CpG')
temp1 <- stack(temp)
Combined_GC$GCcon <- temp1$values

pdf('GCcon_Cover.pdf')
ggplot(aes(x = GCcon, y = values, color = ind, linetype = ind), data = Combined_GC) + 
  stat_smooth(level = 0.99) + 
  scale_colour_manual(values = rep(c("black", "blue", "purple", "gray", "tan3", "red"), each = 2)) +
  scale_linetype_manual(values = rep(c(1,2), times = 6)) +
  guides(colour=guide_legend(''), linetype=guide_legend('')) +
  labs(title='GC content vs Coverage', x='GC content', y='Coverage') +
  ylim(c(0,2))
dev.off()

#############################################################################
############### uctable and use window size 100bp 
#############################################################################
time coverageBed -abam HG_LCH_49T_unique_sorted.bam -b bed100bp.bed > bed100bp.txt

bedtools getfasta -fi ../../genomics/genome.fa -bed bed100bp.bed -fo bed100bp.fa -name

bed100bp_fa <- readDNAStringSet('bam file/bed100bp.fa',use.names = T)
bed100bp_coverage <- read.delim('bam file/bed100bp.txt', sep='\t',stringsAsFactors = F, header=F)
bed100bp_coverage <- mutate(bed100bp_coverage, V10 = V6*50/V8)
bed100bp_bed <- read.delim('bam file/bed100bp.bed', sep='\t', stringsAsFactors = F, header=F)
CGconc_bed100bp <- letterFrequency(bed100bp_fa, letters = 'CG', as.prob = T)
bed100bp_bed <- mutate(bed100bp_bed, CGconc_bed100bp)
colnames(bed100bp_bed) <- c('Chrom', 'Start', 'End', 'Name', 'Strand', 'CG')
temp <- select(bed100bp_bed, Name, CG)
temp <- arrange(temp, Name)
temp1 <- select(bed100bp_coverage, c(V4,V10))
temp1 <- arrange(temp1, V4)
temp2 <- data.frame(CG = temp$CG, 
                    Coverage = temp1$V10,
                    Name = gsub('*-[0-9]*$', '', temp$Name),
                    stringsAsFactors = F)

# CGconc_bedtry <- c()
# for(i in 1:nrow(bedtry_coverage)){
#   index <- with(bedtry_bed, which(V1 == bedtry_coverage$V1[i] &
#                                    V2 == bedtry_coverage$V2[i] &
#                                    V4 == bedtry_coverage$V4[i] &
#                                     V5 == bedtry_coverage$V5[i]))
#   CGconc_bedtry <- c(CGconc_bedtry, bedtry_bed$V6[index])
#   cat(i, '\n')
# }
# 
# temp <- data.frame(GCconc = CGconc_bedtry,
#                    Coverage = bedtry_coverage$V9,
#                    name = sapply(strsplit(bedtry_coverage$V4, '-'), function(i) i[1])
#                    )

pdf('100bpwindow.pdf')
ggplot(aes(x = CG, y = Coverage, color = Name), data = temp2) + 
  stat_smooth(se = FALSE, method = 'lm') + 
  labs(title='CG content vs Coverage(100bp window)', x='CG content', y='Coverage') +
  guides(color=FALSE)
dev.off()

temp3 <- temp2 %>% 
  group_by(Name) %>% 
  dplyr::summarise(correlation = cor(CG, Coverage, method = 'spearman'))

# library(broom)
# temp4 <- tidy(temp3, model)
# slope_pergene <- temp4[seq(from = 2, to = nrow(temp4), by = 2),]


pdf('Corpergene.pdf')
ggplot() +
  geom_density(aes(x = correlation), data = temp3)
dev.off()

################ try 200bp window 

time coverageBed -abam HG_LCH_49T_unique_sorted.bam -b bed200bp.bed > bed200bp.txt

bedtools getfasta -fi ../../genomics/genome.fa -bed bed200bp.bed -fo bed200bp.fa -name

bed200bp_fa <- readDNAStringSet('bam file/bed200bp.fa',use.names = T)
bed200bp_coverage <- read.delim('bam file/bed200bp.txt', sep='\t',stringsAsFactors = F, header=F)
bed200bp_coverage <- mutate(bed200bp_coverage, V10 = V6*50/V8)
bed200bp_bed <- read.delim('bam file/bed200bp.bed', sep='\t', stringsAsFactors = F, header=F)
CGconc_bed200bp <- letterFrequency(bed200bp_fa, letters = 'CG', as.prob = T)
bed200bp_bed <- dplyr::mutate(bed200bp_bed, CGconc_bed200bp)
colnames(bed200bp_bed) <- c('Chrom', 'Start', 'End', 'Name', 'Strand', 'CG')
temp <- select(bed200bp_bed, Name, CG)
temp <- arrange(temp, Name)
temp1 <- select(bed200bp_coverage, c(V4,V10))
temp1 <- arrange(temp1, V4)
temp2 <- data.frame(CG = temp$CG, 
                    Coverage = temp1$V10,
                    Name = gsub('*-[0-9]*$', '', temp$Name),
                    stringsAsFactors = F)

temp3 <- temp2 %>% 
  group_by(Name) %>% 
  dplyr::summarise(correlation = cor(C.G, Coverage, method = 'spearman'))

pdf('Corpergene_200bp.pdf')
ggplot() +
  geom_density(aes(x = correlation), data = temp3)
dev.off()

#####################################################################################
###### use bedtools to count 100bp window across whole genome
#####################################################################################
time for i in $(ls | grep .bam)
do   
echo $i
coverageBed -abam $i -b Wholegenome_100bp.bed > $i-wholegenome-100bp.txt
done
	   
bedtools getfasta -fi ../../genomics/genome.fa -bed Wholegenome_100bp.bed -fo Wholegenome_100bp.fa -name

Wholegenome_100bp_fa <- readDNAStringSet('bam file/Wholegenome_100bp.fa',use.names = T)

#### data.table's fread is 10 times faster than read.table!!!!!
#### and dont need to specify sep or other extra arguments!!
HG_LCH_49T_100bp <- as.data.frame(fread('bam file/HG_LCH_49T_unique_sorted-wholegenome-100bp.txt'))
HG_LCH_53T_100bp <- as.data.frame(fread('bam file/HG_LCH_53T_unique_sorted-wholegenome-100bp.txt'))
HG_LCH_55T_100bp <- as.data.frame(fread('bam file/HG_LCH_55T_unique_sorted-wholegenome-100bp.txt'))
HG_LCH_56T_100bp <- as.data.frame(fread('bam file/HG_LCH_56T_unique_sorted-wholegenome-100bp.txt'))
T_31A1_100bp <- as.data.frame(fread('bam file/T-31A1_unique_sorted-wholegenome-100bp.txt'))
T47K_100bp <- as.data.frame(fread('bam file/T47K_unique_sorted-wholegenome-100bp.txt'))

HG_LCH_49T_100bp <- mutate(HG_LCH_49T_100bp, V8 = V4*50/V6)
HG_LCH_53T_100bp <- mutate(HG_LCH_53T_100bp, V8 = V4*50/V6)
HG_LCH_55T_100bp <- mutate(HG_LCH_55T_100bp, V8 = V4*50/V6)
HG_LCH_56T_100bp <- mutate(HG_LCH_56T_100bp, V8 = V4*50/V6)
T_31A1_100bp <- mutate(T_31A1_100bp, V8 = V4*50/V6)
T47K_100bp <- mutate(T47K_100bp, V8 = V4*50/V6)

HG_LCH_49T_100bp <- arrange(HG_LCH_49T_100bp, V1,V2,V3)
HG_LCH_53T_100bp <- arrange(HG_LCH_53T_100bp, V1,V2,V3)
HG_LCH_55T_100bp <- arrange(HG_LCH_55T_100bp, V1,V2,V3)
HG_LCH_56T_100bp <- arrange(HG_LCH_56T_100bp, V1,V2,V3)
T_31A1_100bp <- arrange(T_31A1_100bp, V1,V2,V3)
T47K_100bp <- arrange(T47K_100bp, V1,V2,V3)

Coverage_combined_100bp <- bind_cols(HG_LCH_49T_100bp[,1:3, with=F], 
                                     HG_LCH_49T_100bp[,8, with = F],
                                     HG_LCH_53T_100bp[,8, with = F],
                                     HG_LCH_55T_100bp[,8, with = F],
                                     HG_LCH_56T_100bp[,8, with = F],
                                     T_31A1_100bp[,8, with = F],
                                     T47K_100bp[,8, with =F])

colnames(Coverage_combined_100bp) <- c('Chrom', 'Start', 'End', 'Cover1',
                                       'Cover2', 'Cover3', 'Cover4', 
                                       'Cover5', 'Cover6')

Wholegenome_100bp_bed <- fread('Wholegenome_100bp.bed')
GCcont <- letterFrequency(Wholegenome_100bp_fa, letters = 'CG', as.prob = T)
ATCGcont <- letterFrequency(Wholegenome_100bp_fa, letters = 'ATCG', as.prob = T)

Wholegenome_100bp_bed <- dplyr::mutate(Wholegenome_100bp_bed, GCcont, ATCGcont)
Wholegenome_100bp_bed <- as.data.frame(Wholegenome_100bp_bed)
colnames(Wholegenome_100bp_bed) <- c('Chrom', 'Start', 'End', 'CG')

Wholegenome_100bp_bed <- arrange(Wholegenome_100bp_bed, Chrom, Start, End)
Coverage_combined_100bp <- arrange(Coverage_combined_100bp, Chrom, Start, End)

temp <- select(Coverage_combined_100bp, 4:9)
temp <- stack(temp)
temp <- bind_cols(temp, as.data.frame(rep(Wholegenome_100bp_bed[,4], times = 6)))
colnames(temp) <- c('values', 'ind', 'GC')
#temp <- filter(temp, GC != 0)


pdf('GCcon_Cover_100bp.pdf')
ggplot(aes(x = GC, y = values, color = ind), data = temp) + 
  stat_smooth(level = 0.99) + 
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  guides(colour=guide_legend('')) +
  labs(title='GC content vs Coverage', x='GC content', y='Coverage')
dev.off()

#### dplyr too slow and give errors on conditional mutate
# temp <- mutate(temp, group = ifelse(GC<=0.25, 'Group1', 
#                                     ifelse(GC <=0.5, 'Group2', 
#                                            ifelse(GC <=0.75, 'Group3', 'Group4')))) 

####### compare data.table and regular data.frame
##### data.table is much faster than base and dplyr http://stackoverflow.com/questions/24459752/can-dplyr-package-be-used-for-conditional-mutating
setDT(temp)
temp[GC <= 0.25, group := 'group1']
temp[GC > 0.25 & GC <= 0.5, group := 'group2']
temp[GC > 0.5 & GC <= 0.75, group := 'group3']
temp[GC >0.75, group := 'group4']

##### examine the distribution of coverage for quantile groups of GC centent
pdf('GCcon_Cover_100bp_quantile.pdf')
ggplot(aes(x = values, color = ind), data = temp) + 
  geom_density() +
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  guides(colour=guide_legend('')) +
  xlim(c(0,2)) + 
  facet_grid(facets = .~group)
dev.off()





############## get the density of GC contents 100bp windown across human genome
temp <- Wholegenome_100bp_bed %>% 
  select(CG) %>%
  filter(CG!=0)

pdf('GCdist_100bp_wholegenome.pdf')
ggplot() + 
  geom_density(aes(x = CG), data = temp, adjust = 2) +
  labs(x = 'GC content')
dev.off()


#pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")

##################################################################################
# calculate RPKM rather than coverage 

HG_LCH_49T_100bp <- mutate(HG_LCH_49T_100bp, V9 = V4*10^9/(sum(V4)*100))
HG_LCH_53T_100bp <- mutate(HG_LCH_53T_100bp, V9 = V4*10^9/(sum(V4)*100))
HG_LCH_55T_100bp <- mutate(HG_LCH_55T_100bp, V9 = V4*10^9/(sum(V4)*100))
HG_LCH_56T_100bp <- mutate(HG_LCH_56T_100bp, V9 = V4*10^9/(sum(V4)*100))
T_31A1_100bp <- mutate(T_31A1_100bp, V9 = V4*10^9/(sum(V4)*100))
T47K_100bp <- mutate(T47K_100bp, V9 = V4*10^9/(sum(V4)*100))

RPKM_combined_100bp <- bind_cols(HG_LCH_49T_100bp[,1:3, with=F], 
                                     HG_LCH_49T_100bp[,9, with = F],
                                     HG_LCH_53T_100bp[,9, with = F],
                                     HG_LCH_55T_100bp[,9, with = F],
                                     HG_LCH_56T_100bp[,9, with = F],
                                     T_31A1_100bp[,9, with = F],
                                     T47K_100bp[,9, with =F])

colnames(RPKM_combined_100bp) <- c('Chrom', 'Start', 'End', 'RPKM1',
                                   'RPKM2', 'RPKM3', 'RPKM4', 
                                   'RPKM5', 'RPKM6')
RPKM_combined_100bp <- mutate(RPKM_combined_100bp, RPKM_average = (RPKM1+RPKM2+RPKM3+RPKM4+RPKM5+RPKM6)/6)
RPKM_combined_100bp <- arrange(RPKM_combined_100bp, Chrom, Start, End)

temp <- select(RPKM_combined_100bp, 4:9)
temp <- stack(temp)
temp <- bind_cols(temp, as.data.frame(rep(Wholegenome_100bp_bed[,4], times = 6)))
colnames(temp) <- c('values', 'ind', 'GC')
temp <- filter(temp, GC != 0)

pdf('GCcon_RPKM_100bp.pdf')
ggplot(aes(x = GC, y = values, color = ind), data = temp) + 
  stat_smooth(level = 0.99) + 
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  guides(colour=guide_legend('')) +
  labs(title='GC content vs RPKM', x='GC content', y='Coverage')
dev.off()

############## distribution of RPKm for different groups of levels of GC content
######### 
setDT(temp)
temp[GC <= 0.25, group := 'group1']
temp[GC > 0.25 & GC <= 0.5, group := 'group2']
temp[GC > 0.5 & GC <= 0.75, group := 'group3']
temp[GC >0.75, group := 'group4']

##### examine the distribution of RPKM for quantile groups of GC centent
pdf('GCcon_RPKM_100bp_quantile.pdf')
ggplot(aes(x = values, color = ind), data = temp) + 
  geom_density() +
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  guides(colour=guide_legend('')) +
  xlim(c(0,2)) + 
  facet_grid(facets = .~group)
dev.off()

##### examine the overall distribtion of RPKM for each sample 
pdf('RPKM_100bp.pdf')
ggplot(aes(x = log2(values), color = ind), data = temp) + 
  geom_density() +
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  guides(colour=guide_legend(''))
dev.off()

################ cut regions into four quantiles and 
################ take a look at the distribution of CG content 
temp1 <- bind_cols(RPKM_combined_100bp[10], Wholegenome_100bp_bed[4])
temp <- quantile(temp1[,1])
cuts <- cut(temp1[,1], breaks = temp, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
temp1 <- bind_cols(temp1, as.data.frame(cuts))

pdf('QuantileGC.pdf')
ggplot(aes(x = CG, color = cuts), data = temp1) + 
  geom_density()+
  scale_colour_manual(values = c("black", "blue", "purple", "gray"))
dev.off()

pdf('QuantileRPKM.pdf')
ggplot() + 
  geom_histogram(aes(x = RPKM_average, fill = cuts), data = temp1, position = 'identity')+
  scale_fill_manual(values = c("black", "blue", "purple", "gray")) +
  xlim(c(0,2)) +
  facet_grid(facets = .~cuts)
dev.off()

############# take a look at the distribution of raw read
Counts_combined_100bp <- bind_cols(HG_LCH_49T_100bp[,1:3, with=F], 
                                 HG_LCH_49T_100bp[,4, with = F],
                                 HG_LCH_53T_100bp[,4, with = F],
                                 HG_LCH_55T_100bp[,4, with = F],
                                 HG_LCH_56T_100bp[,4, with = F],
                                 T_31A1_100bp[,4, with = F],
                                 T47K_100bp[,4, with =F])

colnames(Counts_combined_100bp) <- c('Chrom', 'Start', 'End', 'Count1',
                                   'Count2', 'Count3', 'Count4', 
                                   'Count5', 'Count6')

temp <- select(Counts_combined_100bp, 4:9)
temp <- stack(temp)

pdf('Countdensity.pdf')
ggplot(aes(x = values, color = ind), data = temp) + 
  geom_density()+
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  xlim(c(0,10))
dev.off()

pdf('Counthist.pdf')
ggplot(aes(x = Count1), data = Counts_combined_100bp) + 
  geom_histogram(aes(y = ..density..), binwidth=1, colour = "darkgreen", fill = "white")+
  xlim(c(0,10)) +
  geom_density()
dev.off()

temp <- Counts_combined_100bp[4]
temp1 <- quantile(temp$Count1)
setDT(temp)
temp[Count1 <= temp1[4], group := 'group1']
temp[Count1 > temp1[4], group := 'group2']

pdf('Counthist_model.pdf')
ggplot() + 
  geom_histogram(aes(x = Count1, y = ..density.., fill=group, alpha = group),  data = temp, binwidth=1, position = 'identity')+
  xlim(c(0,10)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.6, 0.6))
dev.off()

###################################### take a look at the distribution of promoter mapped reads
pdf('Promotercount.pdf')
ggplot() + 
  geom_histogram(aes(x = V6, y = ..density..),  data = HG_LCH_49T, binwidth=1, position = 'identity') +
  xlim(c(0,200))
dev.off()

###################################### create 4k window
time for i in *.bam
do   
echo $i
coverageBed -abam $i -b Wholegenome_4k.bed > $i-wholegenome-4k.txt
done

bedtools getfasta -fi ../../genomics/genome.fa -bed Wholegenome_4k.bed -fo Wholegenome_4k.fa -name

Wholegenome_4k_fa <- readDNAStringSet('bam file/Wholegenome_4k.fa',use.names = T)

#### data.table's fread is 10 times faster than read.table!!!!!
#### and dont need to specify sep or other extra arguments!!
HG_LCH_49T_4k <- fread('HG_LCH_49T_unique_sorted-wholegenome-4k.txt')
HG_LCH_53T_4k <- fread('HG_LCH_53T_unique_sorted-wholegenome-4k.txt')
HG_LCH_55T_4k <- fread('HG_LCH_55T_unique_sorted-wholegenome-4k.txt')
HG_LCH_56T_4k <- fread('HG_LCH_56T_unique_sorted-wholegenome-4k.txt')
T_31A1_4k <- fread('T-31A1_unique_sorted-wholegenome-4k.txt')
T47K_4k <- fread('T47K_unique_sorted-wholegenome-4k.txt')

HG_LCH_49T_4k <- mutate(HG_LCH_49T_4k, V8 = V4*50/V6)
HG_LCH_53T_4k <- mutate(HG_LCH_53T_4k, V8 = V4*50/V6)
HG_LCH_55T_4k <- mutate(HG_LCH_55T_4k, V8 = V4*50/V6)
HG_LCH_56T_4k <- mutate(HG_LCH_56T_4k, V8 = V4*50/V6)
T_31A1_4k <- mutate(T_31A1_4k, V8 = V4*50/V6)
T47K_4k <- mutate(T47K_4k, V8 = V4*50/V6)

HG_LCH_49T_4k <- mutate(HG_LCH_49T_4k, V9 = V4*10^9/(sum(V4)*4000))
HG_LCH_53T_4k <- mutate(HG_LCH_53T_4k, V9 = V4*10^9/(sum(V4)*4000))
HG_LCH_55T_4k <- mutate(HG_LCH_55T_4k, V9 = V4*10^9/(sum(V4)*4000))
HG_LCH_56T_4k <- mutate(HG_LCH_56T_4k, V9 = V4*10^9/(sum(V4)*4000))
T_31A1_4k <- mutate(T_31A1_4k, V9 = V4*10^9/(sum(V4)*4000))
T47K_4k <- mutate(T47K_4k, V9 = V4*10^9/(sum(V4)*4000))

HG_LCH_49T_4k <- arrange(HG_LCH_49T_4k, V1, V2, V3)
HG_LCH_53T_4k <- arrange(HG_LCH_53T_4k, V1, V2, V3)
HG_LCH_55T_4k <- arrange(HG_LCH_55T_4k, V1, V2, V3)
HG_LCH_56T_4k <- arrange(HG_LCH_56T_4k, V1, V2, V3)
T_31A1_4k <- arrange(T_31A1_4k, V1, V2, V3)
T47K_4k <- arrange(T47K_4k, V1, V2, V3)

RPKM_combined_4k <- bind_cols(HG_LCH_49T_4k[,1:3, with=F], 
                                 HG_LCH_49T_4k[,9, with = F],
                                 HG_LCH_53T_4k[,9, with = F],
                                 HG_LCH_55T_4k[,9, with = F],
                                 HG_LCH_56T_4k[,9, with = F],
                                 T_31A1_4k[,9, with = F],
                                 T47K_4k[,9, with =F])

colnames(RPKM_combined_4k) <- c('Chrom', 'Start', 'End', 'RPKM1',
                                   'RPKM2', 'RPKM3', 'RPKM4', 
                                   'RPKM5', 'RPKM6')
RPKM_combined_4k <- mutate(RPKM_combined_4k, RPKM_average = (RPKM1+RPKM2+RPKM3+RPKM4+RPKM5+RPKM6)/6)
RPKM_combined_4k <- arrange(RPKM_combined_4k, Chrom, Start, End)

Wholegenome_4k_bed <- fread('bam file/Wholegenome_4k.bed')
GCcont <- letterFrequency(Wholegenome_4k_fa, letters = 'CG', as.prob = T)
ATCGcont <- letterFrequency(Wholegenome_4k_fa, letters = 'ATCG', as.prob = T)
Wholegenome_4k_bed <- dplyr::mutate(Wholegenome_4k_bed, GCcont, ATCGcont)
Wholegenome_4k_bed <- as.data.frame(Wholegenome_4k_bed)
colnames(Wholegenome_4k_bed) <- c('Chrom', 'Start', 'End', 'CG','ATCG')
Wholegenome_4k_bed <- arrange(Wholegenome_4k_bed, Chrom, Start, End)

temp <- select(RPKM_combined_4k, 4:9)
temp <- stack(temp)
temp <- bind_cols(temp, as.data.frame(rep(Wholegenome_4k_bed[,4], times = 6)))
colnames(temp) <- c('values', 'ind', 'GC')

pdf('GCcon_RPKM_4k.pdf')
ggplot(aes(x = GC, y = values, color = ind), data = temp) + 
  stat_smooth(level = 0.99) + 
  scale_colour_manual(values = c("black", "blue", "purple", "gray", "tan3", "red")) +
  guides(colour=guide_legend('')) +
  labs(title='GC content vs RPKM', x='GC content', y='Coverage')
dev.off()

temp1 <- bind_cols(RPKM_combined_4k[10], Wholegenome_4k_bed[4])
temp <- quantile(temp1[,1])
cuts <- cut(temp1[,1], breaks = temp, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
temp1 <- bind_cols(temp1, as.data.frame(cuts))

######### RPKM average for all samples promoter region 
RPKM_combined_promoter <- bind_cols(HG_LCH_49T[,1:3, with=F], 
                              HG_LCH_49T[,11, with = F],
                              HG_LCH_53T[,11, with = F],
                              HG_LCH_55T[,11, with = F],
                              HG_LCH_56T[,11, with = F],
                              T_31A1[,11, with = F],
                              T47K[,11, with =F])
colnames(RPKM_combined_promoter) <- c('Chrom', 'Start', 'End', 'RPKM1',
                                'RPKM2', 'RPKM3', 'RPKM4', 
                                'RPKM5', 'RPKM6')
RPKM_combined_promoter <- mutate(RPKM_combined_promoter, RPKM_average = (RPKM1+RPKM2+RPKM3+RPKM4+RPKM5+RPKM6)/6)
RPKM_combined_promoter <- arrange(RPKM_combined_promoter, Chrom, Start, End)
hg4kpromoter <- arrange(hg4kpromoter, V1, V2, V3)


temp1 <- bind_cols(RPKM_combined_4k[10], Wholegenome_4k_bed[4])
temp <- quantile(temp1[,1])
cuts <- cut(temp1[,1], breaks = temp, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
temp1 <- bind_cols(temp1, as.data.frame(cuts), as.data.frame(rep('Allregions', nrow(temp1))))
colnames(temp1) <- c('RPKM', 'CG', 'cuts', 'type')

temp3 <- bind_cols(RPKM_combined_promoter[10], as.data.frame(hg4kpromoter[,V6]))
temp4 <- quantile(temp3[,1])
cuts <- cut(temp3[,1], breaks = temp4, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
temp3 <- bind_cols(temp3, as.data.frame(cuts), as.data.frame(rep('Promoter', nrow(temp3))))
colnames(temp3) <- c('RPKM', 'CG', 'cuts', 'type')

temp5 <- bind_rows(temp1, temp3)

# quantile of GC
pdf('QuantileGC_4k_combined.pdf')
ggplot() + 
  geom_density(aes(x = CG, color = cuts, linetype = type), adjust = 2, data = temp5)+
  scale_colour_manual(values = c("black", "blue", "purple", "gray")) + 
  scale_linetype_manual(values = c(1,2))
dev.off()

###################### distribution of count in 4k window
pdf('Counthist_4k.pdf')
ggplot(aes(x = V4), data = HG_LCH_49T_4k) + 
  geom_histogram(aes(y = ..density..), binwidth = 1, colour = "darkgreen", fill = "white") +
  xlim(c(0,200))
dev.off()

############### combine the 4k with promoter 4k 
temp <- bind_cols(select(HG_LCH_49T_4k, V4), 
                  as.data.frame(rep('Overall', nrow(HG_LCH_49T_4k))))
colnames(temp) <- c('count', 'group')

temp1 <- bind_cols(select(HG_LCH_49T, V6), 
                   as.data.frame(rep('Promoter', nrow(HG_LCH_49T))))
colnames(temp1) <- c('count', 'group')

temp2 <- bind_rows(temp, temp1)

pdf('Counthist_4k_overallvspromoter.pdf')
ggplot() + 
  geom_histogram(aes(x = count, y = ..density.., fill=group, alpha = group),  data = temp2, binwidth=1, position = 'identity')+
  xlim(c(0,150)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.6, 0.6))
dev.off()

############ take only the 4k windows that are GC low as background noise
############ plot the count reads distribution for GC low (<40%) and Promoter regions
temp <- bind_cols(HG_LCH_49T_4k[4], Wholegenome_4k_bed[4])
temp <- filter(temp, CG < 0.3)
temp <- mutate(temp, group = 'Background')
temp <- select(temp, V4, group)
colnames(temp) <- c('count', 'group')

temp1 <- bind_cols(select(HG_LCH_49T, V6), 
                   as.data.frame(rep('Promoter', nrow(HG_LCH_49T))))
colnames(temp1) <- c('count', 'group')

temp2 <- bind_rows(temp, temp1)

pdf('Counthist_4k_backgroundvspromoter.pdf')
ggplot() + 
  geom_histogram(aes(x = count, y = ..density.., fill=group, alpha = group),  data = temp2, binwidth=1, position = 'identity')+
  xlim(c(0,150)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.6, 0.6))
dev.off()


################################################################################################
#################################################################################################
############# exclude promoter regions and CpG islands regions from 100bp window
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

# temp <- union(Exclude_index, Exclude_index1)
# temp1 <- bind_cols(as.data.frame(RPKM_combined_100bp[-temp,10]), as.data.frame(Wholegenome_100bp_bed[-temp, 4]))
# temp2 <- quantile(temp1[,1])
# cuts <- cut(temp1[,1], breaks = temp2, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
# temp1 <- bind_cols(temp1, as.data.frame(cuts))
# colnames(temp1) <- c('RPKM', 'CG', 'cuts')
# 
# # quantile of GC
# pdf('QuantileGC_100bp_exclude.pdf')
# ggplot() + 
#   geom_density(aes(x = CG, color = cuts), adjust = 2, data = temp1)+
#   scale_colour_manual(values = c("black", "blue", "purple", "gray"))
# dev.off()
# 
# ########## take a look at the promote and CpG islands
# temp <- union(Exclude_index, Exclude_index1)
# temp1 <- bind_cols(as.data.frame(RPKM_combined_100bp[temp,10]), as.data.frame(Wholegenome_100bp_bed[temp, 4]))
# temp2 <- quantile(temp1[,1])
# cuts <- cut(temp1[,1], breaks = temp2, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
# temp1 <- bind_cols(temp1, as.data.frame(cuts))
# colnames(temp1) <- c('RPKM', 'CG', 'cuts')
# 
# # quantile of GC
# pdf('QuantileGC_100bp_promoter_cpgisland.pdf')
# ggplot() + 
#   geom_density(aes(x = CG, color = cuts), adjust = 2, data = temp1)+
#   scale_colour_manual(values = c("black", "blue", "purple", "gray"))
# dev.off()
temp1 <- bind_cols(as.data.frame(RPKM_combined_100bp[,10]), as.data.frame(Wholegenome_100bp_bed[, 4]))
temp2 <- quantile(temp1[,1])
cuts <- cut(temp1[,1], breaks = temp2, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
temp1 <- bind_cols(temp1, as.data.frame(cuts), as.data.frame(rep('Allregions', nrow(temp1))))
colnames(temp1) <- c('RPKM', 'CG', 'cuts', 'type')

temp <- union(Exclude_index, Exclude_index2)
temp3 <- bind_cols(as.data.frame(RPKM_combined_100bp[temp,10]), as.data.frame(Wholegenome_100bp_bed[temp, 4]))
temp4 <- quantile(temp3[,1])
cuts <- cut(temp3[,1], breaks = temp2, labels = c('0-25%', '25-50%', '50-75%', '75-100%'))
temp3 <- bind_cols(temp3, as.data.frame(cuts), as.data.frame(rep('Promoter', nrow(temp3))))
colnames(temp3) <- c('RPKM', 'CG', 'cuts', 'type')

temp5 <- bind_rows(temp1, temp3)

# quantile of GC
pdf('QuantileGC_100bp_combined1.pdf')
ggplot() + 
  geom_density(aes(x = CG, color = cuts, linetype = type), adjust = 2, data = temp5)+
  scale_colour_manual(values = c("black", "blue", "purple", "gray")) + 
  scale_linetype_manual(values = c(1,2))
dev.off()

###################################### use 100bp all-regions as observed signal, promoter regions
###################################### as background noise
pdf('RPKM_allvspromoter_100bp.pdf')
ggplot() + 
  geom_histogram(aes(x = RPKM, y = ..density.., fill=type, alpha = type),  data = temp5, position = 'identity')+
  xlim(c(0,1.5)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.6, 0.6))
dev.off()

############## examine the distribution of count 100bp window
temp <- bind_cols(as.data.frame(HG_LCH_49T_100bp[,V4]), 
                  as.data.frame(rep('Overall', nrow(HG_LCH_49T_100bp))))
colnames(temp) <- c('count', 'group')

temp1 <- bind_cols(as.data.frame(HG_LCH_49T_100bp[unique(Exclude_index), V4]), 
                   as.data.frame(rep('Promoter', length(unique(Exclude_index)))))
colnames(temp1) <- c('count', 'group')

temp2 <- bind_rows(temp, temp1)

pdf('Counthist_100bp_overallvspromoter.pdf')
ggplot() + 
  geom_histogram(aes(x = count, y = ..density.., fill=group, alpha = group),  data = temp2, binwidth=1, position = 'identity')+
  xlim(c(0,15)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.8, 0.5))
dev.off()

########################### Exclude promoter regions from 4k window
#############################################################################################
#############################################################################################
registerDoMC(80)
t1 <- Sys.time()
Exclude_index_4k <- ExcludeIntersection(Wholegenome_4k_bed, hg4kpromoter)
t2 <- Sys.time()
t2 - t1

temp <- unique(Exclude_index_4k)
temp3 <- bind_cols(as.data.frame(RPKM_combined_4k[-temp,10]), as.data.frame(Wholegenome_4k_bed[-temp, 4]))
temp4 <- quantile(temp3[,1])
cuts <- cut(temp3[,1], breaks = unique(temp4), labels = c('0-25%','25-50%', '50-75%', '75-100%'))
temp3 <- bind_cols(temp3, as.data.frame(cuts), as.data.frame(rep('Promoter', nrow(temp3))))
colnames(temp3) <- c('RPKM', 'CG', 'cuts', 'type')

pdf('test.pdf')
ggplot() + 
  geom_density(aes(x = CG, color = cuts), adjust = 2, data = temp3)+
  scale_colour_manual(values = c("black", "blue", "purple", "gray"))
dev.off()








#############################################################################################
######################################### try 200 bp window 
##############################################################################################
write.table(Wholegenome_200bp, 'bam file/Wholegenome_200bp.bed', sep='\t', quote=F, row.names = F, col.names=F)

for i in $(cat samplelist.txt)
do echo $i.bam 
coverageBed -abam $i.bam -b Wholegenome_200bp.bed > $i-Wholegenome_200bp.txt
done
exit 0





#############################################################################################
##################################### identify regions for each gene to identify the 
##################################### background noise
################################## exclude promoter regions, predicted CpG islands regions
############################### and regions with 'N', no ATCG, 
hg4kpromoter <- mutate(hg4kpromoter, V7 = V2+2000)
# hg4kpromoter_unique <- hg4kpromoter %>% 
#   group_by(V1) %>%
#   distinct(V7)
Exclude_index3 <- which(Wholegenome_100bp_bed$ATCG != 1)
#Exclude_index_100bp <- dplyr::union(Exclude_index, Exclude_index3)
Exclude_index_100bp <- Reduce(dplyr::union, list(Exclude_index, Exclude_index2, Exclude_index3))

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

temp <- HG_LCH_49T[,c(4,11), with = F]
temp <- arrange(temp, V4)
temp <- mutate(temp, group = 'Promoter')
colnames(temp) <- c('name', 'RPKM_average', 'group')

temp1 <- bind_rows(Background_region_aggre, temp)

pdf('RPKMhist_promotervsbackground.pdf')
ggplot() + 
  geom_histogram(aes(x = RPKM_average, y = ..density.., fill=group, alpha = group),  data = temp1, position = 'identity')+  scale_fill_manual(values = c('blue', 'red')) +
  xlim(c(0,1)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.8, 0.5))
dev.off()

pdf('RPKMlog2hist_promotervsbackground.pdf')
ggplot() + 
  geom_histogram(aes(x = log2(RPKM_average), y = ..density.., fill=group, alpha = group),  data = temp1, position = 'identity')+  scale_fill_manual(values = c('blue', 'red')) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.8, 0.5))
dev.off()

write.table(select(Background_region, 1:3,7), 'bam file/Background_regions.bed', sep='\t', row.names=F, col.names = F, quote=F)

############## count reads in background region
time for i in $(cat samplelist.txt)
do echo $i.bam 
coverageBed -abam $i.bam -b Background_regions.bed > $i-Background_regions.txt
done
exit 0

########## carry out DE analysis for prostate cancer patients 
############ summerize regions to gene level
colnames(Background_count_combined) <- c(paste('S', 1:6, sep=''))
Background_count_combined <- mutate(Background_count_combined, name = HG_LCH_49T_background$V4)
Background_count_combined <- Background_count_combined %>% group_by(name) %>% summarise_each(funs(sum)) %>% arrange(name) %>% select(-1)
rownames(Background_count_combined) <- rownames(Promoter_Count_Combined)

conditions <- as.factor(c(rep('C1', 2), rep('C2', 3)))
XB <- XBSeqDataSet(Promoter_Count_Combined[,-3], Background_count_combined[,-3], conditions)
XB <- estimateRealCount(XB)
XB <- estimateSizeFactors(XB)
Normalized_mat <- counts(XB, normalized = T)
write.table(Normalized_mat, 'Normalized_mat_prostate_all.txt', sep='\t', row.names = T, col.names = T, quote=F)
XB <- estimateSCV(XB, method = "pooled", sharingMode = "maximum",  fitType = "local")
Teststas <- XBSeqTest(XB, levels(conditions)[1L], levels(conditions)[2L])
write.table(Teststas, 'Teststas_prostate_all.txt', sep='\t', row.names = T, col.names = T, quote=F)
DE_id <- Teststas$id[which(abs(Teststas$log2FoldChange) > 1 & Teststas$pval < 0.01 & Teststas$baseMean > 15)]
DE_id <- sapply(DE_id, function(i) strsplit(i, '#')[[1]][1])
write.table(DE_id, 'DE_id_prostate.txt', sep='\t', row.names = F, col.names = F, quote=F)

######################## 3d mds plot based on DNA methylation 
dataMat <- fread('research/DNAmethylation/Normalized_mat_prostate.txt',data.table = F)
rownames(dataMat) <- dataMat$V1
dataMat <- dataMat[,-1]
colnames(dataMat) <- c(paste('T', 1:2, sep= ''), paste('C', 1:3, sep= ''))
pcaplot(as.matrix(log2(dataMat+1)), method = 'mds', text =T, cv.Th = 1.8, psi=8, color = as.factor(c(rep('T', 2), rep('N', 3))))

rgl.postscript('research/DNAmethylation/pca3d_prostate.pdf','pdf')

### heatmap
mean_gene <- apply(log2(dataMat+1), 1, mean)
var_gene <- apply(log2(dataMat+1), 1, sd)
dataMat_sel <- dataMat[mean_gene > 4 & var_gene > 1,]
pdf('research/DNAmethylation/heatmap_prostate.pdf')
heatmap.my(as.matrix(log2(dataMat_sel+1)),colsidebar = as.factor(c(rep('T', 2), rep('N', 3))), 
           breakratio = c(4,1,4))
dev.off()











############### identify 80 windows for each promoter based on GC content 
registerDoMC(80)
t1 <- Sys.time()
Background_region_80window <- foreach (i=1:nrow(hg4kpromoter), .options.multicore=list(preschedule=FALSE), 
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
                                bgregion <- bgregion %>% mutate(name = hg4kpromoter[i,4]) %>% select(Chrom, Start, End, name)
                                bgregion
                              }
t2 <- Sys.time()
t2 - t1

write.table(Background_region, 'bam file/Background_regions_80.bed', sep='\t', row.names=F, col.names = F, quote=F)

################## simulation studies to compare between MEDIPS vs MBDcap,
ROI <- HG_LCH_49T[,1:4]
ME_group1_1 <- MEDIPS.createROIset(file='MEDIPS/C1/temp.txt', ROI=ROI, BSgenome= 'BSgenome.Hsapiens.UCSC.hg19')
ME_group1_2 <- MEDIPS.createROIset(file='MEDIPS/C1/temp1.txt', ROI=ROI, BSgenome= 'BSgenome.Hsapiens.UCSC.hg19')
ME_group2_1 <- MEDIPS.createROIset(file='MEDIPS/C2/temp2.txt', ROI=ROI, BSgenome= 'BSgenome.Hsapiens.UCSC.hg19')
ME_group2_2 <- MEDIPS.createROIset(file='MEDIPS/C2/temp3.txt', ROI=ROI, BSgenome= 'BSgenome.Hsapiens.UCSC.hg19')

DM_res <- MEDIPS.meth(MSet1 = list(ME_group1_1, ME_group1_2), MSet2 = list(ME_group2_1, ME_group2_2))
DM_res_sel <-  MEDIPS.selectSig(results = DM_res, p.value = 0.1, adj = T, ratio = 2, bg.counts = NULL, CNV = F)


























############ Exclude GC content above 0.4 from wholegenome 100 bp data
index <- which(Wholegenome_100bp_bed$CG > 0.4)
temp <- HG_LCH_49T_100bp[-union(Exclude_index_100bp, index),V4]

temp <- bind_cols(as.data.frame(temp), 
                  as.data.frame(rep('Overall_lowGC', length(temp))))
colnames(temp) <- c('count', 'group')

temp1 <- bind_cols(as.data.frame(HG_LCH_49T_100bp[unique(Exclude_index), V4]), 
                   as.data.frame(rep('Promoter', length(unique(Exclude_index)))))
colnames(temp1) <- c('count', 'group')

temp2 <- bind_rows(temp, temp1)

pdf('Counthist_100bp_overall_lowGCvspromoter.pdf')
ggplot() + 
  geom_histogram(aes(x = count, y = ..density.., fill=group, alpha = group),  data = temp2, binwidth=1, position = 'identity')+
  xlim(c(0,15)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.8, 0.5))
dev.off()

############ Exclude GC content above 0.4 from wholegenome 4kb data
index <- which(Wholegenome_4k_bed$CG > 0.4 | Wholegenome_4k_bed$ATCG != 1)
index <- which(Wholegenome_4k_bed$ATCG != 1)
temp <- HG_LCH_49T_4k[-union(Exclude_index_4k, index), 4]

temp <- bind_cols(as.data.frame(temp), 
                  as.data.frame(rep('Overall_lowGC', length(temp))))
colnames(temp) <- c('count', 'group')

temp1 <- bind_cols(select(HG_LCH_49T, V6), 
                   as.data.frame(rep('Promoter', nrow(HG_LCH_49T))))
colnames(temp1) <- c('count', 'group')

temp2 <- bind_rows(temp, temp1)

pdf('Counthist_4k_overallvspromoter.pdf')
ggplot() + 
  geom_histogram(aes(x = count, y = ..density.., fill=group, alpha = group),  data = temp2, binwidth=1, position = 'identity')+
  xlim(c(0,200)) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_alpha_manual(values = c(0.8, 0.5))
dev.off()












############# Also cut regions into four quantiles



######## examine the distribution of correlation


######## promoter RPKM vs GC content


######## promoter RPKM density (or log2)


########### table low coverage and high coverage 

############ do a test 
############ try to choose 40 100bp window that has close Gc content as one selected promoter
########### examine the distribution








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



simulateNB <- function(n, size1, size2, prob, seed = 1990){
  set.seed(1990)
  background <- data.frame(value = rnbinom(n, size = size1, prob = prob))
  signal <- data.frame(value = rnbinom(n, size = size2, prob = prob))
  observed <- data.frame(value = rnbinom(n, size = size1 + size2, prob=prob))
  sumed <- background + signal
  ggplot() + geom_density(aes(x=value),colour='green', data= background) +
    geom_density(aes(x=value), colour='red',data= signal) +
    geom_density(aes(x=value), colour='blue',data= observed) +
    geom_density(aes(x=value), colour='orange',data= sumed)
}







########################################################test script
temp1 <- bind_cols(as.data.frame(RPKM_combined_100bp[,10]), as.data.frame(Wholegenome_100bp_bed[, 4]))
temp2 <- quantile(temp1[,1])
cuts <- cut(temp1[,1], breaks = unique(temp2), labels = c('0-25%','25-50%', '50-75%', '75-100%'))
temp1 <- bind_cols(temp1, as.data.frame(cuts), as.data.frame(rep('Allregions', nrow(temp1))))
colnames(temp1) <- c('RPKM', 'CG', 'cuts', 'type')

temp <- union(Exclude_index, Exclude_index1)
temp3 <- bind_cols(as.data.frame(RPKM_combined_100bp[-temp,10]), as.data.frame(Wholegenome_100bp_bed[-temp, 4]))
temp4 <- quantile(temp3[,1])
cuts <- cut(temp3[,1], breaks = unique(temp4), labels = c('0-25%','25-50%', '50-75%', '75-100%'))
temp3 <- bind_cols(temp3, as.data.frame(cuts), as.data.frame(rep('Promoter', nrow(temp3))))
colnames(temp3) <- c('RPKM', 'CG', 'cuts', 'type')

temp5 <- bind_rows(temp1, temp3)

# quantile of GC
pdf('test1.pdf')
ggplot() + 
  geom_density(aes(x = CG, color = cuts), adjust = 2, data = temp3)+
  scale_colour_manual(values = c("black", "blue", "purple", "gray"))
dev.off()




