# Code used to make the stackplot of histone coverage per chr in R

# Calculate the coverage of a histone mark per chr

#load the required libraries:
library(dplyr)
library(tidyverse)

#Set the working directory:
setwd("C:/Users/location/of/input/file(s)")

#Read in the input files:
BRD4 <- read.table("MACS2broad_BCO_0.05_trim20_SRR1596612.1_peaks_sorted_norandom_bed", header=FALSE)

#Calculate the size of the peaks:
BRD4$size <- BRD4$V3 - BRD4$V2

#Calculate the coverage in bp per chr:
BRD4_bp_sum_per_chr <- BRD4 %>% group_by(V1) %>% summarise(sum(size))

#Then in excel calculated the % coverage per chr using the mm10 genome size of 2725521370 (this was with random,unplaced & chrM regions removed)- produced one txt file containing information for all 16 histone marks

##################################################################################

#Generate the stackplot in R

#load the required libraries:
library(ggplot2)

#Set the working directory:
setwd("C:/Users/location/of/histone/coverage/input/file")

#Read in the input file chrX was labelled as chr20 & chrY as 21:
hist_cov <- read.table("coverage of 16 histone marks per chr.txt", header=TRUE)


#Set the colours for the different marks:
cols <- c("five_hMC"="firebrick","BRD4"="coral","H2AZ"="pale golden rod","H3K27ac"="yellow","H3K27me3"="dark golden rod",
          "H3K4me1"="forest green", "H3K4me3"="medium aqua marine", "H3K9ac"="lime green", "H3K9me3"="light steel blue", "H4K12ac"="corn flower blue", "H4K16ac"="deep sky blue",
          "H4K5ac"="blue violet", "H4K8ac"="medium orchid", "H4Kac"="purple", "kac"="plum", "kcr"="deep pink")

M <- factor(hist_cov$mark, levels=c("five_hMC","BRD4","H2AZ","H3K27ac","H3K27me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","H4K12ac","H4K16ac","H4K5ac",
                                    "H4K8ac", "H4Kac","kac","kcr"))

chr <- factor(hist_cov$chr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13", "chr14","chr15","chr16","chr17","chr18", "chr19","chr20","chr21"))

#make the plot:
pdf("Stackplot ST 16 histone mark coverage per Chr.pdf", width = 10, height = 10, pointsize = 10)
ggplot(hist_cov, aes(fill=M, y=coverage, x=chr)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("% coverage of each chr by the 16 histone marks") +
  scale_fill_manual(values=cols) +
  theme_classic() +  
  labs(x='Chr', y="% coverage", fill="histone mark") +
  scale_x_discrete(limits =c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13", "chr14","chr15","chr16","chr17","chr18", "chr19","chr20","chr21"))
dev.off() 


