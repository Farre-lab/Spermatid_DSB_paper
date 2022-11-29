#This document describes the R code used to make the correlation plot of chromosome size vs the coverage of spermatid DSBs per Mb (Supplementary Figure 6)

#Load the required libraries:
library(tidyverse)
library(dplyr)
library(ggplot2)

#Set the working directory:
setwd("C:/User/location/of/input/files")

#Read in the rmsk file of only the GSAT_MM regions & calculate the size of each region:
rmsk_GSAT_MM_no_unplaced_norandom <- read.table("rmsk_GSAT_MM_no_unplaced_norandom.bed", header=FALSE)
rmsk_GSAT_MM_no_unplaced_norandom$size <- rmsk_GSAT_MM_no_unplaced_norandom$end- rmsk_GSAT_MM_no_unplaced_norandom$start

#Sum the coverage of GSAT_MM regions per chr, the chromosome column =V1: 
GSAT_MM_bp_per_chr <- rmsk_GSAT_MM_no_unplaced_norandom %>% group_by(V1) %>% summarise(total_GSAT_MM_bp=sum(size)) 
write.table(GSAT_MM_bp_per_chr, "GSAT_MM_bp_per_chr.txt", row.names=FALSE, sep="\t", quote = FALSE)
#This file does not contain data for all chromosomes, as not all chromosomes contained GSAT_MM regions.

#use these summed values of bp covered by GSAT_MM per chromosome to manually correct the UCSC chromosome sizes in excel and added a column for the 
#bp covered by DSB per chr (minus GSAT_MM regions). Then imported this file into R.
#Read in this file:
UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem <- read.table("UCSC_chr_size_minus_GSAT_MM_and_bp_cov_by_SESE_DSB_minus_GSAT_MM.txt", header=TRUE)
#File format:
#V1	UCSC_size	GSAT_MM_bp_total	 UCSC_size_minus_GSAT_MM	GSAT_MM_removed_total_bp_cov_by_DSB_macs2_peaks	chr_type 
#chr10	130694993	0	130694993	1907123	autosome 

#Then calculated chromosome length in Mb (excluding the bases covered by GSAT_MM). this was the scaling factor to scale the coverage of the spermatid DSB data to.
UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$OneMb_scaling <- UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$UCSC_size_minus_GSAT_MM/1000000

#Then used the scaling factor tO calculate the coverage of spermatid DSBs scaled per Mb:
UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$bp_covered_by_DSB_per_Mb_No_GSATMM <-UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$GSAT_MM_removed_total_bp_cov_by_SESE_merged_DSB_macs2_peaks/UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$OneMb_scaling

# Plot a scatter plot of chr size against bp covered by DSBs per Mb. This will take into account if longer chr have more DSB because they are longer:

pdf("GSAT_MM removed BOTH AXIS Scatter plot chr size v PER MB coverage of DSBs.pdf")
ggplot(UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem, aes(x =UCSC_size_minus_GSAT_MM, y =bp_covered_by_DSB_per_Mb_No_GSATMM, colour=chr_type)) +
  geom_point() +
  geom_text(aes(label=V1),hjust=-0.20, vjust=-0.40, size=2) +
  stat_smooth(method = "lm", col = "blue") +
  theme(axis.text = element_text(size =8),axis.title.y = element_text(size = 10), axis.title.x = element_text(size =10),plot.title = element_text(color = "black", size =8, face = "bold")) +
  ggtitle("The relationship between chromosome size & the coverage of DSBs scaled per Mb")
dev.off() 

#################################################################################################################
#Calculating the regression stats:

#load the required package:
library("ggpubr")

#calculating the correlation for the plot of chromosome size vs coverage of DSBs scaled per Mb:
corr_test_size_no_GSAT_MM_v_macs2bp_PER_Mb_removing_GSAT_MM <- cor.test(UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$UCSC_size_minus_GSAT_MM, UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$bp_covered_by_DSB_per_Mb_No_GSATMM, method = "pearson", conf.level = 0.95)
corr_chr_size_no_GSAT_MM_v_macs2_cov_PER_Mb_No_GSAT_MM <- data.frame(unlist(corr_test_size_no_GSAT_MM_v_macs2bp_PER_Mb_removing_GSAT_MM))
write.table(corr_chr_size_no_GSAT_MM_v_macs2_cov_PER_Mb_No_GSAT_MM,"corr_chr_size_no_GSAT_MM_v_macs2_cov_PER_Mb_No_GSAT_MM.txt", sep="\t", quote=FALSE) 




