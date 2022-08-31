# This document describes the R code used to make the correlation plot of chromosome size vs the coverage of spermatid DSBs

#Load the required libraries:
library(tidyverse)
library(dplyr)

#Set the working directory:
setwd("C:/User/location/of/input/files")

#Read in the rmsk file of only the GSAT_MM regions & calculate the size of each region:
rmsk_GSAT_MM_no_unplaced_norandom <- read.table("rmsk_GSAT_MM_no_unplaced_norandom.bed", header=FALSE)
rmsk_GSAT_MM_no_unplaced_norandom$size <- rmsk_GSAT_MM_no_unplaced_norandom$V3- rmsk_GSAT_MM_no_unplaced_norandom$V2

#Sum the coverage of GSAT_MM regions per chr: 
GSAT_MM_bp_per_chr <- rmsk_GSAT_MM_no_unplaced_norandom %>% group_by(V1) %>% summarise(total_GSAT_MM_bp=sum(size)) 
write.table(GSAT_MM_bp_per_chr, "GSAT_MM_bp_per_chr.txt", row.names=FALSE, sep="\t", quote = FALSE) 
 
#Use these summed values to manually correct the UCSC sizes in excel by subtracting the total coverage of GSAT_MM regions - save as .txt file
#calculate the coverage per chr of the ST DSB peaks:
#reading in ST DSB input file with no random, unplaced or peaks in GSAT-MM regions:
DSB_norandom_NO_GSAT_MM <- read.table("ST_DSBs_norandom_NO_GSAT_MM.bed", header=FALSE)

#Making a column for peak size:
DSB_norandom_NO_GSAT_MM$size <- DSB_norandom_NO_GSAT_MM$V3 - DSB_norandom_NO_GSAT_MM$V2

#Summing the number of bp covered by MACS2 peaks per chr:
DSB_bp_per_chr_NO_GSAT_MM <- DSB_norandom_NO_GSAT_MM  %>% group_by(V1) %>% summarise(GSAT_MM_removed_total_bp_cov_by_DSB_macs2_peaks=sum(size))

#Add this data to the table of chr size & generate the correlation plot: 
UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem <- read.table("UCSC_chr_size_minus_GSAT_MM_and_bp_cov_by_DSB_minus_GSAT_MM.txt", header=TRUE)

#The input file for the plot had the following format:
#Input file for plot: 
#V1	UCSC_size	GSAT_MM_bp_total	 UCSC_size_minus_GSAT_MM	GSAT_MM_removed_total_bp_cov_by_DSB_macs2_peaks	chr_type 
#chr10	130694993	0	130694993	1907123	autosome 


#Plotting:
pdf("GSAT_MM removed BOTH AXIS Scatter plot chromosome size v the coverage of DSBs.pdf")
ggplot(UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem, aes(x =UCSC_size_minus_GSAT_MM, y =GSAT_MM_removed_total_bp_cov_by_DSB_macs2_peaks, colour=chr_type)) +
  geom_point() +
  geom_text(aes(label=V1),hjust=-0.20, vjust=-0.40, size=2) +
  stat_smooth(method = "lm", col = "blue") +
  theme(axis.text = element_text(size =8),axis.title.y = element_text(size = 10), axis.title.x = element_text(size =10),plot.title = element_text(color = "black", size =8, face = "bold")) +
  ggtitle("The relationship between chromosome size & the coverage of ST DSBs")
dev.off() 


#Calculating a correlation value for the plot:
corr_test_size_no_GSAT_MM_v_macs2bp_removing_GSAT_MM <- cor.test(UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$UCSC_size_minus_GSAT_MM, UCSC_size_GSAT_MM_rem_and_DSB_cov_GSAT_MM_rem$GSAT_MM_removed_total_bp_cov_by_DSB_macs2_peaks, method = "pearson", conf.level = 0.95)
corr_chr_size_no_GSAT_MM_v_macs2_cov_No_GSAT_MM <- data.frame(unlist(corr_test_size_no_GSAT_MM_v_macs2bp_removing_GSAT_MM))
write.table(corr_chr_size_no_GSAT_MM_v_macs2_cov_No_GSAT_MM,"corr_test_size_no_GSAT_MM_v_macs2bp_removing_GSAT_MM.txt", sep="\t", quote=FALSE) 
 

