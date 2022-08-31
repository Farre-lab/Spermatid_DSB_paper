# This document describes the code used to make the correlation plot of RepeatMasker content per chr vs the coverage of spermatid DSBs in R.

#Load the required libraries:
library(tidyverse)
library(dplyr)

#Set the working directory:
setwd("C:/location/of/input/files")

#Calculate the number of bp covered by RepeatMasker content per chr #(from the rmsk file GSAT_MM regions had been removed)
# The input file is in the format of chr(V1)/start/end
rmsk_no_GSAT_MM <- read.table("GSAT_MM_removed_rmsk.bed", header=FALSE)
rmsk_no_GSAT_MM$size <- rmsk_no_GSAT_MM$V3 -rmsk_no_GSAT_MM$V2

#Calculate the coverage of the different rmsk classes per chr:
rmsk_no_GSAT_MM_bp_per_chr <- rmsk_no_GSAT_MM %>% group_by(V1) %>% summarise(GSAT_MM_removed_rmks_total_bp=sum(size))
write.table(rmsk_no_GSAT_MM_bp_per_chr, "rmsk_no_GSAT_MM_bp_per_chr.txt", row.names=FALSE, sep="\t", quote=FALSE)

#Then manually add the bp covered by DSB peaks (GSAT-MM removed) into this table, a column for the chr number & chromosome type (chr_type) (autosome/chrX,chrY) then re-save
#Read in this new table:  
rmsk_and_DSB_NO_GSAT_MM <- read.table("GSAT_MM_removed_rmsk_bp_per_chr_and_GSAT_MM_rem_DSB_bp_cov_per_chr.txt", header=TRUE)

#Generate the correlation plot:  
pdf("GSAT_MM removed Scatter plot rmsk bp per chr v the coverage of ST DSBs.pdf")
ggplot(rmsk_and_DSB_NO_GSAT_MM, aes(x =GSAT_MM_removed_rmks_total_bp, y =GSAT_MM_removed_total_bp_cov_by_DSB_macs2_peaks, colour=chr_type)) +
  geom_point() +
  geom_text(aes(label=V1),hjust=-0.20, vjust=-0.40, size=2) +
  stat_smooth(method = "lm", col = "blue") +
  theme(axis.text = element_text(size =8),axis.title.y = element_text(size = 10), axis.title.x = element_text(size =10),plot.title = element_text(color = "black", size =8, face = "bold")) +
  ggtitle("The relationship between rmsk coverage per chr & the coverage of ST DSBs")
dev.off()


#Calculating a correlation value for the plot:
corr_test_rmsk_cov_v_macs2bp_removing_GSAT_MM <- cor.test(rmsk_and_DSB_NO_GSAT_MM$GSAT_MM_removed_rmks_total_bp, rmsk_and_DSB_NO_GSAT_MM$GSAT_MM_removed_total_bp_DSB_macs2_peaks, method = "pearson", conf.level = 0.95)
corr_rmsk_cov_no_GSAT_MM_v_macs2_cov_No_GSAT_MM <- data.frame(unlist(corr_test_rmsk_cov_v_macs2bp_removing_GSAT_MM))
write.table(corr_rmsk_cov_no_GSAT_MM_v_macs2_cov_No_GSAT_MM,"corr_rmsk_cov_no_GSAT_MM_v_macs2_cov_No_GSAT_MM.txt", sep="\t", quote=FALSE) 
