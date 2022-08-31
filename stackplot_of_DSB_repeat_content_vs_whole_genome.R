#R code to generate the stackplot of RepeatMasker (rmsk) coverage in the DSBs vs genomic coverage

#Load the required libraries:
library(dplyr)
library(ggplot2)

#Set the working directory:
setwd("C:/Users/path/to/input/files")

#Read in the input files:
rmsk_NoGSAT_MM <- read.table("GSAT_MM_removed_rmsk.bed", header=FALSE)
DSB_NoGSAT_MM <- read.table("DSBs_NO_GSAT_MM.bed", header=FALSE)

#Calculate the total bp covered by the DSB peaks as the 100% value for ST DSB coverage calculations:
DSB_NoGSAT_MM$size <- DSB_NoGSAT_MM$V3 -DSB_NoGSAT_MM$V2
total_bp_in_merged_DSB <- DSB_NoGSAT_MM %>% summarise(sum(size))
#total_bp_in_merged_DSB
#sum(size)
#  40478867


#Calculate the coverage of rmsk classes as a percentage of the whole genome (with random, unplaced & chrM regions removed the mm10 genome =2725521370bp) V4 is the rmsk class column
#random unplaced and chrM regions were also removed from the rmsk input file
rmsk_NoGSAT_MM$size <- rmsk_NoGSAT_MM$V3 - rmsk_NoGSAT$V2
rmsk_NoGSAT_MM_cov <- rmsk_NoGSAT_MM %>% group_by(V4) %>% summarise(sum(size)/2725521370*100)


#Read in the bedtools intersect file that contains the number of bp of each repeat masker class that overlapped the spermatid DSB MACS2 peaks
bed_inter_rmsk_classes_v_DSB <- read.table("bed_inter_wao_GSAT_MM_removed_rmsk_norandomUnChrM_vs_SESE_DSB_NO_GSAT_MM", header=FALSE)

#Group the data by the rmsk class (V4) & calculate the total bp per class (V8) and then the % as a percentage of the total bp covered by ST DSBs
bed_inter_rmsk_classes_v_DSB_cov <- bed_inter_rmsk_classes_v_ST_DSB %>% group_by(V4) %>% summarise(sum(V8)/40478867*100)
write.table(bed_inter_rmsk_classes_v_DSB_cov,"pc_cov_of_ST_DSB_macs2_peaks_with_rmsk_classes.txt", row.names = FALSE, quote = FALSE, sep="\t")



#Make a long format file in R/excel with the layout:

#rmsk_class	coverage	sample
#	0.09	DSB_coverage
#LINE	55.25	DSB_coverage
#Low_complexity	0.95	DSB_coverage
#LTR	3.69	DSB_coverage
#Other	0.1	DSB_coverage
#RNA	0	DSB_coverage
#Satellite	2.34	DSB_coverage
#Simple_repeat	20.98	DSB_coverage
#SINE	2.94	DSB_coverage
#Unknown	0.01	DSB_coverage
#DNA	1.079716282	genome_coverage
#LINE	20.07993091	genome_coverage
#Low_complexity	0.748570245	genome_coverage
#LTR	11.70027058	genome_coverage
#Other	0.295896157	genome_coverage
#RNA	0.004349553	genome_coverage
#Satellite	0.177871326	genome_coverage
#Simple_repeat	2.174054353	genome_coverage
#SINE	7.551749998	genome_coverage

#Read in this file to R:
rmsk_DSB_stackplot <- read.table("coverage of DSBs and genome with rmsk classes.txt", header=TRUE)


#Generate the plot, setting the colours to use:
pdf("Stackplot rmsk cov in whole genome and DSB.pdf", width = 10, height = 10, pointsize = 10)
ggplot(rmsk_DSB_stackplot,            # Create ggplot2 plot scaled to 1.00
       aes(x = sample, y = coverage, fill=rmsk_class)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values=c("coral","firebrick","gold","dodgerblue3","deep pink",
                             "dark golden rod", "turquoise1", "green3", "purple",
                             "black")) +
  theme(text = element_text(size =22))
dev.off()

