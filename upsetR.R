# R Scripts used to create figures 1a/b UpsetR plots & supplementary figure 1

#Load the required libraries:
library(dplyr)
library(UpSetR)

#Set the working directory:
setwd("C:/Users/location/of/input/output/files")

#The total bp covered by MACS2 peaks per genomic window (either 1kb or 5kb were calculated) & the data for each window was binarized.
#Read in the input files for binarization:
onekb_windows_bp_total_DSB <- read.table("mm10_1kb_total_bp_18_19_20.txt_header_t.txt", header=TRUE)

#Making a merged column so that this can be used in the upSetR file
onekb_windows_bp_total_DSB$region <- paste(onekb_windows_bp_total_DSB$Chr, onekb_windows_bp_total_DSB$window_start, onekb_windows_bp_total_DSB$window_end, sep="-")

#make columns of binary data- so if the region has bp covered by DSB peaks- call the region 1 and zero if not


onekb_windows_bp_total_DSB<- onekb_windows_bp_total_DSB %>% mutate(DSB18_binary=case_when(DSB_18 >=1 ~ 1,
                                                                 DSB_18 ==0 ~0)) 

onekb_windows_bp_total_DSB<- onekb_windows_bp_total_DSB %>% mutate(DSB19_binary=case_when(DSB_19 >=1 ~ 1,
                                                                 DSB_19 ==0 ~0))

onekb_windows_bp_total_DSB <- onekb_windows_bp_total_DSB %>% mutate(DSB20_binary=case_when(DSB_20 >=1 ~ 1,
                                                                 DSB_20 ==0 ~0)) 

#Remove certain columns from the file before saving.
upsetR_1kbDSB_nogsat <- subset(onekb_windows_bp_total_DSB, select=-c(Chr, window_start, window_end))

#Now make the columns in the correct order for the output file:
upset_1kbDSB_nogsatMM_input <- upsetR_1kbDSB_nogsat[, c("region","DSB18_binary","DSB19_binary","DSB20_binary","DSB_18","DSB_19","DSB_20")]
write.csv(upset_1kbDSB_nogsatMM_input, "upset_1kb_DSB_norandom_no_gsatMM_input.csv", row.names=FALSE, quote=FALSE)  

#This example shows the 1kb DSB data. The same scripts were applied to make the 5kb upsetR plot (Suppl_Figure_1, but using 5kb input files)
#For Figure 1b, the DSB and sBLISS files were merged into one to create the upsetR plot using 1kb windows of the mm10 genome.


#Plotting:
pdf("DSB_1kb_window_noGSAT_MM_norandom_UpsetR.pdf")
all_DSB_1kb <- read.csv("upset_1kb_DSB_norandom_no_gsatMM_input.csv", header=TRUE)
upset(all_DSB_1kb, point.size=6)
dev.off()


