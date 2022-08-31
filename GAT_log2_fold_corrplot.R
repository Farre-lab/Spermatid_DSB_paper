# R code to generate the corrplot of GAT log2 fold change

#GAT was run and the fold change converted to log2 fold change & added to a data matrix. P-values were added to a separate data matrix

#Loading the required libraries:
library(corrplot)
library(RColorBrewer)

#Set the working directory:
setwd("C:/path/to/input/files")

#Read in the log2 fold change matrix & the p-value matrix:
log2_fold <- read.table("GAT_matrix_log2_fold.txt", header=TRUE)
log2_fold <-data.matrix(log2_fold)

GAT_p_value <- read.table("GAT_matrix_pvalue.txt", header=TRUE)            
GAT_p_value  <-data.matrix(GAT_p_value)

#Set the plot colours:
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(8, "RdBu"))(nb.cols)
rev <-rev(mycolors)

#Plotting:
pdf("corrplot_GAT_fold_log2_matrix.pdf", width=18, height=16)
corrplot(log2_fold, is.corr=FALSE, method ="color", type = c("full"), order = "original", p.mat=GAT_p_value, sig.level = 0.01, insig = "blank",
         tl.col = "black", tl.srt = 90, col=rev, title = "GAT fold log2", mar=c(0,0,1,0), col.lim=c(-15,15))
dev.off()



