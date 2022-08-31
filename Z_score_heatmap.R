# R script used to create figure 2a the heatmap of regioneR Z-scores

# Run regioneR to generate the Z-scores:
#regioneR_BED_vs_BED.sh  (to run the regioneR_BED_vs_BED.R script)

BEDFILE_1="/path/to/bed_file/filename"
BEDlabel_1=$(echo filename)

BEDFILE_2="/path/to/bed_file/filename"
BEDlabel_2=$(echo filename)

Rscript --vanilla /path/to/script/regioneR_BED_vs_BED.R $BEDFILE_1 $BEDFILE_2 ${BEDlabel_1}_vs_${BEDlabel_2} $BEDlabel_1 $BEDlabel_2


#regioneR_BED_vs_BED.R
#!/home/anaconda/envs/r4-base/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

#loading required packages
library(genomation)
library(regioneR)
library(rtracklayer)
library(BiocManager)

#getting the mm10 genome
mm10.genome <- getGenomeAndMask("mm10", mask = NA)$genome

mm10.random_removed <- filterChromosomes(mm10.genome, keep.chr=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                                                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15",
                                                                 "chr16","chr17","chr18", "chr19", "chrX", "chrY"))
# reading in the input BED files
BED_file_1 <- toGRanges(args[1])
BED_file_1 <- keepStandardChromosomes(BED_file_1, pruning.mode="coarse")
#BED_file_1_no_chrM <- dropSeqlevels(BED_file_1, "chrM", pruning.mode="coarse")

BED_file_2 <- toGRanges(args[2])
BED_file_2 <- keepStandardChromosomes(BED_file_2, pruning.mode="coarse")
#BED_file_2 <- dropSeqlevels(BED_file_2, "chrM", pruning.mode="coarse")

#making the permutation test plot
set.seed(123)
pt <- overlapPermTest(A=BED_file_1, B=BED_file_2, ntimes=1000,  mc.set.seed=FALSE, genome=mm10.random_removed, non.overlapping=FALSE, per.chromosome=TRUE)

#making the z-score plot
lz <- localZScore(pt=pt, A=BED_file_1, B=BED_file_2)


x1<-args[4] # label of BED file 1

x2<-args[5] # label of BED file 2

# saving the output to a pdf file
pdf(file=paste(args[3],".pdf"))

par(mfrow=c(2,1))
plot(pt)
title(main="\n\n",
      sub=paste(x1,"vs", x2))

plot(lz)
title(main="\n\n",
      sub=paste(x1,"vs", x2)) 
      
dev.off() 



# Plotting-Used the R-tool pheatmap to generate the heatmap plot:

#Loading the required packages:
library(pheatmap)
library("RColorBrewer")

#Setting the working directory:
setwd("C:/Users/path/to/input/files")

# Manually populated a data matrix with the Z-score values. Values on the diagonal (the same file vs the same file) were given a fixed Z-score of 1000.
# Correlations with a P-value of >= 0.05 (non-significant values) were left empty. These cells will not be filled in the final plot.

#Reading in the input file of regioneR z-scores:
data <- read.table("data_matrix_of_Z-score_values.txt", sep="\t", header=TRUE)
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2 onwards into a matrix
rownames(mat_data) <- rnames                  # assign row names

#Setting the colour palette for the plot:
paletteLength <- 10
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

myBreaks <- c(seq(min(mat_data), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat_data)/paletteLength, max(mat_data), length.out=floor(paletteLength/2)))

#Plot the heatmap:
pdf("pheatmap_RegioneR_z-scores.pdf", height=8, width=8)
pheatmap(mat_data,cex.main=6,color=myColor, breaks=myBreaks, scale ="none", cluster_rows = F, cluster_cols = F, margin = c(5,5), na_col="white", display_numbers=TRUE)            
dev.off()





