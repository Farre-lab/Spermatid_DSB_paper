#R code was used to generate the circos plot Figure 2C.
#circos plot generated using R version 4.0.3 with circlize_0.4.13

#install circlize
install.packages("circlize")
library(circlize)

#Set the working directory:
setwd("C:/Users/path/to/input/files")

#The maximum number of tracks that can easily fit on the plots is 4-5 if the mm10 ideogram is included
#bed files are required to generate the circos plot, with the format of chr/start/end

#Read in the input files:
DSB_merged_NO_GSAT_MM <- read.table("DSB_merged_NO_GSAT_MM_sorted.bed", header=FALSE)
H3C <- read.table("MACS2broad_BCO_0.05_DRR124324.1_trim20_peaks_sorted.bed", header=FALSE)
H3K9me3_ret <- read.table("MACS2broad_BCO_0.05_DRR124334.1_trim20_peaks_sorted.bed", header=FALSE)
H4_ret <- read.table("MACS2broad_BCO_0.05_DRR124336.1_trim20_peaks_sorted.bed", header=FALSE)
BRD4 <- read.table("MACS2broad_BCO_0.05_trim20_SRR1596612.1_peaks_sorted.bed", header=FALSE)


#Making a pdf without zoomed in regions:
pdf("circos merged DSB NO GSAT_MM_red H3C_ret_orange H3k9me3_ret_green H4_ret_blue BRD4_purple.pdf")
circos.clear()
circos.initializeWithIdeogram(species = "mm10")
circos.genomicDensity(DSB_merged_NO_GSAT_MM, col ="red", track.height = 0.1)
circos.genomicDensity(H3C, col ="orange", track.height = 0.09)
circos.genomicDensity(H3K9me3_ret, col ="green", track.height = 0.09)
circos.genomicDensity(H4_ret, col ="blue", track.height = 0.09)
circos.genomicDensity(BRD4, col ="purple", track.height = 0.09)
dev.off()

#Code to zoom in on chr 11 & chrY
extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  rbind(bed, zoom_bed)
}

data = read.cytoband(species = "mm10")
cytoband_df = data$df
chromosome = data$chromosome

xrange = c(data$chr.len, data$chr.len[c("chr11", "chrY")])
normal_chr_index = 1:21
zoomed_chr_index = 22:23


#Making a pdf with zoomed in regions.
pdf("circos merged DSB NO GSAT_MM_red H3C_ret_orange H3k9me3_ret_green H4_ret_blue BRD4_purple_zoom11Y.pdf")
circos.clear()
circos.par(start.degree = 0)
circos.initializeWithIdeogram(extend_chromosomes(cytoband_df, c("chr11", "chrY")), 
                              sector.width = sector.width)

circos.genomicDensity(extend_chromosomes(DSB_merged_NO_GSAT_MM, c("chr11", "chrY")), track.height = 0.1, col="red")
circos.genomicDensity(extend_chromosomes(H3C, c("chr11", "chrY")), track.height = 0.09, col="orange")
circos.genomicDensity(extend_chromosomes(H3K9me3_ret, c("chr11", "chrY")), track.height = 0.09, col="green")
circos.genomicDensity(extend_chromosomes(H4_ret, c("chr11", "chrY")), track.height = 0.09, col="blue")
circos.genomicDensity(extend_chromosomes(BRD4, c("chr11", "chrY")), track.height = 0.09, col="purple")
dev.off()

#The extended regions were manually cut from this plot & added to the non-zoomed plot in Inkscape.

