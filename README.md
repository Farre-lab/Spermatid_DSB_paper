README.txt:

This repository contains the computational scripts (R code) used in the manuscript "A shared "vulnerability code" underpins varying sources of DNA damage throughout paternal germline transmission in mouse." by Frances Burden, et al.

Used to produce Figures 1 & 2 and the supplementary Figures.

code available on Zenodo:
[![DOI](https://zenodo.org/badge/531082884.svg)](https://zenodo.org/badge/latestdoi/531082884)

One script has been included for each separate R plot generated. R scripts requiring more computational power such as RegioneR can be run on a high performance computer (cluster)

This repository contains the following scripts:
1) upsetR.R: Code used to plot Figures 1a/b & supplementary Figure 1.
2) 1kb_corrplot.R:  R Markdown code used to plot the 1kb spermatid DSB corrplot Supplementary Figure 2.
3) Z_score_heatmap.R: Code used to run regioneR on a cluster & the R code used to plot the Z-score heatmap Figure 2a.
4) circos_plot.R: Code used to make the circos plot, Figure 2c.
5) stackplot_of_16_histone_mark_coverage_per_chr.R: Code used to make Supplementary Figure 3, The stackplot of the coverage of the chromosomes with the 16 different histone marks.
6) Stackplot_of_cov_of_states_per_chr.R: Code used to make Supplementary Figure 4, The stackplot of the coverage of the 16 ChromHMM states per chromosome.
7) GAT_log2_fold_corrplot.R: Code used to make the Genomic association tester log2 fold corrplot, Supplementary Figure 5.
8) chr_size_vs_DSB_coverage_per_Mb.R: Code used to produce the correlation plot of chromosome size vs the coverage of spermatid DSBs per Mb, Supplementary Figure 6.
9) stackplot_of_DSB_repeat_content_vs_whole_genome.R: Code used to make Supplementary Figure 7, stackplot of the coverage of RepeatMasker classes in the spermatid DSBs vs genomic coverage.
10) OD_scaling.pl: Perl script used to scale the Oxidative damage data.
