# R markdown code used to generate the 1kb corrplot

---
title: "Correlation plot of the total bp covered by MACS2 peaks per 1kb window for DSB18, DSB19, DSB20"
output: 
  word_document: default
  html_document: default

---

```{r setup, include=TRUE}
setwd("C:/Users/path/to/input/files")
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width=12, fig.height=8)
library(corrplot)
library(RColorBrewer)
library(Hmisc)
library(tidyverse)
library(lattice)
```

```
## Reading in the master data file (Chr, start/end already removed)
```{r}
data_1kb <- read.table("1kb_mm10_sum_bpcov_SESE_macs2_peaks_col_per_DSB18_DSB19_DSB20_noChr.txt", header=TRUE)

```

##Use rcorr to make the correlation matrix then extract the correlation coefficients
```{r}
data_corr_p <- rcorr(as.matrix(data_1kb))
```

```

To display the correlation values: `data_corr_p$r data_corr_p$P'

```
##Making the correlation plot

```{r}
corrplot(data_corr_p$r, method = "color", type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 90, col=brewer.pal(n=10, name="RdBu"), title = "corr plot of 1kb windows of bp covered by MACS2 DSB with peaks overlapping GSAT_MM removed", p.mat=data_corr_p$P, sig.level = 0.01, insig = "blank", mar=c(0,0,1,0))

```
