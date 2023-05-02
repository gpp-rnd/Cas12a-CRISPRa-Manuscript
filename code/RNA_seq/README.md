DESeq2
================

The following code is for Figure 3 and Supplementary 3.

``` r
library("data.table")
library("RColorBrewer")
library("dplyr")
library("tibble")
library("stringr")
library("readxl")
library("reshape")
library("DESeq2")
```

    ## Warning: package 'S4Vectors' was built under R version 4.1.3

``` r
library('biomaRt')
library("ggplot2")
par(font.main=6)
```

``` r
coul <- brewer.pal(8, "Set2") 
```

### Control vs VP64_p65_no_guide

![](README_files/figure-gfm/Control%20vs%20VP64_p65_no_guide%20baseline%20vs%20LFC-1.png)<!-- -->
### Parental Control vs VP64_vp64_no_guide

![](README_files/figure-gfm/Control%20vs%20VP64_vp64_no_guide%20baseline%20vs%20LFC-1.png)<!-- -->

### VP64_p65_with_guide vs VP64_p65_no_guide

![](README_files/figure-gfm/Condition_VP64_p65_with_guide_vs_VP64_p65_no_guide%20baseline%20vs%20LFC-1.png)<!-- -->

### VP64_VP64_with_guide vs VP64_VP64_no_guide

![](README_files/figure-gfm/Condition_VP64_VP64_with_guide_vs_VP64_VP64_no_guide%20baseline%20vs%20LFC-1.png)<!-- -->

# Guide Effect

#### 5xtag-dCas12a-VP64 and p65-Nanobody with a cassette containing three CD4-targeting sgRNAs to 5xtag-dCas12a-VP64 and p65-Nanobody without a guide-containing cassette

![](README_files/figure-gfm/VP64_p65_with_guide%20with%20VP64_p65_no_guide%20gene%2050kbs-1.png)<!-- -->

#### 5xtag-dCas12a-VP64 and VP64-Nanobody with a cassette containing three CD4-targeting sgRNAs to 5xtag-dCas12a-VP64 and VP64-Nanobody without a guide-containing cassette

![](README_files/figure-gfm/VP64_VP64_with_guide%20with%20VP64_VP64_no_guide%20gene%2050kbs-1.png)<!-- -->
