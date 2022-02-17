
<!-- README.md is generated from README.Rmd. Please edit that file -->

# parcutils

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.1.0-blue.svg)](https://github.com/cparsania/parcutils)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

The goal of `parcutils` is to provide day to day bioinformatics utility
functions. Most of the functions in the package are useful for analyzing
and visualizing complex RNA-seq studies.

## Installation

``` r

if(require("devtools") && require("BiocManager")){
  options(repos = BiocManager::repositories() )
  devtools::install_github("cparsania/parcutils")
} else{
  install.packages(c("devtools","BiocManager"))
  options(repos = BiocManager::repositories() )
  devtools::install_github("cparsania/parcutils")
}
```

## Example

### Get intergalactic regions and signal intensity from a .bw file.

Given the .gff file and .bw file, identify genome wide intergalactic
regions and their signal intensity.

``` r
library(parcutils)


bw_file <- system.file("extdata" , "example.bw" , package = "parcutils")
gff_file <- system.file("extdata" , "C_glabrata_CBS138_version_s02-m07-r06_features.gff" ,package = "parcutils")

oo <- parcutils::get_intergenic_signals(bw_file = bw_file,gff_file  = gff_file)

oo 
#> GRanges object with 3715 ranges and 1 metadata column:
#>                        seqnames      ranges strand |     score
#>                           <Rle>   <IRanges>  <Rle> | <numeric>
#>      [1] ChrA_C_glabrata_CBS138      1-1607      + |   1.08453
#>      [2] ChrA_C_glabrata_CBS138   2637-2670      + |  14.50112
#>      [3] ChrA_C_glabrata_CBS138  4810-11696      + |   1.98691
#>      [4] ChrA_C_glabrata_CBS138 13043-14976      + |   7.80337
#>      [5] ChrA_C_glabrata_CBS138 15887-17912      + |   4.28257
#>      ...                    ...         ...    ... .       ...
#>   [3711] mito_C_glabrata_CBS138 17988-18003      + |   6.55200
#>   [3712] mito_C_glabrata_CBS138 18076-18084      + |   6.55200
#>   [3713] mito_C_glabrata_CBS138 18157-18179      + |   5.59057
#>   [3714] mito_C_glabrata_CBS138 18263-18280      + |   4.09500
#>   [3715] mito_C_glabrata_CBS138 18368-18403      + |   2.84375
#>   -------
#>   seqinfo: 14 sequences from an unspecified genome; no seqlengths
```

#### Export results

``` r
## export bed file 

oo %>% rtracklayer::export.bed(con = file(paste("intergenic", ".bed",sep = "")))

## export bdg file

oo %>% rtracklayer::export.bedGraph(con = file(paste("intergenic", ".bdg",sep = "")))
```

#### Visualize results

Below is the IGV snapshot showing genes(blue), intergenic regions (red)
and intergenic intensities (green).

![](man/figures/intergenic_snapshot.png)
