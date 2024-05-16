
# scDotPlot

<!-- badges: start -->

<!-- badges: end -->

 Dot plots of single-cell RNA-seq data allow for an examination of the relationships between cell groupings (e.g. clusters) and marker gene expression. The scDotPlot package offers a unified approach to perform a hierarchical clustering analysis and add multiple annotations to the columns and/or rows of a scRNA-seq dot plot. It works with SingleCellExperiment and Seurat objects as well as data frames.

## Installation

scDotPlot can be installed from [Bioconductor](https://bioconductor.org/packages/scDotPlot):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
  }

BiocManager::install("scDotPlot")
```

To install the development version directly from GitHub:

```r
if(!requireNamespace("remotes", quietly = TRUE)){
    install.packages("remotes")
  }

remotes::install_github("ben-laufer/scDotPlot")
```

## Examples

See the package [vignette](https://bioconductor.org/packages/release/bioc/vignettes/gg4way/inst/doc/scDotPlot.html).

