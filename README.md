
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tRexCAP

tRex-CAP (Cut and Paste) uses a divide and conquer strategy and the tREX
algorithm as the building block for constructing whole genome
three-dimensional structures.

# Installation

You can install tRex-CAP package directly from github by using your
favorite installer (remotes/devtools/pak):

``` r
remotes::install_github("osu-stat-gen/tREX-cap",  build_vignettes = TRUE)
```

Note: For OSX, make sure you have gsl and openMP installed before
installing tRex-CAP.

# How To Run tRexCAP

## Data Format

tRexCAP uses the contact matrix. The matrix should be and symmetric,
where the (i, j)-th entry of the matrix is the count of interactions
between bin i and j. The diagnoal elements of the matrix should be 0. An
simulated contact matrix of n=43 bins can be loaded by `data(sim_hic)`.

The optional bias matrix should be of dimension . The 3 columns describe
effective fragment information, GC content, and mappability for each
locus, respectively.

``` r
library(tRexCAP)
data(sim_hic)
data(sim_bias)
```

# Running Cut-and-Paste

Once the data is loaded, Cut-and-Paste algorithm can run separately by
calling `Cut` and `Paste` in two separate steps.

Cut algorithm also requires block_size and no overlap arguments.

Paste requires the contact matrix and the result of Cut algorithm. It
returns the estimated coordinate matrix.

``` r
CutRes = Cut(contact = sim_hic, bias = sim_bias, block_size = 12, noverlap = 1, CPU = 5)
PasteRes = Paste(contact = sim_hic, cutresult = CutRes, CPU = 5)
```

Alternatively, it can be run in a single step by calling `CutAndPaste`

``` r
CAPRes = CutAndPaste(contact = sim_hic, block_size = 12, noverlap = 1, CPU = 5)
```

# Visualization

Once the coordinates are obtained, the output can be visualized with
`draw.strcut`.

``` r
draw.struct=function(coordinates)
```
