
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tRexCAP

<!-- badges: start -->
<!-- badges: end -->

tRex-CAP uses a divide and conquer strategy and the tREX algorithm as
the building block for constructing whole genome three-dimensional
structures.

## Installation

You can install tRex-CAP package directly from github by using your
favorite installer (remotes/devtools/pak):

``` r
remotes::install_github("osu-stat-gen/tREX-cap@dev")
```

## Example

To run the cut portion of the algorithm, call `Cut` function directly

``` r
library(tRexCAP)
data(sim_hic)
data(sim_bias)
CutRes = Cut(contact = sim_hic, bias = sim_bias, block_size = 12, noverlap = 1, CPU = 5)
```

To run the paste algorithm, use `Paste` function.

``` r
PasteRes = Paste(contact = sim_hic, cutresult = CutRes, CPU = 5)
```

It is also possible to run both sequentially by calling `CutAndPaste`
function.

``` r
CAPRes = CutAndPaste(contact = sim_hic, block_size = 12, noverlap = 1, CPU = 5)
```
