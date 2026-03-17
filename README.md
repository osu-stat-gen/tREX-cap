# tREX-cap
A package that uses a divide and conquer strategy and the tREX algorithm as the building block for constructing whole genome three-dimensional structures.

# Package Creation
```

usethis::create_package(path = "/data/tRexCAP",
                        rstudio = FALSE,
                        check_name = TRUE)
library(usethis)

usethis::use_package('coda')
usethis::use_package('MCMCglmm')
usethis::use_package('Rcpp', min_version = "0.12.0")
usethis::use_package('RcppArmadillo')
usethis::use_package('MASS')

use_package("doParallel", type="Suggest")
use_package("doMC", type="Suggest")


use_rcpp_armadillo('APG_grad.cpp')
use_c("bn.c")
use_c("init.c")
use_c("nz.c")
use_c('pram.c')
use_c('trex')

use_gpl3()

use_author(given = "Jincheol", family = "Park", role = c("aut"))
use_author(given = "Meng", family =  "Wang", email = "cherrywang1006@gmail.com", role = c("aut"))
use_author(given = "Emerson", family = "Webb", email = "webb.944@buckeyemail.osu.edu", role = c("aut"))
use_author(given = "Shili", family = "Lin", email = "shili@stat.osu.edu", role = c("aut", "cre"))


usethis::use_package("RcppGSL", type = "Imports")
usethis::use_package("RcppGSL", type = "LinkingTo")

load('/data/tREX-cap/data/human1422.rda')
usethis::use_data(human1422)

load('/data/tREX-cap/data/sim_bias.rda')
usethis::use_data(sim_bias)

load('/data/tREX-cap/data/sim_bias_human.rda')
usethis::use_data(sim_bias_human)

load('/data/tREX-cap/data/sim_hic.rda')
usethis::use_data(sim_hic)

sim_bias.rda  sim_bias_human.rda  sim_hic.rda


usethis::use_vignette("running-trex-cap")
```
