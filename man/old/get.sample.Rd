\name{get.sample}
\alias{get.sample}
\title{
This function is designed to get posterior samples of coordinates of specified locus.
}
\description{
An output of tRex() is a list pf two matrices of posterior samples: 1.posterior samples of coordinates of loci 2.posterior samples of parameters. The get.param() is supposed to get posterior samples of coordinates from an output of tRex().
}
\usage{
get.sample(psample, i)
}
\arguments{
  \item{psample}{
This is an output of tRex().
}
  \item{i}{
This argument is supposed to use in specifying locus to get postrior samples.
}
}
\value{
Matrix of columns with 10000 rows. 1st, 2nd, and 3rd colums corresponds respectively to x,y, and z coordinates. Each rows corresponds to a posterior sample  of coordinates.
}
\examples{
library(tRex)
data(sim_hic)
data(sim_bias)

psample2=trex(contact=sim_hic, bias=sim_bias, mcmc = 10000, burn = 1000, jump.beta1=0.018, jump.cov1= 0.3, jump.cov2 = 0.3, leapfrog.L = 10, leapfrog.e = 0.0001, method.type ="tPAM")

S2=get.sample(psample2,2) ### Get posterior samples of second locus.
}

