\name{get.param}
\alias{get.param}
\title{
This function is designed to get posterior samples of parameters associated with Poisson models.
}
\description{
An output of tRex() consists of two sets of matrices of posterior samples: 1.posterior samples of coordinates of loci fragments 2.posterior samples of parameters. The get.param() is supposed to get posterior samples of parameters from an output of tRex().
}
\usage{
get.param(psample)
}
\arguments{
  \item{psample}{
This is an output of tRex().
}
}
\value{
It is a matrix of 10000 rows and 5 columns of beta1, cov1, cov2, sigma.x, and sigma.u.
}
\examples{
}

