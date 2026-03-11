\name{summarize.struct}
\alias{summarize.struct}
\title{
This is a function to be used in summarizing posterior samples of loci coordinates.
}
\description{
In summarizing posterior samples of loci coordinates, a user can select one of two methods: posterior mode and posterior mean.
}
\usage{
summarize.struct(psample, method = c("Mode", "Mean"))
}
\arguments{
  \item{psample}{
This is an output of tRex().
}
  \item{method}{
A user is supposed to select between posterior mode and posterior mean.
}
}

\value{
A matrix is yielded of three columns with n rows, where n is a number of loci.
}

\examples{
}
