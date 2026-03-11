\name{trex}
\alias{trex}
\title{
A flagship function.
}
\description{
This is the flagship function of tRex package. Now that Bayesian estimation methods is employed, a user needs to specify jumping rules of proposal distribution in Metropolis-Hastings sampling. For efficiency, Hamiltonian MCMC has been employed in approximating posteriors of coordinates of loci.
}
\usage{
trex(contact, bias, mcmc = 10000, burn = 10000, jump.beta1, jump.x = 0.1, jump.u = 0.1, jump.cov1 = 0.1,jump.cov2 = 0.1,leapfrog.L = 10, leapfrog.e = 0.001, method.type = c("tRex", "tPAM", "PRAM","bn"))
}
\arguments{
  \item{contact}{
This is a nxn matrix typically of HiC library, where n is the number of loci fragments.
}
  \item{bias}{
This is a nx3 matrix in which 1st column should be effective lengths, 2nd be GC contents, and 3rd be mappability. 
}
  \item{mcmc}{
It is a number of MCMC iterations. The default is set to 10000. Thinning procedure is embedded in the code in a way that a number of posterior samples should be 10000. For example, if mcmc=1000000, every sample of 100th iteration is to be collected to make 10000 posterior samples.
}
  \item{burn}{
It is a number of burn-in. The default is set to 10000.
}
  \item{jump.beta1}{
It is a jumping rule of proposal distribution for beta1 which is the coefficient of log distance between loci. 
}
  \item{jump.x}{
It is a jumping rule of proposal distribution for random component X, a random component to induce dependence between loci. It is usually a good practice to set jump.x to make 0.20 acceptance rate.
}
  \item{jump.u}{
It is a jumping rule of proposal distribution for random component U. It is usually a good practice to set jump.u to make 0.20 acceptance rate.
}
  \item{jump.cov1}{
It is a jumping rule of proposal distribution for random component U. It is usually a good practice to set jump.u to make 0.20 acceptance rate.
}
  \item{jump.cov2}{
It is a jumping rule of proposal distribution for random component U. It is usually a good practice to set jump.u to make 0.20 acceptance rate.
}
  \item{leapfrog.L}{
It is a number of steps. This variable is used in leapfrog methods in Hamiltonian MCMC. The default is set 10. A user does not need to modify the defualt value.
}
  \item{leapfrog.e}{
It is a stepsize for leapfrog methods in Hamiltonian MCMC. It is usually a good practice to set this value to make around 0.90 acceptance rate for tRex and PRAM. In the case of tPAM, around 0.75 seems to work well.
}
  \item{method.type}{
A user can select one of three methods: tRex, tPAM, or tRex.
}
}
\value{
  \item{[[1]]}{Posterior samples of coordinates}
  \item{[[2]]}{Posterior samples of parameters associated with the Poisson model}
 }
\examples{
data(sim_hic)
data(sim_bias)
#psample1=trex(contact=sim_hic, bias=sim_bias, mcmc = 10000, burn = 10000, jump.beta1=0.027, jump.x= 0.14, jump.u = 0.14, leapfrog.L = 10, leapfrog.e = 0.001, method.type ="tRex")
psample1=trex(contact=sim_hic, mcmc = 10000, burn = 10000, jump.beta1=0.027, jump.x= 0.14, jump.u = 0.14, leapfrog.L = 10, leapfrog.e = 0.001, method.type ="bn")
}