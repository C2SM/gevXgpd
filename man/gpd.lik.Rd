\name{gpd.lik}
\alias{gpd.lik}
\alias{expon.lik}
\alias{gpd.lik.retv}
\alias{expon.lik.retv}
\title{ Likelihood function of GPD }
\description{
  The negative log likelihood function for the Generalized Pareto 
  Distribution and the exponential distribution.
  These functions are used to calculate maximum likelihood 
  estimates and likelihood profiles for the 
  parameters (\code{gpd.lik}, \code{expon.lik}) and for the return values 
  (\code{gpd.lik.retv}, \code{expon.lik.retv}). This implementation 
  is for the three parameter 
  problem where the frequency of exceedances \code{lamda} 
  is a parameter too. (This is 
  useful for confidence interval calculation.) 

  \code{gpd.lik.retv} and \code{expon.lik.retv} are reparameterisations of
  the likelihood function using \code{lamda}, the return value 
  (at a predefined return period \code{ret}), and shape \code{k} (in case of
  gpd to specify the distribution. 

  This function is usually not called by the user directly.
}
\usage{
gpd.lik(pars, data, chi=0, tim=length(data))
expon.lik(pars, data, chi=0, tim=length(data))
gpd.lik.retv(pars, data, chi=0, tim=length(data), ret)
expon.lik.retv(pars, data, chi=0, tim=length(data), ret)
}
\arguments{
  \item{pars}{ A vector with the parameter values. For \code{gpd.lik} this
      is the exceedance frequency \code{lamda}, 
      the scale parameter \code{alpha}
      and the shape parameter \code{k} (in this order). 
      For \code{gpd.lik.retv}
      \code{pars} contains exceedance frequency, return value, 
      and shape (in this order). }
      For \code{expon.lik} and \code{expon.lik.retv} only 2d without shape. 
  \item{data}{ A vector with the data sample. }
  \item{chi}{ The exceedance threshold. If default \code{data} are 
      exceedances directly.  }
  \item{tim}{ The total time over which the data is spanning. This is 
      needed to calculate the likelihood of the number of exceedances. }
  \item{ret}{ A scalar number giving the return period (in units of the
      average time between exceedances !!! ). }
}
\details{

}
\value{
  A real value giving the negative log likelihood.
}
\references{ see e.g. Coles 2001 }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{fitGPD}}, \code{\link{conf.bounds.pot}}, 
          \code{\link{lik.prof}} }
\examples{
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
