\name{gev.lik}
\alias{gev.lik}
\alias{gumbel.lik}
\alias{gev.lik.retv}
\alias{gumbel.lik.retv}
\title{ Likelihood function of GEV and GUMBEL}
\description{
  The negative log likelihood function for the Generalized Extreme 
  Value Distribution \code{gev.lik} and the Gumbel Distribution 
  \code{gumbel.lik}.
  These functions are used to calculate maximum likelihood estimates
  and likelihood profiles for the 
  parameters (\code{gev.lik}, \code{gumbel.lik}) and for the return values 
  (\code{gev.lik.retv}). 

  \code{gev.lik.retv} and \code{gumbel.lik.retv} are reparameterisations of
  the likelihood function using return value (at a predefined return period
  \code{ret}), scale and shape (the latter only with gev) to specify 
  the distribution.

  These functions are usually not called by the user directly.
}
\usage{
gev.lik(pars, data)
gumbel.lik(pars, data)
gev.lik.retv(pars, data, ret)
gumbel.lik.retv(pars, data, ret)
}
\arguments{
  \item{pars}{ A vector with the parameter values. For \code{gev.lik} this
      is the position parameter \code{chi}, the scale parameter \code{alpha}
      and the shape parameter \code{k} (in this order). 

      For \code{gumbel.lik} this is a two dimensional vector only with
      parameters \code{chi} and \code{alpha}.

      For \code{gev.lik.retv}
      \code{pars} is return value, scale and shape parameters (in this 
      order).}
  \item{data}{ A vector with the data sample. }
  \item{ret}{ A scalar number giving the return period (in units of the
      block size !!! ). }
}
\details{
   
}
\value{
  A real value giving the negative log likelihood.
}
\references{ see e.g. Coles 2001 }
\author{ Christoph Frei }
\note{  }
\seealso{ \code{\link{conf.bounds.xval}}, \code{\link{lik.prof}} }
\examples{
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
