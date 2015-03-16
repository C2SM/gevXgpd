\name{conf.bounds.pot.sim}
\alias{conf.bounds.pot.sim}
\alias{conf.bounds.pot.mle}
\alias{conf.bounds.pot.lprof}
\title{ Confidence bounds for GPD analysis }
\description{
  Functions for estimating confidence intervals of return values 
  and parameters in a peak over threshold analysis. These
  are the functions for individual methods and should not be used by
  the user. The user should invoke \code{conf.bounds.pot} instead.
}
\usage{
conf.bounds.pot.sim(pot, ret = c(5, 10, 50, 100, 500), 
                    probs = c(0.025, 0.975), size = 500)

conf.bounds.pot.mle(pot, ret = c(5, 10, 50, 100, 500), 
                    probs = c(0.025, 0.975))

conf.bounds.pot.lprof(pot, ret = c(5, 10, 50, 100, 500), 
                    probs = c(0.025, 0.975), profs.out = FALSE)
}
\arguments{
  \item{pot}{ An object from \code{fitGPD} containing the data and
         fitted GPD distribution. }
  \item{ret}{ A vector of return periods for which the confidence 
         intervals shall be calculated. }
  \item{probs}{ A vector of probs specifying the quantiles of the 
         confidence intervals. Default is the 95% confidence interval. 
         It is possible (recommended) to calculate several confidence
         bounds at once, e.g. probs = c(0.025,0.25,0.5,0.75,0.975) }
  \item{size}{ The number of samples to simulate when \code{cal="SIM"} i.e. 
         the interval is calculated by parameteric resampling. }
  \item{profs.out}{ TRUE or FALSE, if likelihood profiles should be provided
         in the output. Is only relevant for \code{cal="LPROF"}.}
}
\details{
  These
  are the functions for individual methods and should not be used by
  the user. The user should invoke \code{conf.bounds.pot} instead.
}
\value{

}
\references{  }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{conf.bounds.pot}} }
\examples{
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
