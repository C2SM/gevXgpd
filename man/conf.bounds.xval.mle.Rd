\name{conf.bounds.xval.mle}
\alias{conf.bounds.xval.mle}
\alias{conf.bounds.xval.sim}
\alias{conf.bounds.xval.lprof}
\title{ Confidence intervals for GEV return values }
\description{
  Functions for estimating confidence intervals of return values in
  a block maximum extreme value analysis. These
  are the functions for individual methods and should not be used by
  the user. The user should invoke \code{conf.bounds.xval} instead.
}
\usage{
conf.bounds.xval.mle(xval, ret = c(5, 10, 50, 100, 500)/xval$lamda, 
                     probs = c(0.025, 0.975))

conf.bounds.xval.sim(xval, ret = c(5, 10, 50, 100, 500)/xval$lamda, 
                     probs = c(0.025, 0.975), size = 500)

conf.bounds.xval.lprof(xval, ret = c(5, 10, 50, 100, 500)/xval$lamda, 
                     probs = c(0.025, 0.975), profs.out = FALSE)
}
\arguments{
  \item{xval}{ An object from \code{fitGEV} containing the data and
         fitted GEV distribution. }
  \item{ret}{ A vector of return periods for which the confidence 
         intervals shall be calculated. }
  \item{probs}{ A vector of probs specifying the quantiles of the 
         confidence intervals. Default is the 95% confidence interval. 
         It is possible (recommended) to calculate several confidence
         bounds at once, e.g. probs = c(0.025,0.25,0.5,0.75,0.975) since
         this is more efficient than several calls to this function.}
  \item{size}{ The number of samples to simulate when \code{cal="SIM"} i.e. 
         the interval is calculated by parameteric resampling. }
  \item{profs.out}{ TRUE or FALSE, if likelihood profiles should be provided
         in the output. Is only relevant for \code{cal="LPROF"}.}
}
\details{
  These
  are the functions for individual methods and should not be used by
  the user. The user should invoke \code{conf.bounds.xval} instead.
}
\value{

}
\references{  }
\author{ Christoph Frei }
\note{ 

}
\seealso{ \code{\link{conf.bounds.xval}} }
\examples{

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
