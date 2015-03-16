\name{qqgev}
\alias{qqgev}
\title{ Q-Q-Plot for Generalized Extreme Value Distribution }
\description{
  Makes a Q-Q plot of a sample for the Generalized Extreme Value Distribution.
}
\usage{
qqgev(x, fit, main="GEV Q-Q Plot",
         xlab="Theoretical Quantiles",ylab="Empirical Quantiles",...)
}
\arguments{
  \item{x}{ A sample (numeric vector) for which to depict the Q-Q plot. }
  \item{fit}{ A data object from fit.GEV, or a list with entries \code{chi},
          \code{alpha}, \code{k}, i.e. estimated or stipulated parameters
          of the GEV to use.}
  \item{main}{ Title of the plot. A character. }
  \item{xlab}{ x-axis label. }
  \item{ylab}{ y-axis label. }
  \item{...}{ Arguments passed on to function \code{plot}.}
}
\details{
  Function \code{ppoints} is used for the definition plotting points.

  Note that the GEV fit is not calculated internally. A fit (GEV parameters)
  has to calculated (stipulated) prior to calling the function.
}
\value{
  No output. A plot only. 
}
\references{  }
\author{ Christoph Frei }
\note{ }
\seealso{ }
\examples{

# a simulated sample from the Gumbel distribution
dd <- RGEV(alpha=1,chi=0.0,k=0.0,n=100)

hh <- fitGEV(data.df=data.frame(dd),estim="mlik")

# a q-q plot 
qqgev(dd,fit=hh)
abline(a=0,b=1)

}
\keyword{ distribution }% at least one, from doc/KEYWORDS

