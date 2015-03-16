\name{pretty.logscale}
\alias{pretty.logscale}
\title{ Pretty numbers on a log scale }
\description{
  Calculates a vector of "pretty" values that can be used as
  labels on a log-scale axis. The function is very similar to
  the built-in function pretty() which determines "pretty"
  equidistant labels.
}
\usage{
pretty.logscale(x, ...)
}
\arguments{
  \item{x}{ A vector of values (at least two) indicating the range
      over which the labels should go. }
  \item{\dots}{ Arguments passed on to pretty(), such as \code{n} giving
      the approximate number labels required. }
}
\details{
  
}
\value{
  A vector of "pretty" labels as numbers.
}
\references{ }
\author{ Christoph Frei }
\note{ 

}
\seealso{ \code{\link{pretty}} }
\examples{
pretty.logscale(x=c(1,200))

}
\keyword{ dplot }% at least one, from doc/KEYWORDS
