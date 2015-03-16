\name{Y.of.T}
\alias{Y.of.T}
\alias{T.of.Y}
\alias{F.of.T}
\alias{T.of.F}
\alias{Y.of.F}
\alias{F.of.Y}
\title{ Transformations for the GEV }
\description{
  Transformations between the Gumbel Variate Y, the return period T,
  and the cumulative frequency F of the Generalized Extreme Value
  Distribution.
}
\usage{
Y.of.T(t, lamda = 1)
T.of.Y(y, lamda = 1)
F.of.T(t, lamda = 1)
T.of.F(f, lamda = 1)
Y.of.F(f)
F.of.Y(y)
}
\arguments{
  \item{t}{ The return period in units of the block size. }
  \item{y}{ The Gumbel variate }
  \item{f}{ The cumulative frequency }
  \item{lamda}{ Expected number of data points per block. Values for 
      \code{lamda} different from 1 can be used to express the return period
      in other units than in units of the block size. E.g. if the block 
      size is half a year but the return period should be specified in
      units of years then one would have to set \code{lamda}=2.}
}
\details{
  These functions are rarely needed by the user directly.
  \code{Y.of.T} transforms return period into Gumbel Variate, 
  \code{F.of.T} transforms return period into cumulate frequency, 
  etc.
  Input is be either a single number or a vector.
}
\value{
  A single number of a vector, depending on the input.
}
\references{ }
\author{ Christoph Frei }
\note{ }
\seealso{ \code{\link{XGEV}}, \code{\link{FGEV}} }
\examples{
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
