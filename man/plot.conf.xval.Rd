\name{plot.conf.xval}
\alias{plot.conf.xval}
\title{ Plot confidence bounds for GEV analysis }
\description{
  Adds confidence bounds of the return values into a 
  pre-existing Gumbel diagram. The confidence bounds need be calculated
  first using \code{conf.bounds.xval}. 
}
\usage{
plot.conf.xval(xval.conf.obj, probs=c(0.025,0.975), nspline=3, ...)
}
\arguments{
  \item{xval.conf.obj}{ An object containing the confidence bounds to be
        be plotted. This object can be calculated using 
        \code{conf.bounds.xval}. }
  \item{probs}{ Probabilities for which the confidence lines
          shall be drawn. Note that only those probs can be drawn
          which were calculated in \code{conf.bounds.xval}. }
  \item{nspline}{ The number of subdivisions of the original data points
          to be calculated with spline (see Details). The total length of
          the data points used in plotting is \code{nspline*length(ret)}. }
  \item{\dots}{ Additional graphic arguments passed on to the command
        \code{lines}. E.g line widths and colors.}
}
\details{
  For a smooth apparence of the confidence bounds (relevant in particular
  if calculated by resampling and profile likelihood), the values of the
  quantiles are interpolated between return periods using splines. It is
  therefore not necessary to specify in \code{conf.bounds.xval} very 
  many return periods. A well distributed set of about 10-20 will do.
}
\value{
  No value returned. Just an addition to pre-existing plot.}
\references{  }
\author{ Christoph Frei }
\note{ 
}
\seealso{ \code{\link{conf.bounds.xval}}, \code{\link{plot.xval}} }
\examples{
data(engelberg.annmax)
xval.mlik <- fitGEV(engelberg.annmax,estim="mlik")

# the usual plot with data and fit
plot.xval(xval.mlik,pch=20,cex.axis=1.3,cex=1.2,
          lwd=2,col.fit="blue",ylab="mm",tlim=c(1.1,150))

# calculate the confidence interval
conf.xval.mlik <- conf.bounds.xval(xval.mlik,
                      ret=c(2,3,4,5,7,10,12,15,20,50,70,100,120,150),
                                   cal="MLE",
                                   probs=c(0.025,0.25,0.75,0.975))

# alternative using likelihood profile
conf.xval.lprof <- conf.bounds.xval(xval.mlik,
                      ret=c(2,3,4,5,7,10,12,15,20,50,70,100,120,150),
                                   cal="LPROF",
                                   probs=c(0.025,0.25,0.75,0.975))

# finally plot the confidence interval
plot.conf.xval(conf.xval.mlik,lwd=2,col="red")

# compare to likelihood profile
plot.conf.xval(conf.xval.lprof,lwd=2,col="green")

}
\keyword{ aplot }% at least one, from doc/KEYWORDS
