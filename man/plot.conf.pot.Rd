\name{plot.conf.pot}
\alias{plot.conf.pot}
\title{ Plot confidence bounds for GPD quantiles }
\description{
  Add confidence bounds in a pre-existing return value graph of a 
  peak over threshold analysis. 
}
\usage{
plot.conf.pot(pot.conf.obj, probs=c(0.025,0.975), nspline=3, ...)
}
\arguments{
  \item{pot.conf.obj}{ A confidence bounds object calculated in
          \code{conf.bounds.pot} }
  \item{probs}{ Probabilities for which the confidence lines
          shall be drawn. Note that only those probs can be drawn
          which were calculated in \code{conf.bounds.pot}. }
  \item{nspline}{ The number of subdivisions of the original data points
          to be calculated with spline (see Details). The total length of
          the data points used in plotting is \code{nspline*length(ret)}. }
  \item{\dots}{ Additional arguments for \code{lines}, such as line
          width \code{lwd}, or color \code{col}. }
}
\details{
  For a nice plot, it is necessary to evaluate the confidence bounds
  in sufficient resolution. This requires specifying a high number of
  return periods (argument \code{ret}) in \code{conf.bounds.pot}. 
}
\value{
  No value returned. Just lines add to a plot.
}
\references{ }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{conf.bounds.pot}}, \code{\link{plot.pot}} }
\examples{
data(engelberg.dprec)
threshold <- 30
pks <- extract.peaks(engelberg.dprec,threshold=threshold,
                     method="min.sep",min.sep=5)
npks <- length(pks[,1])
nyrs <- as.numeric(substr(
                   engelberg.dprec[dim(engelberg.dprec)[1],"date"],1,4))-
             as.numeric(substr(engelberg.dprec[1,"date"],1,4))
dat.pot <- fitGPD(pks,lamda=npks/nyrs,chi=threshold)
plot.pot(dat.pot,pch=20,cex.axis=1.3,cex=1.2,
          lwd=2,col.fit="blue",ylab="mm")

# calculate and plot maximum likelihood confidence intervals
tt.seq <- exp(seq(log(2),log(150),length=20))
conf.mle <- conf.bounds.pot(dat.pot,ret=tt.seq,cal="MLE")
plot.conf.pot(conf.mle,lwd=2,col="red")

# calculate and overplot likelihood profile confidence
conf.lprof <- conf.bounds.pot(dat.pot,ret=tt.seq,cal="LPROF")
plot.conf.pot(conf.lprof,lwd=2,col="green")

}
\keyword{ hplot }% at least one, from doc/KEYWORDS
