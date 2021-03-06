\name{plot.pot}
\alias{plot.pot}
\title{ Plot a Pareto diagram }
\description{
  Plot a Pareto diagram (return value as a function of return period) for
  a GPD fit to a dataset of threshold exceedances. The diagram depicts 
  the data (using conventional Weibull plotting points F=rank(x)/(N+1)) 
  and a precalculated
  fit of a GPD to that data. In a second step it is possible to draw 
  confidence intervals. The Pareto diagram uses the Pareto variate as 
  a linear scale for recurrence. Hence the return period is log scale.
}
\usage{
plot.pot(fit.pot, tlim=NULL, tat, tlab="T", ylim=NULL, ylab, 
         col = "black", col.fit=col, plot.data = TRUE, 
         label.data = FALSE, label.cex = 1, label.col = col, 
         label.adj = c(-0.2, -0.2), cex.axis, 
         main="Peaks over Threshold",add=FALSE,...)
}
\arguments{
  \item{fit.pot}{ A data object from fitGPD(). Encompasses the data and 
       fit. A list with elements as described in \code{fitGPD}. }
  \item{tlim}{ The limits of the return period axis (horizontal). If not
       specified \code{tlim} is set such that all data points are contained. }
  \item{tat}{ The values of return period to label on return period axis. 
       Default is "nice values". }
  \item{tlab}{ The label to use for return period axis. Default is "T". }
  \item{ylim}{ The plot range
       for the vertical (i.e. the return value) axis. \code{ylim=NULL} leaves
       the range unspecified. It will then be determined automatically 
       by plot. }
  \item{ylab}{ The label to use for return value axis.  Default is
       the name of the first column in fit.xval.}
  \item{col}{ The color to use for the data points. }
  \item{col.fit}{ The color to use for the GEV fit. Default is: Same as 
       \code{col}. }
  \item{plot.data}{ Whether to plot data points: FALSE or TRUE. Default 
       is TRUE. }
  \item{label.data}{ Whether to label the data points on the graph. The 
       labels are taken from column "lab" in the original data frame. If 
       this column does not exist the plotting fails 
       if \code{label.data}=TRUE . Default is to not plot labels. }
  \item{label.cex}{ The magnification to be used for the plotting labels
       in the graph. }
  \item{label.col}{ The color to be used for the plotting labels
       in the graph.  }
  \item{label.adj}{ The position to be used for the plotting labels
       in the graph. See help(par) under \code{adj} for more details. }
  \item{cex.axis}{ The magnification to be used for the axis annotation. }
  \item{main}{ Title of the plot. }
  \item{add}{ If the points and fit should be added to an already 
       existing plot. Allows to overplot the results from several
       fits, e.g. from the same dataset with varying block size.}
  \item{\dots}{ Additional arguments that are passed to \code{plot()}. 
       Frequently used parameters are: \code{pch}: symbol type 
       (see help of \code{points}), \code{cex}: symbol size, \code{cex.axis}: 
       magnification of axis annotation, \code{lty} and \code{lwd}: line
       type and color of the fit. }
}
\details{
  The axis at the upper margin of the diagram is a transformed variate V
  for the frequency. 
  Return period T is labelled at the bottom axis. The following relation
  holds between the two: V(x) = log(1/(1-F(x))) and T(x) =1/(1-F(x)), i.e.
  V(x)=log(T(x)), F is the cumulative frequency.
  If a fit is a straight line this indicates that the resulting fit is
  an Exponential distribution (i.e. shape=0). 
}
\value{
  No data object is returned. Only a plot.
}
\references{  }
\author{ Christoph Frei }
\note{  }
\seealso{ \code{\link{fitGPD}}, \code{\link{plot.conf.pot}}, 
          \code{\link{conf.bounds.pot}} }
\examples{
data(engelberg.dprec)
threshold <- 30
pks <- extract.peaks(engelberg.dprec,threshold=threshold,
                     method="min.sep",min.sep=5)
npks <- length(pks[,1])
nyrs <- as.numeric(substr(
                  engelberg.dprec[dim(engelberg.dprec)[1],"date"],1,4))-
             as.numeric(substr(engelberg.dprec[1,"date"],1,4)) + 1
dat.pot <- fitGPD(pks,lamda=npks/nyrs,chi=threshold)

# the usual plot with data and fit
plot.pot(dat.pot,pch=20,cex.axis=1.3,cex=1.2,
          lwd=2,col.fit="blue",ylab="mm")

# a plot with labels for the first 4 data points
xx <- pks[,1]
pks[rank(xx)<(length(xx)-3),"date"] <- NA
names(pks)[2] <- "lab"
dat.pot <- fitGPD(pks,lamda=npks/nyrs,chi=threshold)
plot.pot(dat.pot,pch=20,cex.axis=1.3,cex=1.2,
          lwd=2,col.fit="blue",ylab="mm",
          label.data=TRUE,label.adj=c(+1.05,-0.3))

# some additional text description in the graph
par.char <- "Extreme 1-day Precip"
stnam <- "Engelberg"
yy1 <- 1901
yy2 <- 2005
x.l <- 10^par("usr")[1]
x.r <- 10^par("usr")[2]
y.d <- par("usr")[3]
y.u <- par("usr")[4]
dy <- diff(par("usr")[c(3,4)])/15
text(par.char,x=x.r,
           y=y.d+3.1*dy,pos=2,cex=1.5)
text(stnam,x=x.r,
           y=y.d+2.0*dy,pos=2,cex=2.5)
text(paste(as.character(yy1),"-",as.character(yy2)),x=x.r,
           y=y.d+1.0*dy,pos=2,cex=1.5)

}
\keyword{ hplot }% at least one, from doc/KEYWORDS
