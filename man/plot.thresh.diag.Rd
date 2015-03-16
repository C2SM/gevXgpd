\name{plot.thresh.diag}
\alias{plot.thresh.diag}
\title{ Plot Diagnostics for GPD Threshold Choice}
\description{
  Produces a diagnostic plot for choosing appropriate thresholds in a
  peak-over-threshold analysis using GPD. Three types of plots can be
  made: (a) mean residual life plot (showing the threshold dependence
  of the mean of the threshold exceedance). (b) shape as a function of
  threshold. (c) modified scale as a function of threshold.
}
\usage{
plot.thresh.diag(diag, diagnostic = "mean.exceed", 
        print.exceeds = 3, cex.exceeds = 1, 
        xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{diag}{ A diagnostic object calculated by \code{thresh.diag}. }
  \item{diagnostic}{ A character string indicatin the type of diagnostic
         to plot. One of "mean.exceed", "shape", "scale.star". }
  \item{print.exceeds}{ If \code{print.exceeds=1} then the number 
         of exceedances for each threshold are printed at the upper 
         margin of the plot. If there are many thresholds it is better
         to print the number only for every \code{print.exceeds} threshold
         for space reasons. If \code{print.exceeds=0} no numbers are 
         printed. }
  \item{cex.exceeds}{ Size of the text for number exceedances. }
  \item{xlim}{ Threshold range to be plotted. }
  \item{ylim}{ y-range to be plotted. }
  \item{xlab}{ x-label. }
  \item{ylab}{ y-label. }
  \item{\dots}{ Additional arguments to be passed on to \code{plot}. 
        Used for plot color and line width. }
}
\details{
  For details of the calculation of the diagnostics see \code{thresh.diag}.
}
\value{
  A plot.}
\references{ See e.g. Coles 2001. }
\author{ Christoph Frei }
\note{  }

\seealso{ \code{\link{thresh.diag}} }
\examples{
data(engelberg.dprec)
threshold <- 5
pks <- extract.peaks(engelberg.dprec,threshold=threshold,
                     method="min.sep",min.sep=5)
threshs <- seq(5,85,by=2.5)
diag <- thresh.diag(pks,threshs)

# plot shape as function of threshold
X11(width=7,height=5)
plot.thresh.diag(diag,diagnostic="shape",lwd=2,col="blue",
                 cex.lab=1.2,cex.axis=1.2,ylim=c(-0.25,1.0),
                 xlab="Threshold (mm)")
axis(1,at=seq(10,80,by=10),labels=as.character(seq(10,80,by=10)),cex.axis=1.2)
box(lwd=1.5)

# plot modified scale as function of threshold
X11(width=7,height=5)
plot.thresh.diag(diag,diagnostic="scale.star",lwd=2,col="red",
                 cex.lab=1.2,cex.axis=1.2,ylim=c(-20,100),
                 xlab="Threshold (mm)")
axis(1,at=seq(10,80,by=10),labels=as.character(seq(10,80,by=10)),cex.axis=1.2)
box(lwd=1.5)

# mean exceedance as function of threshold. Mean residual life plot.
X11(width=7,height=7)
plot.thresh.diag(diag,diagnostic="mean.exceed",lwd=2,col="green",
                 cex.lab=1.2,cex.axis=1.2,xlab="Threshold (mm)",
                 ylab="Mean Exceedance (mm)")
axis(1,at=seq(10,80,by=10),labels=as.character(seq(10,80,by=10)),cex.axis=1.2)
box(lwd=1.5)

dev.off()
dev.off()
dev.off()

}
\keyword{ hplot }% at least one, from doc/KEYWORDS
