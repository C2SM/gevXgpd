\name{DGEV}
\alias{DGEV}
\alias{FGEV}
\alias{XGEV}
\alias{RGEV}
\title{ The Generalised Extreme Value Distribution }
\description{
  Density (DGEV), Cumulative Frequency (FGEV), 
  Quantile Function (XGEV) and Random
  number generator (RGEV) for the Generalised Extreme Value Distribution.
}
\usage{
DGEV(xx, alpha, chi, k)
FGEV(xx, alpha, chi, k)
XGEV(ff, alpha, chi, k)
RGEV(alpha, chi, k, n)
}
\arguments{
  \item{xx}{ A vector of values (quantiles) for which to calculate 
             the density }
  \item{ff}{ A vector of cumulative frequencies for which to calculate 
             the quantiles }
  \item{alpha}{ The scale parameter of the GEV}
  \item{chi}{ The location parameter of the GEV}
  \item{k}{ The shape parameter of the GEV}
  \item{n}{ The number of random samples to be generated}
}
\details{
  Here we have used the notation of Zwiers and Kharin 1998. 
  The present parameters \code{chi}, 
  \code{alpha}, \code{k} correspond to the parameters mu, sigma 
  and -chi (minus chi) respectively in the notation of Coles (2001).
  
  The functions distinguish between the Gumbel-, the Weibull- and 
  the Frechet distribution depending on the value of k, the shape
  parameter. In case of a very small value of k (abs(k) < 0.000001) 
  the Gumbel distribution is chosen.
}
\value{
  A vector of values. 
}
\references{ see e.g. formulae in Zwiers and Kharin 1998, Palutikof 1999 }
\author{ Christoph Frei }
\note{ }
\seealso{ \code{\link{Y.of.T}}, \code{\link{F.of.T}} }
\examples{

# a histogram of random Gumbel variates
hist(RGEV(alpha=1,chi=0.0,k=0.0,n=1000))

# a density plot (c.f. Zwiers and Kharin 1998)
xx <- seq(-2,5,length=200)
plot(x=xx,y=DGEV(xx=xx,alpha=1.0,chi=0.0,k=0.0),type="l",
     ylim=c(0,0.4),xlab="x",ylab="density",main="GEV",lwd=2,col="green")
lines(x=xx,y=DGEV(xx=xx,alpha=1.0,chi=0.0,k=0.2),lwd=2,col="red")
lines(x=xx,y=DGEV(xx=xx,alpha=1.0,chi=0.0,k=-0.2),lwd=2,col="blue")

# a quantile plot as function of Gumbel variate Y
# (c.f. Zwiers and Kharin 1998)
YY <- seq(-2,5,length=200)
XX1 <- XGEV(F.of.Y(YY),alpha=1.0,chi=0.0,k=0.0)
XX2 <- XGEV(F.of.Y(YY),alpha=1.0,chi=0.0,k=0.2)
XX3 <- XGEV(F.of.Y(YY),alpha=1.0,chi=0.0,k=-0.2)
plot(x=YY,y=XX1,type="l",
     ylim=c(-2.5,6),xlab="Gumbel Variate Y",
     ylab="Quantile X",main="GEV",lwd=2,col="green")
lines(x=YY,y=XX2,lwd=2,col="red")
lines(x=YY,y=XX3,lwd=2,col="blue")
TT <- c(5,10,20,50,100)
axis(3,at=Y.of.T(TT),labels=TT)

}
\keyword{ distribution }% at least one, from doc/KEYWORDS

