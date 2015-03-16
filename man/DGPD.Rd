\name{DGPD}
\alias{DGPD}
\alias{FGPD}
\alias{XGPD}
\alias{RGPD}
\title{ The Generalized Pareto Distribution }
\description{
  Density (DGPD), Cumulative Frequency (FGPD), 
  Quantile Function (XGPD) and Random
  number generator (RGPD) for the Generalized Pareto Distribution.
}
\usage{
DGPD(xx, alpha, chi=0, k)
FGPD(xx, alpha, chi=0, k)
XGPD(ff, alpha, chi=0, k)
RGPD(alpha, chi=0, k, n)
}
\arguments{
  \item{xx}{ A vector of values (quantiles) for which to calculate 
             the density }
  \item{ff}{ A vector of cumulative frequencies for which to calculate 
             the quantiles }
  \item{alpha}{ The scale parameter of the GPD}
  \item{chi}{ The sekected threshold. This is not an estimated 
       parameter but can be chosen differently from 0 to express
       exceedances in terms of original variables.}
  \item{k}{ The shape parameter of the GPD}
  \item{n}{ The number of random samples to be generated}
}
\details{
  Here we have used a notation similar to Zwiers and Kharin 1998. 
  PLEASE NOTE: The present parameters \code{chi}, 
  \code{alpha}, \code{k} correspond to the parameters u, sigma 
  and -chi (minus chi) respectively in the notation of Coles (2001).
  
  In case of a very small value of k (abs(k) < 0.000001) 
  the Exponential distribution is chosen.
}
\value{
  A vector of values. 
}
\references{ see e.g. formulae in Palutikof et al. 1999, Coles 2001 }
\author{ Christoph Frei }
\note{ }
\seealso{ \code{\link{Y.of.T}}, \code{\link{F.of.T}} }
\examples{

# a histogram of random Gumbel variates
hist(RGPD(alpha=1,chi=0.0,k=0.0,n=1000))

# a density plot 
xx <- seq(0,4,length=200)
plot(x=xx,y=DGPD(xx=xx,alpha=1.0,chi=0.0,k=0.0),type="l",
     ylim=c(0,1.0),xlab="x",ylab="density",main="GPD",lwd=2,col="green")
lines(x=xx,y=DGPD(xx=xx,alpha=1.0,chi=0.0,k=0.2),lwd=2,col="red")
lines(x=xx,y=DGPD(xx=xx,alpha=1.0,chi=0.0,k=-0.2),lwd=2,col="blue")

# a quantile plot as function of Exponential variate V
VV <- seq(0,5,length=200)
F.of.V <- function(v) {1-exp(-v)}
V.of.F <- function(f) {-log(1-f)}
XX1 <- XGPD(F.of.V(VV),alpha=1.0,chi=0.0,k=0.0)
XX2 <- XGPD(F.of.V(VV),alpha=1.0,chi=0.0,k=0.2)
XX3 <- XGPD(F.of.V(VV),alpha=1.0,chi=0.0,k=-0.2)
plot(x=VV,y=XX1,type="l",
     ylim=c(0,6),xlab="Exponential Variate V=-log(1-F)",
     ylab="Exceedance Quantile Y",main="GPD",lwd=2,col="green")
lines(x=VV,y=XX2,lwd=2,col="red")
lines(x=VV,y=XX3,lwd=2,col="blue")
TT <- c(5,10,20,50,100)
axis(3,at=V.of.F(F.of.T(TT)),labels=TT)

}
\keyword{ distribution }% at least one, from doc/KEYWORDS

