\name{fitGEV.lmom}
\alias{fitGEV.lmom}
\alias{fitGEV.mlik}
\alias{fitGEV.pmlik}
\alias{fitGUMBEL.lmom}
\alias{fitGUMBEL.mlik}
\title{ Fit GEV/GUMBEL distribution }
\description{
  The functions for fitting GEV or GUMBEL distributions to data. These
  are the functions for individual methods and should not be used by
  the user. The user should invoke \code{fitGEV} instead.
}
\usage{
fitGEV.lmom(data.df, lamda = 1, ret = c(5, 10, 50, 100, 500)/lamda, ...)
fitGEV.mlik(data.df, lamda = 1, ret = c(5, 10, 50, 100, 500)/lamda, 
         method = "Nelder-Mead", maxit = 10000, ...)
fitGEV.pmlik(data.df, lamda = 1, ret = c(5, 10, 50, 100, 500)/lamda, 
         method = "Nelder-Mead", maxit = 10000, ...)
fitGUMBEL.lmom(data.df, lamda = 1, ret = c(5, 10, 50, 100, 500)/lamda, ...)
fitGUMBEL.mlik(data.df, lamda = 1, ret = c(5, 10, 50, 100, 500)/lamda, 
         method = "Nelder-Mead", maxit = 10000, ...)
}
\arguments{
  \item{data.df}{ A data frame with the data to be fitted in the first 
     column. The data frame
     can have more than one column. Additional columns are not used.}
  \item{lamda}{ The number of observations (blocks) per time unit. 
     Values other than
     lamda=1 can be used to express the return period in other units than 
     in units of the block size.}
  \item{ret}{ A vector of return periods for which fitted return levels 
     should be provided in the resulting list (element \code{fit}.}
  \item{method}{ A character string specifying the method of numerical 
     optimisation to be used in the maximum likelihood (\code{estim="mlik"}) 
     and the modified maximum likelihood (\code{estim="pmlik"}) estimations. 
     This parameter is passed on to optim(). }
  \item{maxit}{ The maximum number of iterative steps in optim() for 
     maximum likelihood (\code{estim="mlik"}) and modified maximum 
     likelihood (\code{estim="pmlik"}) estimation. }
  \item{\dots}{ Placeholder for parameters not used by the method. }
}
\details{
  These are the functions for individual methods and should not be used by
  the user. The user should invoke \code{fitGEV} instead.
}
\value{

}
\references{  }
\author{ Christoph Frei }
\note{ 
}
\seealso{ \code{\link{fitGEV}} }
\examples{
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
