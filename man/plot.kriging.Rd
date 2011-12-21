\name{plot.kriging}
\alias{plot.kriging}
\title{Plot Semivariogram}
\description{
Plots distance versus semivariance with a fitted curve indicating the model used.
}
\usage{
\method{plot}{kriging}(x, main = "Semivariogram", xlab = "Distance", ylab = "Semivariance", \dots)
}
\arguments{
  \item{x}{object of class \code{kriging}.}
  \item{main}{See \link{par}.}
  \item{xlab}{See \link{par}.}
  \item{ylab}{See \link{par}.}
  \item{\dots}{arguments, passed to plot.default.} 
}
\references{

}
\author{
Omar E. Olmedo
}
\seealso{
  \code{\link{kriging}.}
}
\keyword{dplot}
\keyword{hplot}
