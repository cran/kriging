\name{image.kriging}
\alias{image.kriging}
\title{Map kriged data}
\description{
Create maps using the coordinates and predicted values in objects of class \code{kriging}.
}
\usage{
  \method{image}{kriging}(x, main = NULL, xlab = "", ylab = "", col = heat.colors(100), \dots)
}
\arguments{
  \item{x}{object of class \code{kriging}.}
  \item{main}{See \link{par}.}
  \item{xlab}{See \link{par}.}
  \item{ylab}{See \link{par}.}
  \item{col}{See \link{par}.}
  \item{\dots}{arguments, passed to image.default.}
}
\author{
Omar E. Olmedo
}
\seealso{
  \code{\link{kriging}.}
}
\keyword{dplot}
\keyword{hplot}
