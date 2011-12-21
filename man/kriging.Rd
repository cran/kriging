\name{kriging}
\alias{kriging}
\title{Ordinary Kriging}
\description{
Simple and highly optimized ordinary kriging algorithm to plot geographical data
}
\usage{
kriging(x, y, response, model = "spherical", lags = 10, pixels = 100, polygons = NULL)
}
\arguments{
  \item{x}{vector of x-axis spatial points.}
  \item{y}{vector of y-axis spatial points.}
  \item{response}{vector of observed values.}
  \item{model}{specification of the variogram model. Choices are \code{"spherical"}, \code{"exponential"} or \code{"gaussian"}. Defaults to \code{"spherical"}.}
  \item{lags}{number of lags. Defaults to \code{10}.}
  \item{pixels}{maximum number of points along either axis. Defaults to \code{100}.}
  \item{polygons}{list of polygons used to grid predicted values on to. The default value of \code{NULL} automatically generates an evenly spaced out rectangular grid of points spanning the range of the data.}
}
\details{
The kriging algorithm assumes a minimum number of observations in order to fit the variogram model. 
}
\value{
An object of class \code{kriging} that inherits from \code{list} and is composed of:
  \item{model}{character; variogram model.}
  \item{nugget}{numeric; value of nugget parameter.}
  \item{range}{numeric; value of range parameter.}
  \item{sill}{numeric; value of sill parameter.}
  \item{map}{data.frame; contains the predicted values along with the coordinate covariates.}
  \item{semivariogram}{data.frame; contains the distance and semivariance values.}
}
\author{
Omar E. Olmedo
}
\seealso{
  \code{\link{image.kriging}, \link{plot.kriging}.}
}
\examples{
# Krige random data for a specified area using a list of polygons
library(maps)
usa <- map("usa", "main", plot = FALSE)
p <- list(data.frame(usa$x, usa$y))

# Create some random data
x <- runif(50, min(p[[1]][,1]), max(p[[1]][,1]))
y <- runif(50, min(p[[1]][,2]), max(p[[1]][,2]))
z <- rnorm(50)

# Krige and create the map
kriged <- kriging(x, y, z, polygons=p, pixels=300)
image(kriged, xlim = extendrange(x), ylim = extendrange(y))

}
\keyword{datagen}
\keyword{models}
