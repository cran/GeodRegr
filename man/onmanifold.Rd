% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeodRegr.R
\name{onmanifold}
\alias{onmanifold}
\title{Manifold check and projection}
\usage{
onmanifold(manifold, y)
}
\arguments{
\item{manifold}{Type of manifold (\code{'euclidean'}, \code{'sphere'},
\code{'hyperbolic'}, or \code{'kendall'}).}

\item{y}{A vector, matrix, or data frame whose columns should represent
points on the manifold.}
}
\value{
A named list containing \item{on}{a logical vector describing whether
  or not each column of \code{y} is on the manifold.} \item{data}{a matrix of
  data frame of the same dimensions as \code{y}; each column of \code{y} has
  been projected onto the manifold.}
}
\description{
Checks whether each data point in \eqn{y} is on the given manifold, and if
not, provides a modified version of \eqn{y} where each column has been
projected onto the manifold.
}
\examples{
y1 <- matrix(rnorm(10), ncol = 2)
y1 <- y1[, 1] + (1i) * y1[, 2]
y2 <- matrix(rnorm(10), ncol = 2)
y2 <- y2[, 1] + (1i) * y2[, 2]
y3 <- matrix(rnorm(10), ncol = 2)
y3 <- y3[, 1] + (1i) * y3[, 2]
y3 <- (y3 - mean(y3)) / norm(y3 - mean(y3), type = '2') # project onto preshape space
y <- matrix(c(y1, y2, y3), ncol = 3)
onmanifold('kendall', y)

}
\author{
Ha-Young Shin
}
