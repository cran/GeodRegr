# Upper regularized gamma function
Q <- function(s, x) {
  return(zipfR::Rgamma(s, x, lower = F))
}

# Inverse of the upper regularized gamma function
Qinv <- function(s, x) {
  return(zipfR::Rgamma.inv(s, x, lower = F))
}

# Lower regularzied gamma function
P <- function(s, x) {
  return(zipfR::Rgamma(s, x, lower = T))
}

# Inverse of the lower regularized gamma function
Pinv <- function(s, x) {
  return(zipfR::Rgamma.inv(s, x, lower = T))
}

# Upper incomplete gamma function
uppergamma <- function(s, x) {
  return(zipfR::Igamma(s, x, lower = F))
}

# Lower incomplete gamma function
lowergamma <- function(s, x) {
  return(zipfR::Igamma(s, x, lower = T))
}

#' Approximate ARE of an M-type estimator to the least-squares estimator
#'
#' Approximate asymptotic relative efficiency (ARE) of an M-type estimator to
#' the least-squares estimator given Gaussian errors, calculated using a tangent
#' space approximation.
#'
#' @param estimator M-type estimator (\code{'l2'}, \code{'l1'}, \code{'huber'},
#'   or \code{'tukey'}).
#' @param k Dimension of the manifold.
#' @param c A positive multiplier, or a vector of them, of \eqn{\sigma}, the
#'   square root of the variance, used in the cutoff parameter for the
#'   \code{'huber'} and \code{'tukey'} estimators; should be \code{NULL} for the
#'   \code{'l2'} or \code{'l1'} estimators.
#' @return Approximate ARE
#' @references Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression.
#' <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @seealso \code{\link{are_nr}}.
#' @examples
#' are('l1', 10)
#'
#' @export
are <- function(estimator, k, c = NULL) {
  if (((estimator == 'huber') | (estimator == 'tukey')) & is.null(c)) {
    stop('a c value must be provided if the chosen m-estimator is huber or tukey')
  }
  if (!is.null(c)) {
    if ((estimator == 'l2') | (estimator == 'l1')) {
      warning('l2 and l1 do not use a c value')
    }
    if (any(c <= 0)) {
      stop('c must be positive')
    }
  }
  if ((k %% 1 != 0) | (k < 1)) {
    stop('k must be a positive integer')
  }
  if (estimator == 'l2') {
    result <- 1
  } else if (estimator == 'l1') {
    result <- ((gamma((k + 1) / 2))^2) / ((gamma(k / 2) * gamma((k + 2) / 2)))
  } else if (estimator == 'huber') {
    if (k == 1) { # k == 1 case must be treated separately
      numfactor <- (k / 2) * lowergamma(k / 2, 0.5 * c^2)
    } else {
      numfactor <- (k / 2) * lowergamma(k / 2, 0.5 * c^2) + c * (k - 1) * (2^-1.5) * uppergamma((k - 1) / 2, 0.5 * c^2)
    }
    denfactor <- lowergamma((k + 2) / 2, 0.5 * c^2) + (0.5 * c^2) * uppergamma(k / 2, 0.5 * c^2)
    result <- (numfactor^2) / (gamma((k + 2) / 2) * denfactor)
  } else if (estimator == 'tukey') {
    numfactor <- (2 * (k + 4) / (c^4)) * lowergamma((k + 4) / 2, 0.5 * c^2) - (2 * (k + 2) / (c^2)) * lowergamma((k + 2) / 2, 0.5 * c^2) + (k / 2) * lowergamma(k / 2, 0.5 * c^2)
    denfactor <- lowergamma((k + 2) / 2, 0.5 * c^2) - (8 / (c^2)) * lowergamma((k + 4) / 2, 0.5 * c^2) + (24 / (c^4)) * lowergamma((k + 6) / 2, 0.5 * c^2) - (32 / (c^6)) * lowergamma((k + 8) / 2, 0.5 * c^2) + (16 / (c^8)) * lowergamma((k + 10) / 2, 0.5 * c^2)
    result <-  (numfactor^2) / (gamma((k + 2) / 2) * denfactor)
  } else {
    stop('the M-estimator must be one of l2, l1, huber, or tukey')
  }
  return(result)
}

# Partial derivative of the are function with respect to c
deriv <- function(estimator, k, c) {
  if (estimator == 'huber') {
    if (k == 1) { # k == 1 case must be treated separately
      factor1 <- (k / 2) * lowergamma(k / 2, 0.5 * c^2)
      factor3 <- (c^(k - 1)) * (2^(-k / 2)) * exp(-0.5 * c^2)
    } else {
      factor1 <- (k / 2) * lowergamma(k / 2, 0.5 * c^2) + c * (k - 1) * (2^-1.5) * uppergamma((k - 1) / 2, 0.5 * c^2)
      factor3 <- (c^(k - 1)) * (2^(-k / 2)) * exp(-0.5 * c^2) + (k - 1) * (2^-1.5) * uppergamma((k - 1) / 2, 0.5 * c^2)
    }
    factor2 <- lowergamma((k + 2) / 2, 0.5 * c^2) + (0.5 * c^2) * uppergamma(k / 2, 0.5 * c^2)
    factor4 <- c * uppergamma(k / 2, 0.5 * c^2)
  } else if (estimator == 'tukey') {
    factor1 <- (2 * (k + 4) / (c^4)) * lowergamma((k + 4) / 2, 0.5 * c^2) - (2 * (k + 2) / (c^2)) * lowergamma((k + 2) / 2, 0.5 * c^2) + (k / 2) * lowergamma(k / 2, 0.5 * c^2)
    factor2 <- lowergamma((k + 2) / 2, 0.5 * c^2) - (8 / (c^2)) * lowergamma((k + 4) / 2, 0.5 * c^2) + (24 / (c^4)) * lowergamma((k + 6) / 2, 0.5 * c^2) - (32 / (c^6)) * lowergamma((k + 8) / 2, 0.5 * c^2) + (16 / (c^8)) * lowergamma((k + 10) / 2, 0.5 * c^2)
    factor3 <- -(8 * (k + 4) / (c^5)) * lowergamma((k + 4) / 2, 0.5 * c^2) + (4 * (k + 2) / (c^3)) * lowergamma((k + 2) / 2, 0.5 * c^2) - c^(k - 1) * 2^(-((k - 2) / 2)) * exp(-0.5 * c^2)
    factor4 <- (16 / (c^3)) * lowergamma((k + 4) / 2, 0.5 * c^2) - (96 / (c^5)) * lowergamma((k + 6) / 2, 0.5 * c^2) + (192 / (c^7)) * lowergamma((k + 8) / 2, 0.5 * c^2) - (128 / (c^9)) * lowergamma((k + 10) / 2, 0.5 * c^2)
  } else {
    stop('the M-estimator must be one of huber or tukey')
  }
  numerator <- 2 * factor1 * factor3 * factor2 - (factor1^2) * factor4
  denominator <- (gamma((k + 2) / 2)) * factor2^2
  result <- numerator / denominator
  return(result)
}

#' Newton-Raphson method for the \code{are} function
#'
#' Finds the positive multiplier of \eqn{\sigma}, the square root of the
#' variance, used in the cutoff parameter that will give the desired
#' (approximate) level of efficiency for the provided M-type estimator. Does so
#' by using \code{are} and its partial derivative with respect to \code{c} in
#' the Newton-Raphson method.
#'
#' As is often the case with the Newton-Raphson method, the starting point must
#' be chosen carefully in order to ensure convergence. The use of the graph of
#' the \code{are} function to find a starting point close to the root is
#' recommended.
#'
#' @param estimator M-type estimator (\code{'huber'} or \code{'tukey'}).
#' @param k Dimension of the manifold.
#' @param startingpoint Initial estimate for the Newton-Raphson method. May be
#'   determined after looking at a graph of the \code{are} function.
#' @param level The desired ARE to the \code{'l2'} estimator.
#' @return Positive multiplier of \eqn{\sigma}, the square root of the variance,
#'   used in the cutoff parameter, to give the desired level of efficiency.
#' @references
#' Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression. <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @seealso \code{\link{are}}.
#' @examples
#' dimension <- 4
#' x <- 1:10000 / 1000
#' # use a graph of the are function to pick a good starting point
#' plot(x, are('huber', dimension, x) - 0.95)
#' are_nr('huber', dimension, 2)
#'
#' @export
are_nr <- function(estimator, k, startingpoint, level = 0.95) {
  if (estimator == 'huber') {
    if (level <= are('l1', k)) {
      stop('the ARE of the L1 estimator is greater than the proposed ARE level')
    }
  }
  if  ((level >= 1) | (level <= 0)) {
    stop('the proposed are level is invalid')
  }
  if (startingpoint < 0) {
    stop('the starting point cannot be negative')
  }
  old_c <- startingpoint + 5
  new_c <- startingpoint
  count <- 0
  while (abs(new_c - old_c) > 0.000001) {
    old_c <- new_c
    new_c <- old_c - (are(estimator, k, old_c) - level) / deriv(estimator, k, old_c)
    if (new_c < 0) {
      stop('try again with a new starting point; consider using a graph of the are function to pick a good starting point')
    }
    count <- count + 1
    if (count > 1000) {
      stop('try again with a new starting point; consider using a graph of the are function to pick a good starting point')
    }
  }
  return(new_c)
}
