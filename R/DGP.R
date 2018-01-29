#' @title Random walk.
#'
#' @description Generate a random walk with drift 1/n.
#'
#' @param n sample size. Number of rows in the generated matrix.
#' @param niter number of columns in the generated matrix.
#' @keywords random walk generation.
#' @export
#'
#' @examples
#' DGP(n = 100, niter = 10)

DGP <- function(n, niter){
  SI <- niter
  set.seed(101)
  rnorm(SI)
  u0 <- 1/n
  rn <- replicate(niter, rnorm(n))
  z <- rn + u0
  y <- apply(z, 2, cumsum)
  return(y)
}
