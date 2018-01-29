#' @title Augmented Dickey-Fuller Statistic
#'
#' @description Calculate the Augmented Dickey-Fuller Statistic with a fixed lag order .
#'
#' @param y the time series to be used.
#' @param adflag is the lag order.
#' @param mflag 1 for ADF with constant and whithout trend, 2 for ADF with constant and trend and 3 for ADF without constant and trend.
#' @keywords AugmentedDickey-FullerTest.
#' @references Phillips, P.C. & Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500". \emph{SSRN Electronic Journal}.
#' @import stats
#' @import foreach
#' @import MASS
#' @export
#' @examples
#' y <- rnorm(10)
#' ADF_FL(y, adflag = 1, mflag = 2)

ADF_FL <- function(y, adflag = 0, mflag = 1){
  t1 <- length(y) - 1
  const <- matrix(1, nrow = t1, ncol = 1)
  trend <- matrix(seq(from = 1, to = t1, by = 1))

  y1 <- matrix(y[(length(y) - t1):(length(y) - 1)])
  dy <- matrix(y[2:(length(y))] - y[1:(length(y) - 1)])
  dy0 <- matrix(dy[(length(dy) - t1 + 1):(length(dy))])

  if (mflag == 1) {
    x <- cbind(y1, const)
  } else if (mflag == 2) {
    x <- cbind(y1, const, trend)
  } else if (mflag == 3) {
    x <- y1
  }

  m0 <- nrow(x)
  t2 <- t1 - adflag
  x20 <- x[(m0-t2+1):m0,]
  n1 <- nrow(x20)
  m1 <- ncol(x20)
  dy01 <- matrix(dy0[(length(dy0) - t2 + 1):length(dy0)])

  if (adflag > 0){
    x2 <- matrix(0, nrow = n1, ncol = (m1 + adflag))
    x2[, 1:m1] <- x20

    for(j in 1:adflag){
      x2[, (m1 + j)] <- dy[(nrow(dy)-t2+1-j):(nrow(dy)-j)]
    }
    } else if (adflag == 0){
    x2 <- x20
    }

  beta <- ginv(t(x2) %*% x2) %*% (t(x2) %*% dy01)
  eps <- as.matrix(dy01) - x2 %*% beta

  if (mflag == 1){
    dof <- t2 - adflag - 2
  }
  if (mflag == 2){
    dof <- t2 - adflag - 3
  }
  if (mflag == 3){
    dof <- t2 - adflag - 1
  }

  c <- t(eps) %*% eps/dof
  A <- t(x2) %*% x2
  sig <- sqrt(diag(as.numeric(c) * solve(A)))
  sig <- as.matrix(sig)
  tvalue <- beta/sig
  estm <- tvalue[1]
  result <- list('Augmented Dickey-Fuller Statistic' = estm)
  return(result)
}

