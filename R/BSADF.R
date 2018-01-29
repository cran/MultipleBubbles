#' @title Critical values for backward SADF statistic sequence.
#'
#' @description Calculate critical value sequences for the backward
#'  sup ADF statistic sequence using Monte Carlo simulations for a sample
#'  generated from a Normal distribution.
#'
#' @param m Number of Monte Carlo Simulations
#' @param t Sample size.
#' @param adflag is the lag order.
#' @param mflag 1 for ADF with constant and whithout trend, 2 for ADF with constant and trend and 3 for ADF without constant and trend.#' @keywords AugmentedDickey-FullerTest backwardSADF MonteCarlo.
#' @references Phillips, P.C. & Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500". \emph{SSRN Electronic Journal}.
#' @import stats
#' @import foreach
#' @import MASS
#' @export
#'
#' @examples
#' foo <- bsadf(m = 20, t = 50, adflag = 1, mflag = 2)
#' plot(foo$quantiles[2,], type = 'l')


bsadf <- function(m, t, adflag = 0, mflag = 1){
  qe <- as.matrix(c(0.9, 0.95, 0.99))

  r0 <- 0.01 + 1.8 / sqrt(t)
  swindow0 <- floor(r0 * t)

  dim <- t - swindow0 + 1
  Msadfs <- matrix(0, nrow = m, ncol = dim)

  for(r2 in swindow0:t){
    set.seed(2)
    e <- replicate(m, rnorm(r2))
    a <- r2^(-1)
    z <- e + a
    y <- apply(z, 2, cumsum)

    badfs <- matrix(0, nrow = (r2 - swindow0 + 1), ncol = m)

    for(j in 1:m){
      for(r1 in 1:(r2 - swindow0 + 1)){
        badfs[r1, j] <- as.numeric(ADF_FL(matrix(y[(r1:r2),j]), adflag, mflag))
      }
    }

    if(r2 == swindow0){
      sadfs <- badfs
    } else {
      sadfs <- apply(badfs, 2, max)
    }
    Msadfs[, (r2 - swindow0 + 1)] <- t(sadfs)
  }

  quantile_bsadfs <- matrix(NA, nrow = nrow(qe), ncol = ncol(Msadfs))

  for(i in 1:ncol(Msadfs)){
    quantile_bsadfs[, i] <- as.matrix(quantile(Msadfs[,i], qe, na.rm = T))
  }

  rownames(quantile_bsadfs) <- c('0.9', '0.95', '0.99')
  result <- list(values = Msadfs, quantiles = quantile_bsadfs)
  return(result)
}
