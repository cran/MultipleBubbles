#' @title Critical values for generalized sup ADF statistic sequence.
#'
#' @description Calculate critical value sequences for the generalized sup ADF statistic sequence
#'  using Monte Carlo simulations for a sample generated from a Normal distribution.
#'
#' @param m Number of Monte Carlo Simulations. Default equals 2000. Must be bigger than 2.
#' @param t Sample size. Default equals 100. Must be bigger than 2.
#' @param adflag Number of lags to be included in the ADF Test. Default equals 0.
#' @param mflag 1 for ADF with constant and whithout trend, 2 for ADF with constant and trend and 3 for ADF without constant and trend.
#' @param swindow0 Minimum window size.
#' @keywords AugmentedDickey-FullerTest GSADFSequence MonteCarlo.
#' @references Phillips, P.C. & Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500". \emph{SSRN Electronic Journal}.
#' @import stats
#' @import foreach
#' @import MASS
#' @export
#'
#' @examples
#' foo <- gsadf(m = 20, t = 50)
#' quant <- rep(foo$quantiles[2], 100)
#' plot(quant, type = 'l')

gsadf <- function(m, t, adflag = 0, mflag = 1, swindow0 = floor(r0*t)){
  qe <- matrix(c(0.9, 0.95, 0.99))
  r0 <- 0.01 + 1.8 / sqrt(t)
  #swindow0 <- floor(r0 * t)
  dim <- t - swindow0 + 1

  y <- DGP(t, m)

  gsadf <- matrix(1, nrow = m, ncol = 1)

  for(j in 1:m){
    sadfs <- matrix(0, nrow = dim, ncol = 1)
    for(r2 in swindow0:t){
      dim0 <- r2 - swindow0 + 1
      rwadft <- matrix(0, nrow = dim0, ncol = 1)

      for(r1 in 1:dim0){
        rwadft[r1] <- as.numeric(ADF_FL(y[(r1:r2), j], adflag, mflag))
      }
      sadfs[r2 - swindow0 + 1] <- max(rwadft)
    }
    gsadf[j] <- max(sadfs)
  }

  quantile_gsadf <- quantile(gsadf, qe)
  result <- list(values = gsadf, quantiles = quantile_gsadf)

  return(result)
}

