#' @title Critical values for sup ADF statistic sequence.
#'
#' @description Calculate critical value sequences for the sup ADF statistic sequence using
#'  Monte Carlo simulations for a sample generated from a Normal distribution.
#'
#' @param m Number of Monte Carlo Simulations. Default equals 2000. Must be bigger than 2.
#' @param t Sample size. Default equals 100. Must be bigger than 2.
#' @keywords AugmentedDickey-FullerTest, supADFSequence MonteCarlo.
#' @references Phillips, P.C. & Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500". \emph{SSRN Electronic Journal}.
#' @import stats
#' @import foreach
#' @import MASS
#' @export
#'
#' @examples
#' foo <- sadf(m = 20, t = 50)
#' quant <- rep(foo$quantiles[2], 100)
#' plot(quant, type = 'l')

sadf <- function(m, t){
  qe <- matrix(c(0.9, 0.95, 0.99))
  r0 <- 0.01 + 1.8 / sqrt(t)
  swindow0 <- floor(r0 * t)

  y <- DGP(t, m)

  badfs <- matrix(0, nrow = (t-swindow0+1), ncol = m)
  sadf <- matrix(NA, nrow = m, ncol = 1)

  for(j in 1:m){
    for(i in swindow0:t){
      badfs[i - swindow0 + 1,j] <- as.numeric(ADF_FL(y[1:i, j], 0, 1))
    }
  }

  sadf[,1] <- apply(badfs, 2, max)
  quantile_sadf <- quantile(sadf, qe)

  result <- list(value = sadf, quantiles = quantile_sadf)
  return(result)
}
