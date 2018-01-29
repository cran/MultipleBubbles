#' @title Backward Augmented Dickey-Fuller Sequence.
#'
#' @description In this program, we calculate critical value sequences for the backward ADF statistic
#' sequence for a matrix generated from a standard Normal distribution.
#'
#' @param m Number of Monte Carlo replications. Must be bigger than 2.
#' @param t Sample size. Must be bigger than 2.
#' @param adflag Number of lags to be included in the ADF Test. Default equals 0.
#' @param mflag 1 for ADF with constant and whithout trend, 2 for ADF with constant and trend and 3 for ADF without constant and trend.
#' @keywords AugmentedDickey-FullerTest BackwardADFSequence MonteCarlo.
#' @references Phillips, P.C. & Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500". \emph{SSRN Electronic Journal}.
#' @import stats
#' @import foreach
#' @import MASS
#' @export
#'
#' @examples
#' foo <- badf(m = 100, t = 50, adflag = 1, mflag = 1)
#' plot(foo$quantiles[2,], type = 'l')

badf <- function(m, t, adflag = 0, mflag = 1){
  qe <- as.matrix(c(0.9, 0.95, 0.99))

  r0 <- 0.01 + 1.8 / sqrt(t);
  swindow0 <- floor(r0 * t)

  dim <- t - swindow0 + 1
  adfs <- matrix(data = 0, nrow = m, ncol = dim)

  for (r2 in swindow0:t){
    set.seed(1)
    e <- replicate(m, rnorm(r2))
    a <- r2^(-1)
    z <- e + a
    y <- apply(z, 2, cumsum)
    for(j in 1:m){
      adfs[j, (r2 - swindow0 + 1)] <- as.numeric(ADF_FL(y[,j], adflag, mflag))
    }
  }

  quantile_badfs <- matrix(NA, nrow = nrow(qe), ncol = ncol(adfs))

  for(i in 1:ncol(adfs)){
    quantile_badfs[, i] <- as.matrix(quantile(adfs[,i], qe, na.rm = T))
  }
  rownames(quantile_badfs) <- c('0.9', '0.95', '0.99')

  result <- list(values = adfs, quantiles = quantile_badfs)
  return(result)
}
