#' @title Sup ADF and generalized sup ADF statistics for a time series.
#'
#' @description Calculate the sup ADF and the generalized sup ADF statistics using the
#'  backward ADF statistic sequence and the backward SADF statistic sequence, respectively.
#'
#' @param y the time series.
#' @param adflag the lag order for the ADF test.
#' @param mflag 1 for ADF with constant and whithout trend, 2 for ADF with constant and trend
#'  and 3 for ADF without constant and trend.
#' @param parallel If TRUE, uses parallel computing for the loop. If the data is large it could be faster,
#'  but usually it is slower for small data.
#' @param IC 1 for AIC and 2 for BIC.
#' @references Phillips, P.C. & Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500". \emph{SSRN Electronic Journal}.
#' @import stats
#' @import foreach
#' @import MASS
#' @export
#'
#' @examples{
#' y <- rnorm(100)
#' foo <- sadf_gsadf(y, adflag = 1, mflag = 1, IC = 2)
#' plot(foo$badfs, type = 'l')
#' plot(foo$bsadfs, type = 'l')
#' /dontrun{
#' data(sp_data)
#' y <- sp_data[1:500]
#' foo <- sadf_gsadf(y, adflag = 1, mflag = 1, IC = 1, parallel = FALSE)
#' as <- badf(m = 200, t = length(y))
#' plot(foo$badfs, type = 'l')
#' lines(as$quantiles[2,], col = 'red')
#' fas <- bsadf(m = 200, t = length(y))
#' plot(foo$bsadfs, type = 'l')
#' lines(fas$quantiles[2,], col = 'red')
#' }

sadf_gsadf <- function(y, adflag, mflag, IC, parallel = FALSE){
  t <- length(y)
  r0 <- 0.01 + 1.8 / sqrt(t)
  swindow0 <- floor(r0*t)
  dim <- t - swindow0 + 1

  badfs <- matrix(0, nrow = (t - swindow0 + 1), ncol = 1)

  for(i in swindow0:t){
    badfs[(i - swindow0 + 1), 1] <- as.numeric(ADF_IC(y[(1:i)], adflag, mflag, IC))
  }
  sadf <- max(badfs)

  r2 <- swindow0:t
  r2 <- t(r2)
  rw <- r2 - swindow0 + 1

  bsadfs <- matrix(0, nrow = 1, ncol = dim)

  for(v in 1:length(r2)){
    swindow <- swindow0:r2[v]
    r1 <- r2[v] - swindow + 1
    rwadft <- matrix(0, nrow = length(swindow), ncol = 1)

    for(i in 1:length(swindow)){
      rwadft[i] <- as.numeric(ADF_IC(y[(r1[i]:r2[v])], adflag, mflag, IC))
    }
    bsadfs[1, v] <- max(rwadft)
  }
  if(parallel == T){
    foreach(v = 1:length(r2)) %dopar% {
      swindow <- swindow0:r2[v]
      r1 <- r2[v] - swindow + 1
      rwadft <- matrix(0, nrow = length(swindow), ncol = 1)
      foreach(i = 1:length(swindow)) %dopar% {
        rwadft[i] <- as.numeric(ADF_IC(y[(r1[i]:r2[v])], adflag, mflag, IC))
    }
    bsadfs[1, v] <- max(rwadft)
  }
}
  gsadf <- max(bsadfs[1,])
  bsadfs <- matrix(bsadfs)

  result <- list(badfs = badfs, bsadfs = bsadfs, sadf = sadf, gsadf = gsadf)
  return(result)
}

