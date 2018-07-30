#' @title Augmented Dickey-Fuller Statistic by AIC or BIC
#'
#' @description Calculate the Augmented Dickey-Fuller Statistic with lag order selected by AIC or BIC.
#'
#' @param y the time series to be used.
#' @param adflag the maximum lag order.
#' @param mflag 1 for ADF with constant and whithout trend, 2 for ADF with constant and trend and 3 for ADF without constant and trend.
#' @param IC 1 for AIC and 2 for BIC.
#' @keywords AugmentedDickey-FullerTest.
#' @references Phillips, P.C. & Shi, S. & Yu, J. (2013). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500". \emph{SSRN Electronic Journal}.
#' @export
#' @examples
#' y <- rnorm(10)
#' ADF_IC(y, adflag = 1, mflag = 2, IC = 1)
#' ADF_IC(y, adflag = 1, mflag = 2, IC = 2)
#'

ADF_IC <- function(y, adflag, mflag, IC){
  t0 <- length(y)
  t1 <- length(y) - 1
  const <- matrix(1, nrow = t1, ncol = 1)
  trend <- as.matrix(1:t1)

  y1 <- y[1:t1]
  dy <- as.matrix(y[2:t0] - y[1:t1])
  x <- as.matrix(y1)

  if (mflag == 1){
    x <- cbind(x, const)
  }
  if (mflag == 2){
    x <- cbind(x, const, trend)
  }
  if(mflag == 3){
    x <- x
  }
  x1 <- as.matrix(x)

  ######### whenk k=0 ############
  beta0 <- solve(t(x1) %*% x1, t(x1)) %*% dy
  eps0 <- dy - x1 %*% beta0
  T0 <- length(eps0)

  ### IC0=-2*sum(log(pdf('norm',eps0)))/T0+length(beta0)*log(T0)/T0;
  ## this command greatly slow down the program

  npdf0 <- sum(-1/2 * log(2 * pi) - 1/2 * (eps0^2))

  if(IC == 1){
    IC0 <- -2 * npdf0 / T0 + 2 * length(beta0) / T0
  } else if(IC == 2){
    IC0 <- - 2 * npdf0 / T0 + length(beta0) * log(T0) / T0
  }

  if (mflag == 1){
    se0 <- t(eps0) %*% eps0 / (t1 - 2)
  }
  if (mflag == 2){
    se0 <- t(eps0) %*% eps0 / (t1 - 3)
  }
  if (mflag == 3){
    se0 <- t(eps0) %*% eps0 / (t1 - 1)
  }

  sig0 <- as.matrix(sqrt(diag(as.numeric(se0) * solve(t(x1) %*% x1))))
  ADF0 <- beta0[1] / sig0[1]

  ####### when k>0 #######
  if(adflag > 0){
    IC1 <- matrix(0, nrow = adflag, ncol = 1)
    ADF1 <- matrix(0, nrow = adflag, ncol = 1)

    for(k in 1:adflag){
      t2 <- t1 - k
      xx <- x1[((k + 1):t1),]
      dy01 <- matrix(dy[(k + 1):t1])

      x2 <- cbind(xx, matrix(0, nrow = t2, ncol = k))

      for(j in 1:k){
        x2[, (ncol(xx) + j)] <- dy[(k + 1 - j):(t1 - j)]
      }

      beta <- solve((t(x2) %*% x2), (t(x2) %*% dy01))
      eps <- dy01 - x2 %*% beta

      t <- length(eps)

      npdf <- sum(-1/2 * log(2 * pi) - 1/2 * (eps^2))
      if(IC == 1){
        IC1[k] <- - 2 * npdf / t + 2 * nrow(beta) / t
      } else if(IC == 2){
        IC1[k] <- - 2 * npdf / t + nrow(beta) * log(t) / t
      }

      if(mflag == 1){
        se <- t(eps) %*% eps / (t2 - adflag - 2)
      }
      if(mflag == 2){
        se <- t(eps) %*% eps / (t2 - adflag - 3)
      }
      if(mflag == 3){
        se <- t(eps) %*% eps / (t2 - adflag - 1)
      }
      sig <- sqrt(diag(as.numeric(se) * solve(t(x2) %*% x2)))
      sig <- matrix(sig)
      ADF1[k] <- beta[1] / sig[1]
    }
  }

  ICC <- as.matrix(c(IC0, IC1))
  ADF <- as.matrix(c(ADF0, ADF1))

  lag <- which(ICC == min(ICC))[1]

  ADFlag <- ADF[lag]

  if(IC == 1){
    result <- list('ADF Statistic using AIC' = ADFlag)
  }
  if(IC == 2){
    result <- list('ADF Statistic using BIC' = ADFlag)
  }

    return(result)
}

