#' Generate a time series of fractional Brownian motion.
#'
#' This function generatea a time series of one dimension fractional Brownian motion.
#'
#' @param hurst the hurst index, with the default value 0.71
#' @param n the number of points between 0 and 1 that will be generated, with the default value 100
#'
#' @return numeric vector
#'
#' @export
#'
#' @examples
#' d <- fbm(hurst=0.5, n=500)
#' plot(d)
fbm <- function(hurst=0.7, n=100){
  delta <- 1/n
  r <- numeric(n+1)
  r[1] <- 1
  for(k in 1:n)
    r[k+1] <- 0.5 * ((k+1)^(2*hurst) - 2*k^(2*hurst) + (k-1)^(2*hurst))
  r <- c(r, r[seq(length(r)-1, 2)])
  lambda <- Re((fft(r)) / (2*n))
  W <- fft(sqrt(lambda) * (rnorm(2*n) + rnorm(2*n)*1i))
  W <- n^(-hurst) * cumsum(Re(W[1:(n+1)]))
  X <- ts(W, start=0, deltat=delta)
  return(X)
}
