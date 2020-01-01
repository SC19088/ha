#' Compare different orders k in B-spline regression functions with any m equidistant inner knots
#'
#' This function is B-spline regression functions with any m equidistant inner knots and different orders k = 1 to 6.
#'
#' @param m regression functions with m equidistant inner knots
#' @param n the number of data which are simulated by f=sin(12*(x+0.2))/(x+0.2)
#'
#' @return numeric vector
#'
#' @export
#'
#' @examples
#' Bs(100,6)

Bs <- function(n,m){
  f <- function(x) sin(12*(x+0.2))/(x+0.2)
  knots_eq <- function(x, k, m)
  {
    c(min(x) - ((k-1):0) * (max(x)-min(x))/(m+1),
      seq(from=min(x), to=max(x), length.out=m+2)[-c(1,m+2)],
      max(x) + (0:(k-1)) * (max(x)-min(x))/(m+1))
  }
  x=runif(n,0,1)
  e=rnorm(n,0,1)
  y=f(x)+e
  ks <- 1:6
  cols_k <- rainbow(length(ks))
  par(mai=c(0.65,0.6,0.1,0.1), mgp=c(2,1,0))
  plot(x, y, xlab="x", ylab="y", pch=16)
  abline(v=knots_eq(x, 1, m), lty=c(1,rep(2,m),1), col="grey60")
  for (kk in 1:length(ks))
  {
    k <- ks[kk]
    est <- lm(y ~ -1 +
                splineDesign(x=x, knots=knots_eq(x, k, m), ord=k))
    plot(function(x) splineDesign(x=x, knots=knots_eq(x, k, m), ord=k) %*%
           est$coef,
         from=min(x), to=max(x), n=1001, lwd=2, col=cols_k[kk], add=TRUE)
  }
  legend("bottomright", inset=0.02, legend=c(paste("k = ",ks,sep=""), "knots"),
         lwd=2, lty=c(rep(1,length(ks)),2), col=c(cols_k, "grey60"), bg="white")
}






