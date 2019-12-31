## ----eval=FALSE---------------------------------------------------------------
#  function(hurst=0.7, n=100){
#    delta <- 1/n
#    r <- numeric(n+1)
#    r[1] <- 1
#    for(k in 1:n)
#      r[k+1] <- 0.5 * ((k+1)^(2*hurst) - 2*k^(2*hurst) + (k-1)^(2*hurst))
#    r <- c(r, r[seq(length(r)-1, 2)])
#    lambda <- Re((fft(r)) / (2*n))
#    W <- fft(sqrt(lambda) * (rnorm(2*n) + rnorm(2*n)*1i))
#    W <- n^(-hurst) * cumsum(Re(W[1:(n+1)]))
#    X <- ts(W, start=0, deltat=delta)
#    return(X)
#  }

## ----eval=TRUE----------------------------------------------------------------
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
d <- fbm(hurst=0.5, n=500) 
plot(d)

## ----eval=FALSE---------------------------------------------------------------
#  Bs <- function(n,m){
#    f <- function(x) sin(12*(x+0.2))/(x+0.2)
#    knots_eq <- function(x, k, m)
#    {
#      c(min(x) - ((k-1):0) * (max(x)-min(x))/(m+1),
#        seq(from=min(x), to=max(x), length.out=m+2)[-c(1,m+2)],
#        max(x) + (0:(k-1)) * (max(x)-min(x))/(m+1))
#    }
#    x=runif(n,0,1)
#    e=rnorm(n,0,1)
#    y=f(x)+e
#    ks <- 1:6
#    cols_k <- rainbow(length(ks))
#    par(mai=c(0.65,0.6,0.1,0.1), mgp=c(2,1,0))
#    plot(x, y, xlab="x", ylab="y", pch=16)
#    abline(v=knots_eq(x, 1, m), lty=c(1,rep(2,m),1), col="grey60")
#    for (kk in 1:length(ks))
#    {
#      k <- ks[kk]
#      est <- lm(y ~ -1 +
#                  splineDesign(x=x, knots=knots_eq(x, k, m), ord=k))
#      plot(function(x) splineDesign(x=x, knots=knots_eq(x, k, m), ord=k) %*%
#             est$coef,
#           from=min(x), to=max(x), n=1001, lwd=2, col=cols_k[kk], add=TRUE)
#    }
#    legend("bottomright", inset=0.02, legend=c(paste("k = ",ks,sep=""), "knots"),
#           lwd=2, lty=c(rep(1,length(ks)),2), col=c(cols_k, "grey60"), bg="white")
#  }
#  

## ----eval=TRUE----------------------------------------------------------------
library("splines")
library("nlme")
library("MASS")
library("mgcv")
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
Bs(100,6)

## -----------------------------------------------------------------------------
summary(mtcars)

## -----------------------------------------------------------------------------
plot(mtcars$mpg,mtcars$disp)

## ----echo=TRUE----------------------------------------------------------------
### Nf is the inverse function. 
Nf=function (Sigma,x)
{
  FN=((-2*Sigma^2*log(1-x)))^0.5
  FN
}

### pdf is the  Rayleigh density function.  
pdf=function (Sigma,x)
{
  f=x/Sigma^2*exp(-x^2/(2*Sigma^2))
  f
}
n <- 1000
u <- runif(n)
Sigma=0.5 #Here we can change Sigma 
x <- Nf(Sigma,u)
hist(x, prob = TRUE)
y <- seq(0, 2, .005)
lines(y, pdf(Sigma,y))



## ----echo=TRUE----------------------------------------------------------------
n=1000
x1=rnorm(n,0,1)
x2=rnorm(n,3,1)
choose=function(p)
{
  p.new=sample(c(0,1),n,replace = TRUE,prob = c(1-p,p))
  ### The probility that p.new is 1 is p,the probility that p.new is 0 is 1-p.
  X=p.new*x1+(1-p.new)*x2
  X
}

f=function(x,mu)
{
  f=exp(-(x-mu)^2/2)*(2*pi)^(-0.5)
  f
}


hist(choose(0.75))
hist(choose(0.5))
hist(choose(0.85))
hist(choose(0.95))

## -----------------------------------------------------------------------------
n=1e4
x=runif(n,min=0,max=pi/3)
theta_hat1=mean(sin(x))*(pi/3)
theta_hat2=mean(sin(pi/3-x))*(pi/3)
theta_hat=mean(c(theta_hat1,theta_hat2))
print(c(theta_hat1,theta_hat2))
print(theta_hat)

## -----------------------------------------------------------------------------
theta_hat_MC=theta_hat_AhV=0
for (i in 1:1000){
  n=1e4
  x=runif(n,min=0,max=1)
  theta=exp(-x)/(1+x^2)
  theta_hat_MC[i]=mean(theta)
  theta_NI=exp(-(1-x))/(1+(1-x)^2)
  theta_hat_AhV[i]=mean(c(theta[1:(n/2)],theta_NI[1:(n/2)]))
}
c(sd(theta_hat_MC)^2,sd(theta_hat_AhV)^2,sd(theta_hat_AhV)^2/sd(theta_hat_MC)^2)


## -----------------------------------------------------------------------------
MC.Phi <- function(x, R = 10000, antithetic = TRUE) {
  u <- runif(R/2,0,x)
  if (!antithetic) v <- runif(R/2) else v <- 1 - u
  u <- c(u, v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g <- exp(-u)/(1+u^2)
    cdf[i] <- mean(g)
  }
  cdf
}
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1
for (i in 1:m) {
  MC1[i] <- MC.Phi(x, R = 1000, anti = FALSE)
  MC2[i] <- MC.Phi(x, R = 1000)
}

c(sd(MC1)^2,sd(MC2)^2,sd(MC2)^2/sd(MC1)^2)

## -----------------------------------------------------------------------------
M <- 10000; k <- 5
r <- M/k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
est <- 0
g<-function(x) exp(-x - log(1+x^2)) * (x > 0) * (x < 1) / (exp(-x) / (1 - exp(-1)))
for (i in 1:N) {
  for(j in 1:k) 
    {
      u=runif(r,(j-1)/k,j/k)
      x = -log(1 - u * (1 - exp(-1)))
      fg <- g(x)
      T[j]<-mean(fg)
    }
  est[i] <- mean(T)
}
mean(est)
sd(est)

## -----------------------------------------------------------------------------
n <- 20
alpha <- 0.05
UCL <- sigma_hat <- theta_hat <- 0

### repete 1000 times
for (i in 1:1000)
{
  x<- rchisq(n,df=2)
  theta_hat[i]=mean(x)
  sigma_hat[i]=sd(x)
  UCL[i] <- (n-1) * var(x) / qchisq(alpha, df=n-1)
}

### t-interval method
mean(theta_hat-qt(1-alpha/2,n-2)*sigma_hat/sqrt(n)<=2 & theta_hat+qt(1-alpha/2,n-2)*sigma_hat/sqrt(n)>=2)

### method in example 6.4
mean(UCL>=4)


## -----------------------------------------------------------------------------
n <- 20
rep=1000
skewness <- 0
q=c(0.025,0.05,0.95,0.975)

sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

for (i in 1:rep){
  x <- rnorm(n,0,1)
  skewness[i]=sk(x)
}

### the estimated quantiles
cv_hat=as.vector(quantile(skewness,q))
cv_hat
### the real data
cv <- qnorm(q,0, sqrt(6/n))
cv
### the error is: 
1-cv_hat/cv
### Because we get a small number of skewness, we change the variance of the skewness to 6*(n-2) / ((n+1)*(n+3)
cv <- qnorm(q,0, sqrt(6*(n-2)/((n+1)*(n+3))))
### Revised cv is
cv
### the revised error is: 
1-cv_hat/cv
### The result is much better.

### When q=c(0.025,0.05,0.95,0.975),the standard error is
q*(1-q)/(n*(exp(-cv_hat^2)/sqrt(2*pi)))

## -----------------------------------------------------------------------------
alpha <- .1
n <- seq(100,2000,100)
m <- 2500
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
power=0
for (j in 1:length(n)){
  #critical value for the skewness test
  cv <- qnorm(1-alpha/2, 0, sqrt(6*(n[j]-2) / ((n[j]+1)*(n[j]+3))))
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
  x <- rbeta(n[j],alpha,alpha)
  sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
power[j]=mean(sktests)
}
power

## -----------------------------------------------------------------------------
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
e <- epsilon[j]
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
sigma <- sample(c(1, 10), replace = TRUE,
size = n, prob = c(1-e, e))
## x <- rnorm(n, 0, sigma)
x <- rbeta(n,sigma,sigma)
sktests[i] <- as.integer(abs(sk(x)) >= cv)
}
pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------

n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values
for (j in 1:m) {
x <- rchisq(n,df=1)
ttest <- t.test(x, mu = mu0,alternative = c("two.sided"))
p[j] <- ttest$p.value
}
p.hat1 <- mean(p < alpha)


for (j in 1:m) {
x <- runif(n, 0,2)
ttest <- t.test(x,  mu = mu0)
p[j] <- ttest$p.value
}
p.hat2 <- mean(p < alpha)


for (j in 1:m) {
x <- rexp(n, mu0)
ttest <- t.test(x, mu = mu0,alternative = c("two.sided"))
p[j] <- ttest$p.value
}
p.hat3 <- mean(p < alpha)
print(c(p.hat1,p.hat2,p.hat3))


## -----------------------------------------------------------------------------
# get data
library("boot")
data(scor,package = "bootstrap")
# the scatter plots
pairs(scor)
# the sample correlation matrix
scm=cor(scor)
scm
##  Question 1:Compare the plot with the sample correlation matrix.
##  We can find that the larger the correlation coefficient in the sample matrix, the denser the point in the plot.

# For convenience, I used the boot function directly.
r= function(x,i){
 c(cor(x[i,1], x[i,3]),cor(x[i,3], x[i,4]),cor(x[i,3], x[i,5]),cor(x[i,4], x[i,5]))
}

obj <- boot(data = scor, statistic = r, R = 2000)
# The result of bootstrap.
obj
# The standard error of four correlation is 0.08901265,0.04783641,0.06002932,0.06718623.rho12 is the maximal,and the rho 34 is the minimum.


## -----------------------------------------------------------------------------
library("boot")

func.boot <-function(x,i) {
  #computes the sample skewness coeff.
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3)
  m2 <- mean((x[i] - xbar)^2)
  return( m3/m2^1.5 )
}
norm=0
norm.left=norm.right=basic.left=basic.right=perc.left=perc.right=0

# get three  proportions of normal distribution.
m=0
for (i in 1:500){
  nor=rnorm(10)
    obj.nor = boot(data = nor, statistic = func.boot, R = 999)
  CI=boot.ci(boot.out = obj.nor, type = c("norm","basic", "perc"))
  norm[i]<-m>=CI$normal[2]&m<=CI$normal[3]
  norm.left[i]=m<CI$normal[2]
  norm.right[i]=m>CI$normal[3]
  basic.left[i]=m<CI$basic[4]
  basic.right[i]=m>CI$basic[5]
  perc.left[i]=m<CI$perc[4]
  perc.right[i]=m>CI$perc[5]
}
left.norm=c(mean(norm.left),mean(basic.left),mean(perc.left))
right.norm=c(mean(norm.right),mean(basic.right),mean(perc.right))

# get three proportions of Chi-square distribution
m=4/sqrt(10)
for (i in 1:500){
  chi=rchisq(10,df=5,ncp=0)
  obj.nor = boot(data = chi, statistic = func.boot, R = 999)
  CI=boot.ci(boot.out = obj.nor, type = c("norm","basic", "perc"))
  norm[i]=m>=CI$normal[2]&m<=CI$normal[3]
  norm.left[i]=m<CI$normal[2]
  norm.right[i]=m>CI$normal[3]
  basic.left[i]=m<CI$basic[4]
  basic.right[i]=m>CI$basic[5]
  perc.left[i]=m<CI$perc[4]
  perc.right[i]=m>CI$perc[5]
}
left.chi=c(mean(norm.left),mean(basic.left),mean(perc.left))
right.chi=c(mean(norm.right),mean(basic.right),mean(perc.right))

# the answer of normal distribution
data.frame(left.norm,time.norm=1-left.norm-right.norm,right.norm,row.names=c("norm","basic", "perc"))

# the answer of chi-square distribution
data.frame(left.chi,time.chi=1-left.chi-right.chi,right.chi,row.names=c("norm","basic", "perc"))



## -----------------------------------------------------------------------------
# get data
library("boot")
data(scor,package = "bootstrap")
# the correlation matrix
sigma=cov(scor)
# the function to caculate theta
gettheta=function(sigma){
  eigenvalues=eigen(sigma)$values
  theta=max(eigenvalues)/sum(eigenvalues)
  theta
}
# the sample estimate of theta is 
theta.hat=gettheta(sigma)
print(theta.hat)
# do the jackknife
n=length(scor[,1])
theta <- scor.jack <- 0
for (i in 1:n)
{
  scor.jack=scor[-i,]
  theta[i]=gettheta(cov(scor.jack))
}
# the bias is
bias=(n-1)*(mean(theta)-theta.hat)
print(round(bias,6))
# the standard error
se <- sqrt((n-1) * mean((theta.hat-mean(theta))^2))
se <- (n-1)*sqrt(var(theta)/n)
print(round(se,4))

## -----------------------------------------------------------------------------
library(lattice)
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
r1 <- r2 <- r3 <- r4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1
r1[k] <- summary(J1)$adj.r.squared
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
r2[k] <- summary(J2)$adj.r.squared
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
r3[k] <- summary(J3)$adj.r.squared
J4 <- lm(y ~ x + I(x^2)+I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
e4[k] <- magnetic[k] - yhat4
r4[k] <- summary(J4)$adj.r.squared
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
# From the result,we can find the $e_{2}^{2}$ is the minimum.According to the prediction error criterion, Model 2, the quadratic modelwould be the best fit for the data.
c(mean(r1), mean(r2), mean(r3), mean(r4))
# From the result,we can find the$R_{1}^{2}$ is the maximum.According to maximum adjusted, the quadratic modelwould be the best fit for the data.

## -----------------------------------------------------------------------------

f <- function(x)   0.5*exp(-abs(x))  

rw.Metropolis <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (f(y) / f(x[i-1])))
x[i] <- y else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}
n <- 4 #degrees of freedom for target Student t dist.
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 10
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
# print the reject rate
print(c(rw1$k, rw2$k, rw3$k, rw4$k)/2000)
no.reject <- data.frame(sigma=sigma,no.reject=c(rw1$k, rw2$k, rw3$k, rw4$k))
knitr::kable(no.reject,format='latex')
# plot
   #display 4 graphs together
    refline <- qt(c(.025, .975), df=n)
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }
     #reset to default






## -----------------------------------------------------------------------------
# the function f is use"==" to judge if log(exp x) = exp(logx) = x
f=function(x){
  a=log(exp(x))==exp(log(x))
  b=log(exp(x))==x
  c=exp(log(x))==x
  a+b+c==3
}
f(1:10)
log(exp(10))-exp(log(10))
log(exp(10))-10
exp(log(10))-10
# We can find the example show the caculation about exp(log(x)) may be wrong.

# the function g is use all.equal to judge if log(exp x) = exp(logx) = x
g=function(x){
  all.equal(log(exp(x)),exp(log(x)),x)
}
g(1:10)
# There are no errors.



## -----------------------------------------------------------------------------
root.114 <- root.115 <- 0
k<-c(4:25,100,500,1000)

# caculate the root in 11.4
f1=function(a,k){
  pt(sqrt(a^2*k/(k+1-a^2)),df=k)-pt(sqrt(a^2*(k-1)/(k-a^2)),df=k-1)
}
for (i in 1:25){
  root.114[i]=uniroot(f1,interval = c(-1,sqrt(k[i])-0.01),k=k[i])$root
}



# caculate the root in 11.5
ck=function(a,k){
  sqrt(a^2*k/(k+1-a^2))
}

# equation function in 11.5
f2=function(a,k){
  log(2)+lgamma(k/2)-0.5*log(pi*(k-1))-lgamma((k-1)/2)+log(integrate(function(u){(1+u^2/(k-1))^(-k/2)},0,ck(a,k-1))$value)-(log(2)+lgamma((k+1)/2)-0.5*log(pi*k)-lgamma(k/2)+log(integrate(function(u){(1+u^2/k)^(-(k+1)/2)},0,ck(a,k))$value))
}

for (i in 1:22){
  root.115[i]=uniroot(f2,interval = c(-1,root.114[i]+0.0001),k=k[i])$root
}

## Note: When k is large, the value range of the uniroot function must be fine, otherwise the two ends of the interval are likely to have the same number. So I set the upper bound separately for k = 100, 500, 1000 
root.115[23]=uniroot(f2,interval = c(-1,9.98),k=k[23])$root
root.115[24]=uniroot(f2,interval = c(-1,sqrt(k[24])-0.03),k=k[24])$root
root.115[25]=uniroot(f2,interval = c(-1,sqrt(k[25]-1)),k=k[25])$root

root.114 #the result in 11.4

root.115 #the result in 11.5




## -----------------------------------------------------------------------------

set.seed(12345)
p <- q <- 0.1
r=1-p-q
lf <- 0

LF=function(AA,BB,p,q){ #log-likelyhood function
  a=(28-AA)*log(2*p*r)+AA*log(p^2)+(24-BB)*log(2*q*r)+BB*log(q^2)+70*log(2*p*q)+2*41*log(r)
  return(a)
}

iter<-function(p,q){
  r<-1-p-q
  a<-2*28*(p^2)/(p^2+2*p*r)+2*p*r/(p^2+2*p*r)*28+70
  b<-2*24*(q^2)/(q^2+2*q*r)+2*q*r/(q^2+2*q*r)*24+70
  c<-2*41+2*p*r/(p^2+2*p*r)*28+2*q*r/(q^2+2*q*r)*24
  return(c(a/(a+b+c),b/(a+b+c)))
}

AA=sum(rbinom(28,1,(p^2/(p^2+2*p*r))))
BB=sum(rbinom(24,1,(q^2/(q^2+2*q*r))))
lf=LF(AA,BB,p,q)

for (i in 1:100){
  m=iter(p[i],q[i])
  p[i+1]=m[1];q[i+1]=m[2]
  AA=sum(rbinom(28,1,(p[i+1]^2/(p[i+1]^2+2*p[i+1]*r))))
  BB=sum(rbinom(24,1,(q[i+1]^2/(q[i+1]^2+2*q[i+1]*r))))
  lf[i+1]=LF(AA,BB,p[i+1],q[i+1])
  if(abs(p[i+1]-p[i])<=1e-4) break
}

p

q

plot(lf,type="o",xlab = "n",ylab = "LF",main="The plot of log-likelyhood function")


## -----------------------------------------------------------------------------
formulas =list(mpg ~ disp,mpg ~ I(1 / disp),mpg ~ disp + wt,mpg ~ I(1 / disp) + wt)
# lapply() to fit linear models
lapply(formulas,lm,data=mtcars)
a<-list(NULL)
length(a)<-4
# for loop to fit linear models
for (i in 1:4){
  a[[i]]=lm(formula = formulas[[i]],data=mtcars)
}
a

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
lapply(bootstraps,lm,formula=mpg~disp)


## -----------------------------------------------------------------------------
# the function
rsq <- function(mod) summary(mod)$r.squared
formulas =list(mpg ~ disp,mpg ~ I(1 / disp),mpg ~ disp + wt,mpg ~ I(1 / disp) + wt)
# lapply() to fit linear models
a=lapply(formulas,lm,data=mtcars)
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
b=lapply(bootstraps,lm,formula=mpg~disp)
# r^2 in Question 1
lapply(a, rsq)
# r^2 in Question 2
lapply(b, rsq)


## -----------------------------------------------------------------------------

trials <- replicate(100,t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
sapply(trials, '[[', 3)


## -----------------------------------------------------------------------------
# We can use the function parSapply() to do this work.
library(parallel) 
cores <- detectCores()
cl <- makePSOCKcluster(cores-1)
#parSapply(cl = NULL, X, FUN, ..., simplify = TRUE,          USE.NAMES = TRUE, chunk.size = NULL)
parSapply(cl,1:10,FUN=sin)



## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
set.seed(12345)

# Use Rcpp to solve the problem
cppFunction('List Metropolis(double sigma, double x0 , int N){
  NumericVector x(N);
  NumericVector u=runif(N);
	x[0]=x0;
	int k=0;
  for(int i = 1; i < N; ++i) {
    	double y=rnorm(1,x[i-1],sigma)[0];
  	  double p=exp(-abs(y))/exp(-abs(x[i-1]));
	    if(u[i-1]<= p) {x[i]=y;}
      else{x[i]=x[i-1]; k=k+1; } 
  	}
  List L=List::create(Named("x")=x,
                      Named("k")=k);
  return L;
}')

N = 2000
sigma = c(.05, .5, 2, 16)
x0 = 25

cpp1=Metropolis(sigma[1],x0,N)
cpp2=Metropolis(sigma[2],x0,N)
cpp3=Metropolis(sigma[3],x0,N)
cpp4=Metropolis(sigma[4],x0,N)

#number of candidate points rejected
Rej = cbind(cpp1$k, cpp2$k, cpp3$k, cpp4$k)
Acc = round((N-Rej)/N,4)
rownames(Acc) = "Accept rates"
colnames(Acc) = paste("sigma",sigma)
knitr::kable(Acc)
#plot
 #display 4 graphs together
    rw = cbind(cpp1$x, cpp2$x, cpp3$x,  cpp4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
    }





# Use R to solve the problem
lap_f = function(x) exp(-abs(x))
rw.Metropolis = function(sigma, x0, N){
 x = numeric(N)
 x[1] = x0
 u = runif(N)
 k = 0
 for (i in 2:N) {
  y = rnorm(1, x[i-1], sigma)
  if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
  else {
  x[i] = x[i-1]
  k = k+1
  }
 }
 return(list(x = x, k = k))
}

rw1 = rw.Metropolis(sigma[1],x0,N)
rw2 = rw.Metropolis(sigma[2],x0,N)
rw3 = rw.Metropolis(sigma[3],x0,N)
rw4 = rw.Metropolis(sigma[4],x0,N)

#number of candidate points rejected
Rej = cbind(rw1$k, rw2$k, rw3$k, rw4$k)
Acc = round((N-Rej)/N,4)
rownames(Acc) = "Accept rates"
colnames(Acc) = paste("sigma",sigma)
knitr::kable(Acc)
#plot
  #display 4 graphs together
    rw = cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
    }
    

## Compare the generated random numbers by the two functions using qqplot.
## Acturally, I don't know what this problem means.I compared x in two methods,but it is doesn't make sense.
    

qqplot(cpp1$x,rw1$x)    
qqplot(cpp2$x,rw2$x)   
qqplot(cpp3$x,rw3$x)
qqplot(cpp4$x,rw4$x)
    
    
    
    
## Campare the computation time of the two functions with microbenchmark.
ts1=microbenchmark(rpp=Metropolis(sigma[1],x0,N),r=rw.Metropolis(sigma[1],x0,N))
summary(ts1)[,c(1,3,5,6)]
ts2=microbenchmark(rpp=Metropolis(sigma[2],x0,N),r=rw.Metropolis(sigma[2],x0,N))
summary(ts2)[,c(1,3,5,6)]
ts3=microbenchmark(rpp=Metropolis(sigma[3],x0,N),r=rw.Metropolis(sigma[3],x0,N))
summary(ts3)[,c(1,3,5,6)]
ts4=microbenchmark(rpp=Metropolis(sigma[4],x0,N),r=rw.Metropolis(sigma[4],x0,N))
summary(ts4)[,c(1,3,5,6)]




