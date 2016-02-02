# ---------------------------------------------------------------------------------------
# This function generates and fit the Wald distribution
# using GAMLSS framewok as function of
#
# mu > 0 = mu (rate of evidence accrual)
# a  > 0 = sigma (response threshold)
# We use the parameterization given by Heathcote, A. (2004). Fitting Wald and
# ex-Wald distributions to response time data: An example using functions for
# the S-PLUS package. Behavior Research Methods, Instruments, &
# Computers,36,678-694.
#
# Create by Raydonal Ospina, 2016
# Modified by:
# Raydonal  2016
# ---------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
WALD <-function (mu.link = "log", sigma.link = "log")
{
    mstats <- checklink("mu.link", "Wald", substitute(mu.link), c("1/mu^2", "inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Wald", substitute(sigma.link),  c("inverse", "log", "identity", "own"))

    structure(
          list(family = c("WALD", "Wald distribution"),
           parameters = list(mu=TRUE, sigma=TRUE),
                nopar = 2,
                 type = "Continuous",

              mu.link = as.character(substitute(mu.link)),
           sigma.link = as.character(substitute(sigma.link)),

           mu.linkfun = mstats$linkfun,
        sigma.linkfun = dstats$linkfun,

           mu.linkinv = mstats$linkinv,
        sigma.linkinv = dstats$linkinv,

                mu.dr = mstats$mu.eta,
             sigma.dr = dstats$mu.eta,

              #first derivate of log-density respect to mu
                dldm = function(y, mu, sigma) {
                                dldm = sigma-(y*mu)
                                dldm
                                }, #OK

              #second derivate of log-density respect to mu
              d2ldm2 = function(y, mu,sigma)  {
                                dldm2 = -y #sigma/mu
                                dldm2
                                }, # -y, #OK

              #first derivate of log-density respect to sigma
                dldd = function(y,mu,sigma) {
                                dldd = mu + (1/sigma) - (sigma/y)
                                dldd
                                }, # OK

              #second derivate of log-density respect to sigma
              d2ldd2 = function(y, mu, sigma) {
                                d2ldd2 = (-1/sigma^2) - (1/y) #1 #(2/sigma^2) +  (mu/sigma)
                                d2ldd2
                                },   # (), # OK,

              # partial derivate of log-density respect to mu and sigma
             d2ldmdd = function(y, mu, sigma){
                                d2ldmdd =  rep(-1,length(y))
                                d2ldmdd
             },

         G.dev.incr  = function(y, mu, sigma,...) -2*dWALD(y,mu,sigma,log=TRUE),

                rqres = expression(rqres(pfun="pWALD", type="Continuous", y=y, mu=mu, sigma=sigma)),

#            mu.initial = expression( mu <- (y+mean(y))/2) ,
#         sigma.initial = expression(sigma <- rep(1, length(y))), #sd(y)/(mean(y))^1.5 ),

#         # Initial values based on moments
        mu.initial = expression( mu <- rep(sqrt(mean(y)/var(y)), length(y))),
        sigma.initial = expression(sigma <- rep(sqrt(mean(y)/var(y))*mean(y), length(y))),

              mu.valid = function(mu) all(mu > 0),
           sigma.valid = function(sigma)  all(sigma > 0),
               y.valid = function(y)  all(y > 0)
            ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
dWALD<-function(x, mu = 1, sigma = 1, log=FALSE) # OK
 {        if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))
 log.lik <- (-0.5*log(2*pi))+log(sigma)-((3/2)*log(x)) - ((sigma - (mu*x))^2 / (2*x))
           if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
      fy
  }
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
pWALD <- function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) # OK
  {    #  browser()
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
    if (any(q < 0))  stop(paste("y must be positive", "\n", ""))
      lq <- length(q)
   sigma <- rep(sigma, length = lq)
      mu <- rep(mu, length = lq)

    cdf1 <- pnorm(((mu*q) - sigma)/sqrt(q))

   lcdf2 <- ( 2*mu*sigma)+pnorm(- ( ((mu*q) + sigma)/sqrt(q)),log.p=TRUE)
        cdf <- cdf1+ exp(lcdf2)

    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf)
    cdf
   }
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
qWALD <- function(p, mu=1, sigma=1,  lower.tail = TRUE, log.p = FALSE)
 {
    #---functions--------------------------------------------
       h1 <- function(q)
       {
     pWALD(q , mu = mu[i], sigma = sigma[i])-p[i]
       }
       h <- function(q)
       {
     pWALD(q , mu = mu[i], sigma = sigma[i])
       }
     #-------------------------------------------------------
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", ""))
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))
        lp <-  max(length(p),length(mu),length(sigma))
          p <- rep(p, length = lp)
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
          q <- rep(0,lp)
         for (i in  seq(along=p))
         {
         if (h(mu[i])<p[i])
          {
           interval <- c(mu[i], mu[i]+sigma[i])
           j <-2
           while (h(interval[2]) < p[i])
              {interval[2]<- mu[i]+j*sigma[i]
              j<-j+1
              }
           }
          else
           {
           interval <-  interval <- c(.Machine$double.xmin, mu[i])
           }
        q[i] <- uniroot(h1, interval)$root
         }
    q
   }
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
rWALD <- function(n, mu=1, sigma=1, ...) # OK
  {
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))
    n <- ceiling(n)
    p <- runif(n)
    r <- qWALD(p,mu=mu,sigma=sigma, ...)
    r
  }
#----------------------------------------------------------------------------------------

