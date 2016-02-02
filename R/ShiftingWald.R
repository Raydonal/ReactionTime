# ---------------------------------------------------------------------------------------
# This function generates and fit the Shifted Wald distribution
# using GAMLSS framewok as function of
#
# m > 0 = mu (rate of evidence accrual)
# a > 0 = sigma (response threshold)
# s > 0 = nu (the shift)
# We use the parameterization given by Heathcote, A. (2004). Fitting Wald and 
# ex-Wald distributions to response time data: An example using functions for 
# the S-PLUS package. Behavior Research Methods, Instruments, & 
# Computers,36,678-694.
#
# Create by Raydonal Ospina, Oct 11, 2014 
# Modified by:
# Raydonal  11/10/2014
# ---------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------  
shiftWALD <-function (mu.link = "log", sigma.link = "log", nu.link ="log") 
{   
    mstats <- checklink("mu.link", "Shift-Wald", substitute(mu.link), c("1/mu^2", "inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Shift-Wald", substitute(sigma.link),  c("1/mu^2", "inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "Shift-Wald", substitute(nu.link), c("1/mu^2", "inverse", "log", "identity", "own"))
    
    structure(
          list(family = c("shiftWALD", "Shift-Wald"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
           
           nopar = 3, 
           
           type = "Continuous",
           
           mu.link = as.character(substitute(mu.link)),  
           
           sigma.link = as.character(substitute(sigma.link)), 
           
           nu.link = as.character(substitute(nu.link)), 
           
           mu.linkfun = mstats$linkfun, 
           
           sigma.linkfun = dstats$linkfun, 
           
           nu.linkfun = vstats$linkfun,
           
           mu.linkinv = mstats$linkinv, 
           
           sigma.linkinv = dstats$linkinv,
           
           nu.linkinv = vstats$linkinv,
           
           mu.dr = mstats$mu.eta, 
           
           sigma.dr = dstats$mu.eta, 
           
           nu.dr = vstats$mu.eta,
        
              #first derivate of log-density respect to mu
                dldm = function(y, mu, sigma, nu) {
                                dldm = sigma-mu*(y-nu)
                                dldm
                                }, #OK
        
              #second derivate of log-density respect to mu
              d2ldm2 = function(y,nu) {
                                dldm2 = nu-y #sigma/mu
                                dldm2
                                }, # -y, #OK
        
              #first derivate of log-density respect to sigma
                dldd = function(y, mu, sigma, nu) {
                                dldd = mu + (1/sigma) - (sigma/(y-nu))                                   
                                dldd
                                }, # OK
        
              #second derivate of log-density respect to sigma
              d2ldd2 = function(y, sigma, nu) {
                                d2ldd2 =   (-1/sigma^2) - (1/(y-nu))
                                d2ldd2
                                },   # (-1/sigma^2 - 1/y), # OK,
        
           
              #first derivate log-density respect to nu
              dldv = function(y, mu, sigma, nu) {       
                              dldv = 3/(2*(y-nu)) - (sigma^2/ (2*(y-nu)^2)) + (mu^2/2) 
                              dldv
                              },
           
              #second derivate log-density respect to nu
              d2ldv2 = function(y, sigma, nu) {         
              d2ldv2 = (3/2 - (sigma^2/(y-nu)))/ (y-nu)^2
              d2ldv2
              },
           
           
              #partial derivate of log-density respect to mu and sigma  
              d2ldmdd = function(y,mu,sigma) {   
              d2ldmdd = rep(1,length(y))
              d2ldmdd
              },
           
              #partial derivate of log-density respect to mu and nu
              d2ldmdv = function(mu) {  
              d2ldmdv = mu
              d2ldmdv
              },
           
              #partial derivate of log-density respect to sigma and nu
              d2ldddv = function(y, sigma, nu) {   
              d2ldddv = -sigma/(y-nu)^2 
              d2ldddv
              },
           
         G.dev.incr  = function(y, mu, sigma, nu, ...) -2*dshiftWALD(y, mu, sigma, nu, log=TRUE),
        
          rqres = expression(rqres(pfun="pshiftWALD", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
        
#         # Initial values based on moments
        mu.initial = expression( mu <- rep(sqrt(mean(y)/var(y)), length(y))),
        sigma.initial = expression(sigma <- rep(sqrt(mean(y)/var(y))*mean(y), length(y))), 
        nu.initial = expression(nu <- rep(min(y), length(y))), 

              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0), 
              nu.valid = function(nu) all(nu > 0), 
               y.valid = function(y)  all(y > 0)
            ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------  

#----------------------------------------------------------------------------------------  
# Shifted Wald density
dshiftWALD <-function(x, mu = 1, sigma = 1, nu=0, log=FALSE) # OK
{        if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(nu < 0))  stop(paste("nu must be positive", "\n", "")) 
         if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
         x = x-nu
         log.lik <- log(sigma) -(sigma-mu*x)^2/(2*x) -0.5*log(2*pi) -((3/2)*log(x))
         if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
         fy 
}



#----------------------------------------------------------------------------------------  
# Shifted Wald cumulative density 
pshiftWALD <- function(q, mu = 1, sigma = 1, nu=0, lower.tail = TRUE, log.p = FALSE) # OK
  {    #  browser() 
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
    if (any(nu < 0))  stop(paste("nu must be positive", "\n", "")) 
    if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
      lq <- length(q)                                                                    
   sigma <- rep(sigma, length = lq)
      mu <- rep(mu, length = lq)           

    q = q-nu
    cdf1 <- pnorm(((mu*q) - sigma)/sqrt(q)) 
   
   lcdf2 <- ( 2*mu*sigma)+pnorm(- ( ((mu*q) + sigma)/sqrt(q)),log.p=TRUE)
        cdf <- cdf1+ exp(lcdf2)
   
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
   }
#----------------------------------------------------------------------------------------  



#----------------------------------------------------------------------------------------  
# Shifted Wald Quantile function
qshiftWALD <- function(p, mu=1, sigma=1, nu=0,  lower.tail = TRUE, log.p = FALSE) 
 {
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pshiftWALD(q , mu = mu[i], sigma = sigma[i], nu=nu[i])-p[i]   
       }
       h <- function(q)
       { 
     pshiftWALD(q , mu = mu[i], sigma = sigma[i], nu=nu[i])   
       }
     #-------------------------------------------------------
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))      
    if (any(nu < 0))  stop(paste("nu must be positive", "\n", "")) 
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
# Shifted Wald random function adapted from pp. 79-80, Dagpunar, J. (1988). 
# Principles of Random Variate Generation. Clarendon Press, Oxford.
rwald2 <- function(n,mu,sigma,nu=0) {
  if(length(n)>1) n <- length(n);
  if(length(mu)>1 && length(mu)!=n) mu <- rep(mu,length=n)
  if(length(sigma)>1 && length(sigma)!=n) lambda <- rep(sigma,length=n)
  y2 <- rchisq(n,1); y2onm <- y2/mu; u <- runif(n)
  r1 <- (2*sigma + y2onm - sqrt(y2onm*(4*sigma+y2onm)))/(2*mu)
  r2 <- (sigma/mu)^2/r1; ifelse(u < sigma/(sigma+mu*r1), nu+r1, nu+r2)
}


#----------------------------------------------------------------------------------------  
rshiftWALD <- function(n, mu=1, sigma=1, nu=0, ...) # OK
  { 
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", "")) 
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qWALD(p,mu=mu,sigma=sigma, nu=nu, ...)
    r
  }
#----------------------------------------------------------------------------------------  

