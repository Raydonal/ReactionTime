# ---------------------------------------------------------------------------------------
# This function generates and fit the Ex-Wald distribution
# using GAMLSS framewok as function of
#
# mu > 0 = mu (rate of evidence accrual)
# a  > 0 = sigma (response threshold)
# t  > 0 = nu (rate of exponential distribution)
#
# We use the approach given by
# Heathcote, A. (2004). Fitting Wald and ex-Wald distributions to response time 
# data: An example using functions for the S-PLUS package. Behavior Research 
# Methods, Instruments, & Computers,36,678-694.
# Create by Raydonal Ospina, Aug 17, 2014 
# Modified by Raydonal  03/09/2014
# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
EXWALD <- function (mu.link="log", sigma.link="log", nu.link = "log") 
{
  
  # Check Link Fuctions
  mstats <- checklink("mu.link", "Ex-Wald", substitute(mu.link), c("inverse", "log", "identity", "own" ))
  
  dstats <- checklink("sigma.link", "Ex-Wald", substitute(sigma.link), c("inverse", "log", "identity", "own" ))
  
  vstats <- checklink("nu.link", "Ex-Wald", substitute(nu.link), c("inverse", "log", "identity", "own" ))
  
  # Structure for fit using GAMLSS
  structure(
    list(family = c("EXWALD", "Ex-Wald"),
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
         nu.dr = vstats$nu.eta, 
         
         dldm = function(y,mu,sigma) ((y/mu)^sigma -1)*(sigma/mu),
         d2ldm2 = function(mu,sigma) - sigma^2/mu^2,
         dldd = function(y,mu,sigma) 1/sigma - log(y/mu)*((y/mu)^sigma-1) ,
         d2ldd2 = function(sigma) -1.82368/sigma^2 ,
         d2ldmdd = function(mu) 0.422784/mu,
         G.dev.incr  = function(y,mu,sigma,...) -2*dWEI(y, mu ,sigma, log=TRUE), 
         rqres = expression(rqres(pfun="pWEI", type="Continuous", y=y, mu=mu, sigma=sigma)),
         mu.initial = expression( { log.Y.m <- log(y) 
                                    var.Y.v <- var(log(y))
                                    sd.Y.s <- 1.283/sqrt(var.Y.v)
                                    mu <- exp(log.Y.m + 0.5772/sd.Y.s)
         }),
         sigma.initial = expression({  var.logY <- var(log(y))
                                       s.Y.s <- 1.283/sqrt(var.logY) 
                                       sigma <-  rep(s.Y.s,length(y))}),
         mu.valid = function(mu) all(mu > 0) , 
         sigma.valid = function(sigma)  all(sigma > 0), 
         y.valid = function(y)  all(y > 0)
    ),
    class = c("gamlss.family","family"))
}



############################## AUXILIARY FUNCTIONS ######################################
# ---------------------------------------------------------------------------------------
# Series approximation to the real (u) and imaginary (v) parts of complex 
# error function, erf(x + iy). Approximation can fail if firstblock > 20 
# but arguments allow forcing. Blocks loop will only process up to  
# maxseries terms.
uandv <- function(x,y,firstblock=20,block=0,tol=.Machine$double.eps^(2/3),maxseries=20)
{
  twoxy <- 2*x*y; xsq <- x^2; iexpxsqpi <- 1/(pi*exp(xsq))
  sin2xy <- sin(twoxy); cos2xy <- cos(twoxy)
  nmat <- matrix(rep((1:firstblock),each=length(x)),nrow=length(x))
  nsqmat <- nmat^2; ny <- nmat*y; twoxcoshny <- 2*x*cosh(ny)
  nsinhny <- nmat*sinh(ny); nsqfrac <- (exp(-nsqmat/4)/(nsqmat + 4*xsq))
  u <- (2*pnorm(x*sqrt(2))-1)+iexpxsqpi*(((1-cos2xy)/(2*x))+2*
                                           ((nsqfrac*(2*x-twoxcoshny*cos2xy+nsinhny*sin2xy))%*%rep(1,firstblock)))
  v <-iexpxsqpi*((sin2xy/(2*x))+2*((nsqfrac*(twoxcoshny*sin2xy+
                                               nsinhny*cos2xy))%*%rep(1,firstblock)))
  n <- firstblock;	converged <- rep(F,length(x))
  repeat {
    if ((block < 1) || (n >= maxseries)) break
    else {
      if ((n + block) > maxseries) block <- (maxseries - n)
      nmat <-matrix(rep((n+1):(n+block),each=sum(!converged)),
                    nrow=sum(!converged))
      nsq <- nmat^2; ny <- nmat*y[!converged];
      twoxcoshny <- 2*x[!converged]*cosh(ny); nsinhny <- nmat*sinh(ny)
      nsqfrac <- (exp(-nsq/4)/(nsq + 4*xsq[!converged]))
      du <- iexpxsqpi[!converged]*((2*nsqfrac*(2*x[!converged]-
                                                 twoxcoshny*cos2xy[!converged]+nsinhny*sin2xy[!converged]))
                                   %*%rep(1,block))
      dv <-iexpxsqpi[!converged]*((2*nsqfrac*(twoxcoshny*sin2xy[!converged]+
                                                nsinhny*cos2xy[!converged]))%*%rep(1,block))
      u[!converged] <- u[!converged] + du;
      v[!converged] <- v[!converged] + dv
      converged[!converged] <- ((du < tol) & (dv < tol))
      if (all(converged)) break
    }
  }
  cbind(u,v)
}
# ---------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------
# real part of w(z) = exp(-z^2)[1-erf(-iz)]
rew <- function(x,y,...) {
  uv <- uandv(y,x,...)
  exp(y^2-x^2)*(cos(2*x*y)*(1-uv[,1])+sin(2*x*y)*uv[,2])
}
# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
# Chi square with bins of width bin, for dist = "wald" (shifted wald),  
# "exw" (Ex-Wald) and "exg" (ex-gaussian). Makes a pass upward through 
# half of bins aggregating if < minn then makes a backward pass until 
# all bins have f >=minn, rt is data, p is paramter set
chisq <- function(rt,p,bin=25, minn = 5, dist = "wald"){
  switch(dist, wald ={pfun <- pwald},exw={pfun <- pexw},exg={pfun <- pexg})
  if (length(p) == 2) p <- c(p,0)
  cutrt <- seq((floor(min(rt)/bin)*bin+bin),(ceiling(max(rt)/bin)*bin-bin),bin)
  fitp <- c(pfun(cutrt,p[1],p[2],p[3]),1)
  fitn <- (fitp-c(0,fitp[-length(fitp)]))*length(rt)
  cumn <- 0; ncut <- ceiling(length(fitn)/2)
  cutok <- vector(length=length(cutrt)); cutok <- rep(T,length(cutok))
  for (i in 1:ncut) {
    cumn <- cumn+fitn[i]
    if (cumn >= 5) cumn <- 0
    else cutok[i] <- F
  }
  cutrt <- cutrt[cutok]
  fitp <- c(pfun(cutrt,p[1],p[2],p[3]),1)
  fitn <- (fitp-c(0,fitp[-length(fitp)]))*length(rt)
  cutok <- cutok[(cutok!=F)]
  cumn <- 0; ncut <- min(c(ncut,(length(fitn)-1)))
  for (i in (length(fitn)):(length(fitn)-ncut)) {
    cumn <- cumn+fitn[i]
    if (cumn >= 5) cumn <- 0
    else cutok[i-1] <- F
  }
  cutrt <- cutrt[cutok]
  df <- length(cutrt)-length(p)
  if (df < 1) c(NA,NA,NA) 
  else {
    fitp <- c(pfun(cutrt,p[1],p[2],p[3]),1)
    fitn <- (fitp-c(0,fitp[-length(fitp)]))*length(rt)
    obs <- tabulate(cut(rt,c(-Inf,cutrt,Inf)))
    chival <- c(sum((fitn-obs)^2/fitn))
    c(chisq=chival,df=df,p=(1-pchisq(chival,df)))
  }
}
# ---------------------------------------------------------------------------------------

############################# WALD #############################################

# Shifted Wald density
dwald <- function(w,m,a,s=0) {
  w <- w - s; a*exp(-(a-m*w)^2/(2*w))/sqrt(2*pi*w^3)
}

# Shifted Wald cumulative density with protection against numerical error
pwald <- function(w,m,a,s=0) {
  w <- w - s; sqrtw <- sqrt(w); k1 <- (m*w-a)/sqrtw; k2 <- (m*w+a)/sqrtw
  p1 <- exp(2*a*m); p2 <- pnorm(-k2); bad <- (p1==Inf) | (p2==0); p <- p1*p2
  p[bad] <- (exp(-(k1[bad]^2)/2 - 0.94/(k2[bad]^2))/(k2[bad]*((2*pi)^.5)))
  p + pnorm(k1)
}

# Shifted Wald random function adapted from pp. 79-80, Dagpunar, J. (1988). 
# Principles of Random Variate Generation. Clarendon Press, Oxford.
rwald <- function(n,m,a,s=0) {
  if(length(n)>1) n <- length(n);
  if(length(m)>1 && length(m)!=n) m <- rep(m,length=n)
  if(length(a)>1 && length(a)!=n) lambda <- rep(a,length=n)
  y2 <- rchisq(n,1); y2onm <- y2/m; u <- runif(n)
  r1 <- (2*a + y2onm - sqrt(y2onm*(4*a+y2onm)))/(2*m)
  r2 <- (a/m)^2/r1; ifelse(u < a/(a+m*r1), s+r1, s+r2)
}

# Start point estimate for the Wald and Shifted Wald, based on first two moments 
# For Shifted Wlad assumes s = p*min(x), where x is a data vector
waldstpt <- function(x,shift=T,p = 0.9) {
  if (shift) {
    s <- p*min(x); x <- x - s;	m <- sqrt(mean(x)/var(x))
    a <- m*mean(x); c(m,a,s)	
  } else {
    m <- sqrt(mean(x)/var(x))
    a <- m*mean(x); c(m,a)
  }	
}







#----------------------------------------------------------------------------------------
dWEI<-function(x, mu=1, sigma=1, log=FALSE)
{ 
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
  fy <- dweibull(x, scale=mu, shape=sigma, log =log)
  fy 
}
#----------------------------------------------------------------------------------------
pWEI <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
{     
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
  cdf <- pweibull(q, scale=mu, shape=sigma, lower.tail = lower.tail, log.p = log.p)
  cdf
}
#----------------------------------------------------------------------------------------
qWEI <- function(p, mu=1, sigma=1,  lower.tail = TRUE, log.p = FALSE)
{ if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  q <- qweibull(p, scale=mu, shape=sigma, lower.tail = lower.tail)
  q
}
#----------------------------------------------------------------------------------------
rWEI <- function(n, mu=1, sigma=1)
{ if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
  r <-  rweibull(n, scale=mu, shape=sigma)
  r
}
#----------------------------------------------------------------------------------------
