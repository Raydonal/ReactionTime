############################# AUXILIARY FUNCTIONS ######################################

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

# real part of w(z) = exp(-z^2)[1-erf(-iz)]
rew <- function(x,y,...) {
	uv <- uandv(y,x,...)
	exp(y^2-x^2)*(cos(2*x*y)*(1-uv[,1])+sin(2*x*y)*uv[,2])
}

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

# -loglikelihood of the Shifted Wald, parm = (m,a,s), x is a raw data vector
negllswald <- function(p,x) {
	-sum(log(dwald(x,p[1],p[2],p[3])))
}

# -loglikelihood of the Wald, parm = (m,a), x is a raw data vector
negllwald <- function(p,x) {
	-sum(log(dwald(x,p[1],p[2])))
}


# gradient and hessian of -loglikelihood of shifted Wald
# p = (m,a,s), x is a raw data vector
negllswaldg <- function(p,x){
	sw <- x - p[3]
	grad <- cbind(p[2] - p[1]*sw,                 # dm
					1/p[2] - p[2]/sw + p[1],         # da
					(p[1]^2 + (3 - p[2]^2/sw)/sw)/2) # ds
	-colSums(grad)
}


# gradient and hessian of -loglikelihood of shifted Wald
# p = (m,a,s), x is a raw data vector
negllswaldh <- function(p,x){
	sw <- x - p[3]
	hes <- cbind(-sw,                                    # hessian
	             1, -1/p[2]^2 - 1/sw,                     # lower
	             p[1], -p[2]/sw^2, (1.5-p[2]^2/sw)/sw^2) # triangle
 matrix(-colSums(hes)[c(1,2,4,2,3,5,4,5,6)],ncol=3,byrow=T)
}


# gradient and hessian of -loglikelihood of Wald
# p = (m,a), x is a raw data vector
negllwaldg <- function(p,x){
	grad <- cbind(p[2] - p[1]*x,               # dm
					1/p[2] - p[2]/x + p[1])       # da
  -colSums(grad)
}

# gradient and hessian of -loglikelihood of Wald
# p = (m,a), x is a raw data vector
negllwaldh <- function(p,x){
	hes <- cbind(-x,                           # hessian lower
	             1, -1/p[2]^2 - 1/x)           # triangle
 matrix(-colSums(hes)[c(1,2,2,3)],ncol=2,byrow=T)
}



# Fit the Wald and Shifted Wald and calculate parameter standard errors 
# correlations, and Chisquare, rt is a raw data vector
fitwald <- function(rt,shift=T,p = 0.9,start=waldstpt(rt,shift,p)){
	if (shift) fit <- nlminb(start=start,lower=c(0,0,0),upper=c(Inf,Inf,min(rt)),
					objective=negllswald,gradient=negllswaldg,hessian=negllswaldh,x=rt)
	else fit <- nlminb(start=start,lower=c(0,0),
					objective=negllwald,gradient=negllwaldg,hessian=negllwaldh,x=rt)
	cat(paste(fit$message,"\n"))
	fit$chisq <- chisq(rt=rt,p=fit$par,dist="wald")
	fit
}

############################# EX WALD #############################################

# Ex-Wald density 
dexw <- function(r,m,a,t) {
	k <- (m^2 - (2/t))
	if (k < 0) {
		return(exp(m*a - (a^2)/(2*r) - r*(m^2)/2)*
              rew(sqrt(-r*k/2),a/sqrt(2*r))/t)
	} else {
		k <- sqrt(k); return(pwald(r,k,a)*exp(a*(m-k) - (r/t))/t)
	}
}

# Ex-Wald cumulative density
pexw <- function(r,m,a,t) {
	pwald(r,m,a) - t*dexw(r,m,a,t)
}

# Ex-Wald random function
rexw <- function(n,m,a,t) {
	rwald(n,m,a) + rexp(n,1/t)
}

# Start point estimate for the Ex-Wald, based on first two moments 
# assuming t = p*stdev(x), where x is a raw data vector
exwstpt <- function(x,p = 0.5) {
	t <- p*sd(x); m <- sqrt((mean(x)-t)/(var(x)-t^2))
	a <- m*(mean(x)-t); c(m,a,t)	
}

# -loglikelihood of Ex-Wald, param = (m,a,t), x is a data vector
negllexw <- function(p,x) {
	-sum(log(dexw(x,p[1],p[2],p[3])))
}

# Fit the Ex-Wald and calculate parameter standard errors and
# correlations, rt is a raw data vector
fitexw <- function(rt,p=0.5,start=exwstpt(rt,p),scaleit=T){
	if (scaleit) fit <- nlminb(start=start,lower=c(1e-8,1e-8,1),objective=negllexw,x=rt,scale=(1/start),
									control = list(eval.max=400,iter.max=300))
	else fit <- nlminb(start=start,lower=c(1e-8,1e-8,1),objective=negllexw,x=rt,
							control = list(eval.max=400,iter.max=300))
	cat(paste(fit$message,"\n"))
	fit$chisq <- chisq(rt=rt,p=fit$par,dist="exw")
	fit
}

############################# EX-GAUSSIAN #############################################

# ex-gaussian density, m s t = mu, sigma, tau
dexg <- function(x,m,s,t){
   exp(((m-x)/t)+0.5*(s/t)^2)*pnorm(((x-m)/s)-(s/t))/t
}

# ex-gaussian cumulative density at x, m,s,t = mu, sigma, tau
pexg <- function(x,m,s,t){
	rtsu <- (x-m)/s
	pnorm(rtsu) - exp((s^2/(2*t^2))-((x-m)/t))*pnorm(rtsu - (s/t))
}

# generates n random samples from the ex-gaussian, p = mu, sigma, tau
rexg <- function(n,m,s,t){
	rexp(n,rate=t) + rnorm(n,mean=m,sd=s)	
}

# Three moment start point estimate for the ex-gaussian
# If this fails both tau and sig are estimated heuristically 
# and mu determined by equating the mean.
exgstartpt <- function(rt, p = 0.8){
	m1 <- mean(rt)
	m2 <- var(rt)
	m3 <- sum((rt-m1)^3)/length(rt-1)
	tau <- m3^(1/3)/2
	sig <- (m2 - tau^2)^.5
	mu <- m1 - tau
	if (any(c(mu,sig,tau) <= 0)) {
		tau <- p*(m2^.5)
		sig <- tau*(1-p^2)^.5
		mu <- m1 - tau
	}
	c(mu,sig,tau)
}

negllexg <- function(p,x){
   -sum((((p[1]-x)/p[3])+0.5*(p[2]/p[3])^2) + log(pnorm(((x-p[1])/p[2])-(p[2]/p[3]))/(p[3]*(2*pi)^.5)))
}


# Gradient of -loglikelihood for ex-gaussian at t, x = mu, sigma, tau
negllexggrad <- function(p,x){
	exggrad <- function(p,x){
		exg <-	exp(((p[1]-x)/p[3])+0.5*(p[2]/p[3])^2)*pnorm(((x-p[1])/p[2])-(p[2]/p[3]))/p[3]
   		rtpiex <- (1/(p[3]*sqrt(2*pi)))*exp(-0.5*((x-p[1])/p[2])^2)
		mug <- exg/p[3] -rtpiex/p[2]
   		sigg <- exg*(p[2]/p[3]^2) - ((1/p[3])+((x-p[1])/p[2]^2))*rtpiex
	   taug <- ((-1 - (p[2]/p[3])^2 + (x-p[1])/p[3])/p[3])*exg + (p[2]/(p[3]^2))*rtpiex
   		cbind(mug,sigg,taug)
	}
	-1*colSums(exggrad(p,x)/(exp(((p[1]-x)/p[3])+0.5*(p[2]/p[3])^2)*pnorm(((x-p[1])/p[2])-(p[2]/p[3]))/p[3]))
}

	
fitexg <- function(rt,p=0.8,start=exgstartpt(rt,p)){
	fit <- nlminb(start=start,lower=c(0,0,0),objective=negllexg,gradient=negllexggrad,x=rt)
	cat(paste(fit$message,"\n"))
	fit$chisq <- chisq(rt=rt,p=fit$par,dist="exg")
	fit
}

#### paper >>>>>>>>>>>>>


#Heathcote, A. (2004) Fitting the Wald and Ex-Wald Distributions to Response Time Data. 
#Behaviour Research Methods, Instruments & Computers, 36, 678-694.

