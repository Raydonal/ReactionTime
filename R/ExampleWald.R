#----------------------------------------------------------------------------------------  
# Example
#
#
# WALD()# gives information about the default links for the normal distribution
library(gamlss)

# Load the implementation
source("WALD.R")

# Load the implementation
source("ShiftWald.R")


# Graphics of WALD distribution
par(mfrow=c(2,2), pty="s")
# Density Function
plot(function(x)dWALD(x, mu=0.18856,sigma=70.711), xlim=c(0,900), ylim=c(0, 0.005), main="Wald Density", xlab="x", ylab="f(x)")
# Cumulative Function Distribution
plot(function(x)pWALD(x, mu=0.18856,sigma=70.711), xlim=c(0,900), ylim=c(0, 1), main="Wald CDF", xlab="x", ylab="F(x)")
# Quantile Function
plot(function(x)qWALD(x, mu=0.18856,sigma=70.711), xlim=c(0,1), ylim=c(0, 900), main="Wald Quantile", xlab="x", ylab="Q(x)")
# Random Number generation
n=1000
y=rWALD(n, mu=0.18856,sigma=70.711)
library(MASS)
hist.scott(y, main="Random y - Wald")
lines(density(y))
dev.off()


# Graphics of shift WALD distribution
par(mfrow=c(2,2), pty="s")
# Density Function
plot(function(x)dshiftWALD(x, mu=0.18856,sigma=70.711, nu=625), xlim=c(700,1400), ylim=c(0, 0.005), main="ShiftWald Density", xlab="x", ylab="f(x)")
# Cumulative Function Distribution
plot(function(x)pshiftWALD(x, mu=0.18856,sigma=70.711, nu=625),xlim=c(700,1400), ylim=c(0, 1), main="ShiftWald CDF", xlab="x", ylab="F(x)")
# Quantile Function
plot(function(x)qshiftWALD(x, mu=0.18856,sigma=70.711, nu=625), xlim=c(0,1), ylim=c(700,1400), main="ShiftWald Quantile", xlab="x", ylab="Q(x)")
# Random Number generation
n=1000
y=rshiftWALD(n, mu=0.18856,sigma=70.711, nu=625)
library(MASS)
hist.scott(y, main="Random y - ShiftWald")
lines(density(y))
dev.off()








# Example with real data and covariates
data(rent)        

# Comparing WALD and IG for rent data (response)
par(mfrow=c(1,2))
histDist(R,family=WALD, data=rent, ylim=c(0,0.0015), main="Wald") 
histDist(R,family=IG, data=rent, ylim=c(0,0.0015),main="IG") 

histDist(R,family=shiftWALD, data=rent, ylim=c(0,0.0015),main="IG") 


# add some covariates
# Using IG
 fit.ig = gamlss(R~Fl+Sp+Sm+B+loc, family=IG, data=rent)

# Using Wald

#Not run: con<-gamlss.control(c.crit = 0.0001, n.cyc = 50, mu.step=0.1, iter=10000, gd.tol=10)

# Here broken. I don't understarnd
fit.wald = gamlss(R~1, family=WALD, data=rent)

#histDist(y, family=WALD)


# Graphics of the Shift WALD distribution
par(mfrow=c(2,2), pty="s")
# Density Function
plot(function(x)dshiftWALD(x, mu=0.18856,sigma=70.711, nu=-0.5), xlim=c(0.,900), ylim=c(0, 0.005), main="Shift Wald Density", xlab="x", ylab="f(x)")
# Cumulative Function Distribution
plot(function(x)pshiftWALD(x, mu=0.18856,sigma=70.711, nu=-0.5), xlim=c(0,900), ylim=c(0, 1), main="Shift Wald CDF", xlab="x", ylab="F(x)")
# Quantile Function
plot(function(x)qshiftWALD(x, mu=0.18856,sigma=70.711, nu=-0.5), xlim=c(0,1), ylim=c(0, 900), main="Shift Wald Quantile", xlab="x", ylab="Q(x)")
# Random Number generation
n=1000
y=rshiftWALD(n, mu=0.18856,sigma=70.711, nu=-0.5)
library(MASS)
hist.scott(y, main="Random y")
lines(density(y))



fit.wald = gamlss(R~1, family=shiftWALD, data=rent)


fit.wald = gamlss(R~Fl+Sp+Sm+B+loc, family=WALD, data=rent)
