# required libraries
require(gamlss.dist)
require(gamlss)

### data given by authors
ex.data =  c(0.542182227221597, 0.085864266966629, 0.999625046869142, 0.638545181852268,
            0.194225721784777, 2.22347206599175, 0.092613423322085, 0.999250093738283,
            2.40119985001875, 0.854143232095988, 0.075740532433446, 2.61267341582302)



### Fitting a Normal, Weibull, Log-Normal, and Ex-Gaussian model to the data and plotting the fittings
op = par(mfrow = c(2, 2))
mod.NO = histDist(ex.data, "NO", density = TRUE, main = "Normal-fitted")
mod.WEI = histDist(ex.data, "WEI", density = TRUE, main = "Weibull-fitted")
mod.LOGNO = histDist(ex.data, "LOGNO", density = TRUE, main = "Log-Normal-fitted")
mod.ExG = histDist(ex.data, "exGAUS", density = TRUE, main = "EGd-fitted")

### Looking at the AIC of each model
GAIC(mod.NO,mod.WEI,mod.LOGNO)

### Extracting the parameter estimates of the ExGaussian fitted model
ex.mu = mod.ExG$mu.fv[1]
ex.sigma = mod.ExG$sigma.fv[1]
ex.nu = mod.ExG$nu.fv[1]


