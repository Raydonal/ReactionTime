#########################
######## ----- THE SLASH-WEIBULL DISTRIBUTION
#########################



# this function generates the Slash-Weibull distribution
# proposed by J. F. Olivares-Pacheco, H. C. Cornide-Reyes
# and M. Monasterio in
# "Una extensión de la distribución Weibull de dos parámetros"
# Revista Colombiana de Estadistica, (2010), 33(2), 219-231


#the parameters are n= number of obervations to be generated
# alpha = ?????
# beta = ?????
# and q = ????? kurtosis???

rsw = function(n = 1, alpha = 1, beta = 0.5, q = 1) 
{
  x = rep(0,n)
  for (i in 1:n) {
    u = runif(1, min=0, max=1)
    w = rweibull(1, shape=beta, scale = alpha)
    x[i] = w / u^(1/q)
  }
  return(x)
}

# example
SW = rsw(100, 10, .5, 10)
plot(density(SW), bty="n")

# future work:
# to fit SW to reaction time data and compare it against other candidates
# (e.g., Ex-Gaussian, Weibull, Gumbel, etc)
# estimation of fit can be assessd by AIC
