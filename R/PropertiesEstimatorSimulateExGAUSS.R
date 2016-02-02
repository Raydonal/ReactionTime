
# Simulating exGAUSS distribution to evaluate Estimators
# General program to split covariates (Future)
# Raydonal Ospina 2014


####
require(gamlss.dist)

## The next situations ha a problem estimation in GAMLSS
## We use moment method estimation
# EGd_1 = rexGAUS(n, 300, 20, 300)      ### very skewed
# EGd_2 = rexGAUS(n, 400, 20, 200)      #### mildly skewed
# EGd_3 = rexGAUS(n, 500, 20, 100)       ### low skewed
# EGd_3 = rexGAUS(n, 600, 20, 50)        ### more normal-like


# small size we have problems. Heuristic method as Swarm 
# optimization can be a solution.


 require(gamlss)

simulation <- function(
  
  nrep = 500,				# Monte Carlo Replics
  
  nsample = 100, #c(5, 10, 30, 50, 100, 300, 500),		# N sample
  
  # To compare  other models c("exGAUSS","LOGNO", "WEI", "NO"),  	# Family distributions
  type = "exGAUS", # More information ?exGAUs
  
  # Parameters to exGAUSS
  par.mu = 600, #c(1,1,-0.5,0.5),   #c(1,1,-0.5,0.5,0.25, -0.25),	# Parameters to mu
  par.sigma =20,  #c(2,1,0.5,0.5),   #0.25, 0.25),	# Parameters to sigma
  par.nu =50,  # c(-1,-1,-0.5,0.5), #,0.25, -0.25),# Parameters to nu
    
  # Link function to parameters submodels    
  family.mu = "identity",	   		#link to mu
  family.sigma ="log",			#Link to sigma
  family.nu = "log"  		#link to nu
)
{
  
  
  pm = proc.time()
  
  time.ini = format(Sys.time(), " %b %d %Y  %H:%M:%S ")
  cat("-------------------------------------------------------------------\n")
  
  cat("\nData and Hour Initial: ",time.ini)
  cat("\n")
  
  cat("R Version:", version$version.string)
  cat("\n")
  
  cat("N REP: ", nrep)
  cat("\n")
  
  cat("-------------------------------------------------------------------\n")
  
  cat("PARAMETERS for nu: ", par.nu)
  cat("\n")
  
  cat("PARAMETERS for mu:    ", par.mu)
  cat("\n")
  
  cat("PARAMETERS for sigma:   ", par.sigma)
  cat("\n")
  
  
  type = match.arg(type)
  
  param = c(par.mu, par.sigma, par.nu)
  
  
  for(i in 1:length(nsample))
  {
    cat("-------------------------------------------------------------------\n")
    
    cat("N SAMPLE: ", nsample[i])
    cat("\n")
    
    tam.sample = nsample[i] #tamanho de amostra
    
    #-----------------------------------------------------------------------------------------
    #Generate covariate to linear preditor to mu
    x0 = rep(1,tam.sample)		# Intercepto
    x1 = rnorm(tam.sample)		# covariada X1 para mu uma Normal padrão
    x2 = rpois(tam.sample,1)	# covariada X2 para mu uma Poisson(1)
    x3 = rbinom(n=tam.sample, prob=0.2, size=5)	# covariada X3 para mu uma Binomial(5,0.2)
    x4 = rbinom(n=tam.sample, prob=0.3, size=5)	# covariada X4 para mu uma Binomial(5,0.3)
    x5 = runif(tam.sample)		# covariada X5 para mu uma Uniforme padrão
    
    # matrix design to mu
    x.mu = as.matrix(x0) #cbind(x0,x1, x2, x3) 
    
    
    # linear predictor to mu
    pred.mu = x.mu %*% par.mu	 
    
    #link function to mu
    link.mu = make.link.gamlss(family.mu)
    
    #mu as function to predictor 
    mu = link.mu$linkinv(pred.mu)
    
    cat("mean mu:    ", mean(mu))
    cat("\n")	 
    
    #-----------------------------------------------------------------------------------------
    
    #Generate covariate to linear preditor to sigma
    v0 = rep(1,tam.sample)		# Intercepto
    v1 = rnorm(tam.sample)		# covariada V1 para nu uma Normal padrão
    v2 = rpois(tam.sample,1)	# covariada V2 para numu uma Poisson(1)
    v3 = rbinom(n=tam.sample, prob=0.2, size=5)	# covariada V3 para nu uma Binomial(5,0.2)
    v4 = rbinom(n=tam.sample, prob=0.3, size=5)	# covariada V4 para nu uma Binomial(5,0.3)
    v5 = runif(tam.sample)		# covariada V5 para nu uma Uniforme padrão
    
    # matrix design to sigma	 
    #	v.sigma = cbind(v0, v1, v2, v3, v4,v5) 
    
    v.sigma = as.matrix(v0) #cbind(v0, v1, v2, v3)
    
    #linear predictor to sigma
    pred.sigma = v.sigma %*% par.sigma
    
    #link function to sigma
    link.sigma = make.link.gamlss(family.sigma)
    
    #sigma as function to predictor #link.sigma$linkinv(pred.sigma)
    sigma = pred.sigma 
    
    cat("mean sigma:  ", mean(sigma))
    cat("\n")	 
    
    
    #-----------------------------------------------------------------------------------------
    
    #Generate covariate to linear preditor to nu
    z0 = rep(1,tam.sample)  	# Intercepto
    z1 = rnorm(tam.sample)		# covariada Z1 para nu uma Normal padrão
    z2 = rpois(tam.sample,1)	# covariada Z2 para numu uma Poisson(1)
    z3 = rbinom(n=tam.sample, prob=0.2, size=5)	# covariada Z3 para nu uma Binomial(5,0.2)
    z4 = rbinom(n=tam.sample, prob=0.3, size=5)	# covariada Z4 para nu uma Binomial(5,0.3)
    z5 = runif(tam.sample)		# covariada Z5 para nu uma Uniforme padrão
    
    # matrix design to nu
    z.nu = as.matrix(z0) #cbind(z0, z1, z2, z3)
    
    #preditor linear de nu
    pred.nu = z.nu %*% par.nu	 
    
    #link function to nu
    link.nu = make.link.gamlss(family.nu)
    
    #nu as function to predictor  link.nu$linkinv(pred.nu)
    nu = pred.nu 
    
    cat("mean nu: ", mean(nu))
    cat("\n")	 
    
    #-----------------------------------------------------------------------------------------
    # Save the true parameters in file .Rdata
    vaucc = cbind(as.vector(mu), as.vector(nu), as.vector(sigma))
    save(vaucc,    file=paste("fv.true",    as.character(nsample[i]), ".Rdata", sep=""))
    
    
    #-----------------------------------------------------------------------------------------
    
    set.seed(121) # fixed the seed
    
    
    # Auxiliar arrays
    
    # Estimators MLE
    mle.mu = NULL
    mle.nu = NULL
    mle.sigma = NULL
    
    
    # Fitted values
    fv.mu    = NULL
    fv.nu = NULL
    fv.sigma   = NULL
    fv.mean  = NULL
    
    
    falla = 0
    
    for(j in 1:nrep+falla)
    {
      
      #	if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      #	runif(1)
      #	seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      
      # Responde variable
      ydata = y = mapply(rexGAUS,  1, mu, sigma, nu) 
      
      #	cat("-------------------------------------------------------------------\n")
      
      
      cond.rep = glim.control(trace=FALSE) # not impression GLIM iteration		
      cont.fit.rep = gamlss.control(trace=FALSE, mu.step = 0.1, sigma.step =0.1, nu.step =0.1)
      
      #		fit.rep=gamlss(ydata~(-1+x0+x1+x2+x3+x4+x5), sigma.formula=~(-1+v0+v1+v2+v3+v4+v5), nu.formula=~(-1+z0+z1+z2+z3+z4+z5), family=type, 					control=cont.fit.rep, i.control = cond.rep)
      
      # Fit model using GAMLSS - MLE are performed via RS() or GC() algorithms ( Score Fisher method using backfitting)
      fit.rep=gamlss(ydata~-1+x.mu, sigma.formula=~-1+v.sigma, nu.formula=~-1+z.nu, family=type, control=cont.fit.rep, i.control = cond.rep)
      
      
      
      if(fit.rep$converged =="TRUE")
      {
        
        #Estimates of parameters
        mle.est.mu 	    = fit.rep$mu.fv[1] #fit.rep$mu.coefficients
        mle.est.nu 	    = fit.rep$nu.fv[1] #fit.rep$nu.coefficients
        mle.est.sigma 	= fit.rep$sigma.fv[1] #fit.rep$sigma.coefficients
        
        
        mle.fv.mu 	= fit.rep$mu.fv
        mle.fv.nu 	= fit.rep$nu.fv
        mle.fv.sigma 	= fit.rep$sigma.fv
        
                
        
        #Saving the estimates of parameters in arrays (used for grasigmacs)
        mle.mu	  = cbind(mle.mu, mle.est.mu)
        mle.nu    = cbind(mle.nu, mle.est.nu)
        mle.sigma	= cbind(mle.sigma, mle.est.sigma)
        
        
        #Saving the estimates of parameters in arrays (used for grasigmacs)
        fv.mu	  = cbind(fv.mu, mle.fv.mu)
        fv.nu  = cbind(fv.nu, mle.fv.nu)
        fv.sigma	  = cbind(fv.sigma, mle.fv.sigma)
        
        # Bias of estimators
        biasmle.mu    = mle.mu - par.mu
        biasmle.nu    = mle.nu - par.nu
        biasmle.sigma = mle.sigma - par.sigma
        
        
        # EQM of estimators
        eqmmle.mu = (mle.mu - par.mu)^2
        eqmmle.nu = (mle.nu - par.nu)^2
        eqmmle.sigma= (mle.sigma - par.sigma)^2         
        
        
      }
      else
      {
        falla = falla+1
        j=j-1 
      }
      
      # Relative Bias
      biasrel.mu     = biasmle.mu  *  (1/par.mu)
      biasrel.nu  = biasmle.nu * (1/par.nu)
      biasrel.sigma    = biasmle.sigma * (1/par.sigma)
      
    }
    
    # Save estimates to graphics depending the nsample
    save(mle.mu,    file=paste("mle.mu",    as.character(nsample[i]), ".Rdata", sep=""))
    save(mle.nu, file=paste("mle.nu", as.character(nsample[i]), ".Rdata", sep=""))
    save(mle.sigma,   file=paste("mle.sigma",   as.character(nsample[i]), ".Rdata", sep=""))
    
    
    save(fv.mu,    file=paste("mle.fv.mu",    as.character(nsample[i]), ".Rdata", sep=""))
    save(fv.nu, file=paste("mle.fv.nu", as.character(nsample[i]), ".Rdata", sep=""))
    save(fv.sigma,   file=paste("mle.fv.sigma",   as.character(nsample[i]), ".Rdata", sep=""))
    
    
    #-----------------------------------------------------------------------------------------------
    # Mean 
    mean.mle.mu	= apply(mle.mu,1, mean)
    
    # Bias
    bias.mle.mu	= apply(biasmle.mu ,1, mean) 
    
    # Standar error
    std.mle.mu	= apply(mle.mu,1, sd)
    
    # EQM
    eqm.mle.mu	=  sqrt(apply(eqmmle.mu,1, mean))  
    
    # Relative bias
    biasrel.mle.mu	= apply(biasrel.mu ,1, mean)
    
    #-----------------------------------------------------------------------------------------------------	
    mean.mle.nu      = apply(mle.nu,1, mean)
    
    bias.mle.nu      = apply(biasmle.nu ,1, mean)  
    
    std.mle.nu      = apply(mle.nu,1, sd)
    
    eqm.mle.nu      =   sqrt(apply(eqmmle.nu,1, mean))  
    
    biasrel.mle.nu      = apply(biasrel.nu ,1, mean)
    
    #-----------------------------------------------------------------------------------------------------	
    mean.mle.sigma      = apply(mle.sigma,1, mean)
    
    bias.mle.sigma      = apply(biasmle.sigma ,1, mean)  
    
    std.mle.sigma      = apply(mle.sigma,1, sd)
    
    eqm.mle.sigma      =   sqrt(apply(eqmmle.sigma,1, mean))  
    
    biasrel.mle.sigma      = apply(biasrel.sigma ,1, mean)
    
    
    #-----------------------------------------------------------------------------------------------------	
    
    par.chose.mu = names(mle.est.mu)
    
    p.mle.mu  = as.data.frame(cbind( round(mean.mle.mu    , digits=5),
                                     round(bias.mle.mu    , digits=5),
                                     round(std.mle.mu     , digits=5),
                                     round(eqm.mle.mu     , digits=5),
                                     round(biasrel.mle.mu , digits=5)
    ))
    #-----------------------------------------------------------------------------------------------------
    
    par.chose.nu = names(mle.est.nu)
    
    p.mle.nu  = as.data.frame(cbind( round(mean.mle.nu    , digits=5),
                                        round(bias.mle.nu    , digits=5),
                                        round(std.mle.nu     , digits=5),
                                        round(eqm.mle.nu     , digits=5),
                                        round(biasrel.mle.nu , digits=5)
    ))
    #-----------------------------------------------------------------------------------------------------
    par.chose.sigma = names(mle.est.sigma)
    
    p.mle.sigma  = as.data.frame(cbind( round(mean.mle.sigma    , digits=5),
                                      round(bias.mle.sigma   , digits=5),
                                      round(std.mle.sigma    , digits=5),
                                      round(eqm.mle.sigma    , digits=5),
                                      round(biasrel.mle.sigma , digits=5)
    ))
    
    #-----------------------------------------------------------------------------------------------------
    
    cat("No. de falhas:", falla)
    cat("\n")
    
    cat("-------------------------------------------------------------------\n")	
    cat("Mle estimators for mu:\n")
    cat("\n")
    colnames(p.mle.mu) <-  c("Mean", "Bias", "Std. Error", "EQM", "Bias Rel.")
    print(p.mle.mu, print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-------------------------------------------------------------------\n")	
    
    cat("Mle estimators for nu:\n")
    cat("\n")
    colnames(p.mle.nu) <-  c("Mean", "Bias", "Std. Error", "EQM", "Bias Rel.")
    print(p.mle.nu, print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-------------------------------------------------------------------\n")	
    
    
    cat("Mle estimators for sigma:\n")
    cat("\n")
    colnames(p.mle.sigma) <-  c("Mean", "Bias", "Std. Error", "EQM", "Bias Rel.")
    print(p.mle.sigma, print.gap = 2, quote = FALSE)
    cat("\n")
    cat("-------------------------------------------------------------------\n")	
    
    time.sample = format(Sys.time(), " %b %d %Y  %H:%M:%S ")
    cat("\nDate and parcial hour :" ,time.sample , "\n")
    cat("\n")
    
    gc() # clean memory
    
    
  }
  time.fin = format(Sys.time(), " %b %d %Y  %H:%M:%S ")
  cat("\nDate and final Hour : ",time.fin , "\n")
  
  cat("\nTotal time  in minutes:", ((proc.time()-pm)[1])/(3600), "\n")
}


simulation()
