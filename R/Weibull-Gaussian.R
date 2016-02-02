## x: Reaction time as integer (vector is acceptable). 
## m: Mean of Weibull and Gaussian component.
## s: Size parameter of Gaussian component.
## n: Size parameter of Weibull component.
## k: Shape parameter of Weibull component.


dweibulgaus <- function(x, m, s, n, k){
  ord <- x
    if(max(x) > 4000){
    x <- 1:(max(x)+100)
  }else{
    x <- 1:4000
  }

  tl <- length(x)
  part1 <- dweibull(x, k, n)
  part2 <- dnorm(x, m, s)
  convolve(part1, part2, conj=F)[ord]
}


pweibulgaus <- function(t, m, s, n, k){
  ord <- t
  if(max(t) > 4000){
    t <- 1:(max(t)+100)
  }else{
    t <- 1:4000
  }

  tl <- length(t)
  part1 <- dweibull(t, k, n)
  part2 <- dnorm(t, m, s)
  cumsum(convolve(part1, part2, conj=F))[ord]
}


rweibulgaus <- function(x, m, s, n, k){
  rweibull(x, k, n)+rnorm(x, m, s)
}

