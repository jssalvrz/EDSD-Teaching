#------------------------------------------------------------------------------------------------------------
# Functions to fit a Gamma - Gompertz
#-----------------------------------------------------------------------------------------------------------

# define functions
# force of mortality
mu <- function(pars,x) {
  a <- pars[1]
  b <- pars[2]
  out <- a*exp(b*x)
  return(out)
}

# cumulative baseline hazard
M0 <- function(pars,x){
  a <- pars[1]
  b <- pars[2]
  out <- a/b*(exp(b*x)-1)
  return(out)
}

# population force of mortality, gamma denotes the
# variance of the Gamma distribution
muGammaNoCov <- function(pars,x) {
  gamma <- pars[3]
  num <- mu(pars,x)
  denom <- 1+gamma*M0(pars,x)
  out <- num/denom
  return(out)
}

lnL <- function(pars,x,D,E,ageInterval) {
  print(pars)
  mubar <- muGammaNoCov(pars,x=(x-ageInterval/2))
  out <- sum(D*log(mubar)-E*mubar)
  print(out)
  return(out)
}



# derivatives of the force of mortality function
dMuGamma <- function(pars,x) {
  a <- pars[1]
  b <- pars[2]
  gamma <- pars[3]
  da <- (b^2*exp(b*x))/(b + a*(exp(b*x)-1)*gamma)^2
  db <- (a*exp(b*x)*(b^2*x + a*gamma*(exp(b*x)-1-b*x)))/
    (b + a*(exp(b*x)-1)*gamma)^2
  dgamma <- -((a^2*b*exp(b*x)*(exp(b*x)-1))/(b +
                                               a*(exp(b*x)-1)*gamma)^2)
  return(as.matrix(data.frame(da=da,db=db,dgamma=dgamma)))
}



# get standard error by the delta method
getSE <- function(pars,hessian,derivFunc,x) {
  sqrt(derivFunc(pars,x)%*%hessian%*%t(derivFunc(pars,x)))
}


GG <- function(D,E){
  # estimate parameters
  opt <-  optim(lnL,x=1:length(D),D=D,E=E,control=list(fnscale=-1),
                ageInterval=1,
                par=c(0.06,0.14,0.2), hessian = T)
  
  # calculate standard errors
  matVC <- solve(-opt$hessian) # variance-covariance matrix
  ages  <-  1:length(D)-0.5  # vector of age groups
  ses <- getSE(pars=opt$par,hessian=matVC,derivFunc=dMuGamma,x=ages)
  SEs   <- rep(NA,times=length(ages)) # standard errors
  for(i in 1:length(SEs)) {
    SEs[i] <- getSE(pars=opt$par,hessian=matVC,derivFunc=dMuGamma,x=ages[i])
  }
  
  estimate <- muGammaNoCov(pars=opt$p,x=ages)
  lower    <- estimate-qnorm(0.99)*SEs
  upper    <- estimate+qnorm(0.99)*SEs
  a <- opt$par[1]
  b <- opt$par[2]
  gamma <- opt$par[3]
  ll    <- opt$value
  AIC   <- 2 * 3 - 2 *ll #3 parameters
  
  return(data.frame(Age = ages, hx = estimate, hx.low =lower,hx.upper=upper, a, b, gamma, AIC, loglik = ll))
  
}


#------------------------------------------------------------------------------------------------------------
# Functions to fit a Gamma - Gompertz via the Mode
#-----------------------------------------------------------------------------------------------------------

# force of mortality
Mmu <- function(pars,x) {
  M <- pars[1]
  b <- pars[2]
  out <-b*exp(b*(x-M)) #a*exp(b*x)
  return(out)
}

# cumulative baseline hazard
MM0 <- function(pars,x){
  M <- pars[1]
  b <- pars[2]
  out <- exp(-b*M)*(exp(b*x)-1)   #a/b*(exp(b*x)-1)
  return(out)
}

# heterogeneous (marginal) force of mortality, gamma denotes the
# variance of the Gamma distribution
MmuGammaNoCov <- function(pars,x) {
  gamma <- pars[3]
  num <- Mmu(pars,x)
  denom <- 1+gamma*MM0(pars,x)
  out <- num/denom
  return(out)
}

# likelihood
MlnL <- function(pars,x,D,E,ageInterval) {
  print(pars)
  mubar <- MmuGammaNoCov(pars,x=(x-ageInterval/2))
  out <- sum(D*log(mubar)-E*mubar)
  print(out)
  return(out)
}
