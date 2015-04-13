GenerateSurvData <- function(n=5000,p=100,s=5){
  ### Generate data to test the Cox function
  
  ### INPUTS
  # n: sample size
  # p: number of covariates
  # s: sparsity index (number of relevant variables)
  
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  beta <- matrix(c(runif(s,min=.5,max=1),rep(0,p-s)))
  
  y <- mapply(function(x) rexp(1, rate=x), exp(-X%*%beta))
  censure <- rbinom(n,1,prob=.5)
  return(X=X,y=y,cens=censure)
}