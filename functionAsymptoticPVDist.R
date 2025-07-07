# Surival prob function with Markham law
tpx <- function(t,x,A = .0001, B = .0001, c =1.089){
  exp(-A*t-(B/log(c))*c^x*(c^t-1))
}

# AnnuityFunction
annuity <- function(age, rate){
  years <- 0:(150-age) #probability after age 150 usually negligible
  exp(-years*rate)%*%tpx(years, age)
}

# Compute Variance of Asymptotic gain distribution based on Eq.26 in Overleaf Report
asymptoVariance <- function(currentAge=65, rate = .02, innitialAsset = 100000,
                            maxAge=120){
  survProb <- c(tpx(0:(maxAge-currentAge),currentAge),0)
  B0 <- innitialAsset/annuity(currentAge,rate)
  deathProb <- 1-survProb
  discount <- exp(-rate*((currentAge:(maxAge+1))-currentAge))
  
  #Compute sum matrix where we are only interested in lower triangle which 
  #represent the min in the variance eq.
  #and the diag which is variance
  sumMatrix <- (discount*survProb)%*%t(discount*deathProb)
  B0^2*sum(diag(sumMatrix),2*sumMatrix[lower.tri(sumMatrix)])
}


# Quick Test:
# Should be around 48 000 for initial age and asset of 65 and $100 000 
# (based on simul with 100 000 participants )
sqrt(asymptoVariance(65))
