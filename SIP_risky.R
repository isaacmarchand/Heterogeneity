{
  # 't' year survival probability for member age 'x' with Gompertz-Makeham law
  tpx <- function(t,x,lambda = 0, b = 10, m = 85){
    exp(-lambda*t-exp((x-m)/b)*(exp(t/b)-1))
  }
  
  # Annuity Function for someone age 'x' with Gompertz-Makeham law
  annuity <- function(age, rate = .02, b = 10, m = 85){
    years <- 0:(150-age) #probability after age 150 usually negligible
    v <- exp(-years*rate)
    v%*%tpx(years, age, b = b, m = m)
  }
  
  # For a pool with 2 groups, one scenario of all the cumulative adjustment (\Delta) 
  # and the approximation cumulative adjustment (\Delta) simultaneously until 
  # every member dies when the pool is invested in a balanced portfolio (50/50).
  # INPUT: nb of members in group 1 and 2, age of members in group 1 and 2, 
  # initial benefit of group 1 and 2, lower bound of adjustment for stability (\epsilon1),
  # pre simulated vector of return for the risky asset until youngest group reach 
  # 120 years old, risk free return rate and Gompertz law parameters for group 1 and 2.
  # OUTPUT: integer giving the number of years the cumulative adjustement remained stable.
  StabilityCalc2PopRisky <- function(nb1, nb2, age1, age2, benefit1,
                                     benefit2, epsilonDown = .1,
                                     portReturn = rep(.02,120), rfrate = .02,
                                     m1 = 85, b1 = 10,
                                     m2 = 85, b2 = 10){
    #compute hurdle rate as expected return
    hrate <- (1-.5)*rfrate+.5*.07
    
    #simul death
    survProb1 <- c(tpx(1,age1:120, m=m1, b=b1),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #no one survive to age 121
    survProb2 <- c(tpx(1,age2:120, m=m2, b=b2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #no one survive to age 121
    surv1 <- c(nb1,numeric(length(survProb1)))
    for (t in 1:length(survProb1)) {
      surv1[t+1] <- rbinom(1,surv1[t],survProb1[t])
    }
    surv2 <- c(nb2,numeric(length(survProb2)))
    for (t in 1:length(survProb2)) {
      surv2[t+1] <- rbinom(1,surv2[t],survProb2[t])
    }
    
    #preCompute the annuities
    annuities1 <- c(sapply(age1:120,function(a) annuity(a,rfrate, m=m1, b=b1))
                    ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    annuities2 <- c(sapply(age2:120,function(a) annuity(a,rfrate, m=m2, b=b2))
                    ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    
    #compute first benefit
    B1 <- benefit1
    B2 <- benefit2
    y <- B2/B1
    
    #create adjustment vector
    lenDelta <- max(length(survProb2),length(survProb1))
    delta <- rep(1,lenDelta)
    
    for (t in 1:lenDelta) {
      death1 <- surv1[t]-surv1[t+1]
      death2 <- surv2[t]-surv2[t+1]
      
      #if no one left in pool set stability as max value
      if((surv2[t+1]+surv1[t+1])==0){
        stability <- t-1
        break
      }
      
      #compute adjustment
      delta[t] <- ((surv1[t]*survProb1[t]*annuities1[t+1]
                    +surv2[t]*y*survProb2[t]*annuities2[t+1])
                   /(surv1[t+1]*annuities1[t+1]
                     +surv2[t+1]*y*annuities2[t+1]))*exp(portReturn[t+1]-hrate)
      
      #check if benefits are still stable
      if(prod(delta)<(1-epsilonDown)){
        stability <- t-1
        break
      }
    }
    stability
  }
  
  # Gives the SIP based on the probability of not being stable every year using 
  # linear interpolation.
  # INPUT: 'x' a vector of a ll the number of stable years simulated by 
  # replicating StabilityCalc2PopRisky()
  # OUTPUT: (double) the SIP value 
  smoothedQuantiles <- function(x,p){
    pn <-  sapply(1:max(x), function(i) sum(x>=i)/length(x)) ## good one
    y <- seq_along(pn)
    z <- 1-pn
    i <- min(which(z>=p))
    
    pStar <- (p - z[i-1])/(z[i]-z[i-1])
    SIP <- (y[i] - y[i-1])*pStar+ y[i-1]
    SIP
  }
  
  # Estimate the SIP for a pool with 2 groups or 1 groups if 'nb2' = 0
  # INPUT: nb of members in group 1 and 2, age of members in group 1 and 2, 
  # initial benefit of group 1 and 2, lower bound of adjustment for stability (\epsilon1),
  # pre simulated vector of return for the risky asset until youngest group reach 
  # 120 years old, risk free return rate and Gompertz law parameters for group 1 and 2.
  # OUTPUT: (double) the SIP value
  SIP2Pop_risky <- function(nb1, nb2, age1, age2, benefit1,
                            benefit2, epsilonDown=.3, beta=.95, nbSimul=10000,
                            avgRisky = .07, volRisky = .15, rfrate = .02,
                            m1 = 85, b1 = 10, m2 = 85, b2 = 10){
    
    #preSimulate portfolio returns
    normRates <- replicate(nbSimul,rnorm(120-min(age1, age2), avgRisky, volRisky))
    portReturn <- (1-.5)*rfrate+.5*normRates
    
    
    #simulate 'nbSimul' stability scenarios
    stability <- sapply(1:nbSimul, 
                        function(i) StabilityCalc2PopRisky(nb1,nb2,age1,age2,
                                                           benefit1,benefit2,
                                                           epsilonDown,
                                                           portReturn[,i],rfrate,
                                                           b1 = b1,b2 = b2,
                                                           m1 = m1,m2 = m2))
    SIP <- smoothedQuantiles(stability, 1-beta)
    
  }
}


####### Example Simulation #####
# compute the SIP of 2 groups of 100 members both age 65, as y goes from 0.1 to 10
# (using base parameter of function for everything else)

nbSimul <- 1000 # for fast computation but very noisy
epsilon1 <- .3 #to account for more volatile returns

nb1 <- 100
nb2 <- 100
age1 <- 65
age2 <- 65

y <- exp(seq(log(1/10),log(10), length.out = 101))

SIP <- sapply(y, function(x) SIP2Pop_risky(nb1,nb2,age1, age2, 1,x,
                                           nbSimul = nbSimul,
                                           epsilonDown = epsilon1))

# would then apply some sort of GAM smoothing on SIP vector
plot(log(y), SIP, type= "l")

