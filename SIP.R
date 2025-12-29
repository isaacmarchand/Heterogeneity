####### Create functions #####
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
  # every member dies.
  # INPUT: nb of members in group 1 and 2, age of members in group 1 and 2, 
  # initial benefit of group 1 and 2, lower bound of adjustment for stability (\epsilon1),
  # and Gompertz law parameters for group 1 and 2.
  # OUTPUT: list containing the number of stable years before first unstable \Delta ($stability)
  # and vactor of which year were stable based on \DeltaStar
  Stability2Gr2Variate <- function(nb1, nb2, age1, age2, benefit1,
                                        benefit2, epsilonDown = .1,
                                        m1 = 85, b1 = 10,
                                        m2 = 85, b2 = 10){
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
    annuities1 <- c(sapply(age1:120,function(a) annuity(a, m=m1, b=b1))
                    ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    annuities2<- c(sapply(age2:120,function(a) annuity(a, m=m2, b=b2))
                   ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    
    #compute first benefit
    B1 <- benefit1
    B2 <- benefit2
    
    y <- B2/B1
    
    #create adjustment vector
    lenDelta <- max(length(survProb2),length(survProb1))
    delta <- rep(1,lenDelta)
    DeltaStar <- rep(1,lenDelta)
    
    for (t in 1:lenDelta) {
      death1 <- surv1[t]-surv1[t+1]
      death2 <- surv2[t]-surv2[t+1]
      
      #stop when no one left in pool
      if((surv2[t+1]+surv1[t+1])==0){
        #joint stability (exact SIP measure)
        stability <- ifelse(sum((cumprod(delta)<(1-epsilonDown)))>0,
                            min(which(cumprod(delta)<(1-epsilonDown)))-1,
                            t-1) 
        
        #set all delta at 0 when everyone is dead
        DeltaStar[t:lenDelta] <- 0
        #all years the cumulative delta is stable (for approx measure)
        stableDeltaTild <- which(DeltaStar>=(1-epsilonDown))
        break
      }
      
      #compute adjustment
      delta[t] <- ((surv1[t]*survProb1[t]*annuities1[t+1]
                    +surv2[t]*y*survProb2[t]*annuities2[t+1])
                   /(surv1[t+1]*annuities1[t+1]
                     +surv2[t+1]*y*annuities2[t+1]))
      
      DeltaStar[t] <- ((surv1[1]*prod(survProb1[1:t])*annuities1[t+1]
                              +surv2[1]*y*prod(survProb2[1:t])*annuities2[t+1])
                             /(surv1[t+1]*annuities1[t+1]
                               +surv2[t+1]*y*annuities2[t+1]))
    }
    list(stability = stability,stableDeltaTild=stableDeltaTild)
  }
  
  
  # CDF of \Delta^* at 'timeK' for 2 groups evaluated at (1-'epsilon') 
  # OUTPUT: Pr[\Delta^star < (1-\epsilon)]
  pDeltaStar <- function(timeK, nb1=100,nb2=0,age1=65,age2=65,
                     ben1=1000, ben2=1000, epsilon = .1,
                     m1 = 85, b1 = 10,
                     m2 = 85, b2 = 10){
    kp1 <- tpx(timeK, age1, m=m1, b=b1)
    if(nb2==0){
      return(pbinom(nb1*kp1/(1-epsilon), nb1, kp1, lower.tail = F))
    }
    
    rfrate <- .02
    kp2 <- tpx(timeK, age2, m=m2, b=b2)
    ak1 <- as.vector(annuity(age1,rfrate, m=m1, b=b1))
    ak2 <- as.vector(annuity(age2,rfrate, m=m2, b=b2))
    y <- ben2/ben1
    
    c <- (nb1*kp1*ak1+y*nb2*kp2*ak2)/(1-epsilon)
    
    nk1 <- 0:nb1
    
    biggerCst <- (c-nk1*ak1)/(y*ak2)
    probC1 <- dbinom(nk1,nb1,kp1)
    probC2 <- pbinom(biggerCst, nb2, kp2, lower.tail = F)
    
    return(sum(probC1*probC2))
  }
  
  
  # Compute the estimates linear correlation in the CV estimator for every time
  # INPUT: N simulations of the list output of 'Stability2Gr2Variate'
  # OUTPUT: \hat{\beta}_n for n = {1,..., max n simulated}
  betaHatFn <- function(simuls){
    nsimul <- ncol(simuls)
    data <- unlist(simuls[1,])
    valueToCompare <- 1:max(data)
    xBar <- sapply(valueToCompare, function(i) sum(data>=i)/nsimul) ## good one
    yBar <- table(unlist(simuls[2,]))/nsimul
    
    betaHat <- numeric(length(valueToCompare))
    for (i in valueToCompare) {
      numerator <- 0
      denom <- 0
      for (j in 1:nsimul) {
        dataCheck <- simuls[,j]
        x <- sum(1:dataCheck$stability==i)
        y <- sum(dataCheck$stableAlphaTild==i)
        numerator <- numerator+(x-xBar[i])*(y-yBar[i])
        denom <- denom+(y-yBar[i])^2
      }
      betaHat[i] <- numerator/denom
    }
    betaHat[is.na(betaHat)] <- 1 #NA comes from those that are constant
    betaHat
  }
  
  # Control variate estimates of probabilities of stability P_n(\epsilon1,\infinity)
  # INPUT: N simulations of the list output of 'Stability2Gr2Variate', the vector
  # of \hat{\beta} form 'betaHatFn', and the 1-CDF of \Delta^* from 1-'pDeltaStar'.
  # OUTPUT: Vector of CV estimates for the stability probabilities P_n(\epsilon1,\infinity)
  CVestimates <- function(simuls,betaHat, survivalDeltaStar){
    controlled <- numeric(length(betaHat))
    nsimul <- ncol(simuls)
    for (i in seq_along(betaHat)) {
      sumValue <- 0
      for (j in 1:nsimul) {
        dataCheck <- simuls[,j]
        x <- sum(dataCheck$stability>=i)
        y <- sum(dataCheck$stableAlphaTild==i)
        sumValue <- sumValue + (x-betaHat[i]*y)
      }
      controlled[i] <- betaHat[i]*survivalDeltaStar[i]+sumValue/nsimul
    }
    controlled
  }
  
  # Give the SIP('epsilon', infinity, beta) using the the CV estimates with 
  # 'nbSimul' simulations and using a linear interpolation between the stability 
  # probabilities. Does it for heterogeneous pool with 2 groups or homogeneous 
  # pool if nb2=0 (still need to put arbitrary value for age2 and benefit2)
  SIP2Pop_CV <- function(nb1,nb2,age1,age2,benefit1,benefit2,
                         epsilon=.1, beta = .95, nbSimul = 10000,
                         m1 = 85, b1 = 10,
                         m2 = 85, b2 = 10){
    simuls <- replicate(nbSimul, Stability2Gr2Variate(nb1,nb2,age1,age2,
                                                           benefit1,benefit2,
                                                           m1 = m1, b1 = b1,
                                                           m2 = m2, b2 = b2))
    betaHat <- betaHatFn(simuls)
    survivalDeltaStar <- sapply(1:length(betaHat),
                             function(x) 1-pDeltaStar(x, nb1,nb2,age1,age2,
                                                      benefit1,benefit2,
                                                      m1 = m1, b1 = b1,
                                                      m2 = m2, b2 = b2))
    pnControlled <- CVestimates(simuls, betaHat, survivalDeltaStar)
    
    #take smoothed quantile
    p <- 1-beta
    y <- seq_along(pnControlled)
    z <- 1-pnControlled
    i <- min(which(z>=p))
    
    pStar <- (p - z[i-1])/(z[i]-z[i-1])
    SIP <- (y[i] - y[i-1])*pStar+ y[i-1]
    SIP
  }
}


####### Example Simulation #####
# compute the SIP of 2 groups of 100 members both age 65, as y goes from 0.1 to 10
# (using base parameter of function for everything else)

nbSimul <- 1000 # for fast computation but very noisy
nb1 <- 100
nb2 <- 100
age1 <- 65
age2 <- 65

y <- exp(seq(log(1/10),log(10), length.out = 101))

SIP <- sapply(y, function(x) SIP2Pop_CV(nb1,nb2,age1, age2, 1,x,
                                        nbSimul = nbSimul))

# would then apply some sort of GAM smoothing on SIP vector
plot(log(y), SIP, type= "l")
