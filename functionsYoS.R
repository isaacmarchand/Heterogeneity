#######Create Base functions
{
  # Surival prob function with Markham law
  tpx <- function(t,x,A = .0001, B = .0001, c =1.089){ #base parameter female
    exp(-A*t-(B/log(c))*c^x*(c^t-1))
  }
  
  # AnnuityFunction
  annuity <- function(age, rate, sexe = "F"){
    years <- 0:(150-age) #probability after age 150 usually negligible
    v <- exp(-years*rate)
    if(sexe == "F"){
      v%*%tpx(years, age)
    }else{
      v%*%tpx(years, age, A_male, B_male)
    }
  }
  
  # Function that simulate the nb of years of stability for 1 scenario for
  StabilityCalc2Pop <- function(nbPoor, nbRich, agePoor, ageRich, assetPoor,
                                assetRich, epsilonDown = .1){
    rate <- .02 #fix interest
    #simul death
    survProbPoor <- c(tpx(1,agePoor:120),0) #no one survive to age 121
    survProbRich <- c(tpx(1,ageRich:120),0) #no one survive to age 121
    survPoor <- c(nbPoor,numeric(length(survProbPoor)))
    for (t in 1:length(survProbPoor)) {
      survPoor[t+1] <- rbinom(1,survPoor[t],survProbPoor[t])
    }
    survRich <- c(nbRich,numeric(length(survProbRich)))
    for (t in 1:length(survProbRich)) {
      survRich[t+1] <- rbinom(1,survRich[t],survProbRich[t])
    }
    
    #preCompute the annuities
    annuitiesPoor <- sapply(agePoor:120,function(a) annuity(a,rate))
    annuitiesRich<- sapply(ageRich:120,function(a) annuity(a,rate))
    
    #compute first benefit
    BRich <- assetRich/annuitiesRich[1]
    BPoor <- assetPoor/annuitiesPoor[1]
    
    #create adjustment vector
    alpha <- rep(1,max(length(survProbRich),length(survProbPoor)))
    
    for (t in 1:max(length(survProbRich),length(survProbPoor))) {
      deathPoor <- survPoor[t]-survPoor[t+1]
      deathRich <- survRich[t]-survRich[t+1]
      
      #if no one left in pool set stability as max value
      if((survRich[t+1]+survPoor[t+1])==0){
        stability <- t-1
        break
      }
      
      #compute adjustment
      alpha[t] <- ((survPoor[t]*BPoor*survProbPoor[t]*annuitiesPoor[t+1]
                    +survRich[t]*BRich*survProbRich[t]*annuitiesRich[t+1])
                   /(survPoor[t+1]*BPoor*annuitiesPoor[t+1]
                     +survRich[t+1]*BRich*annuitiesRich[t+1]))
      
      #check if benefits are still stable
      if(prod(alpha)<(1-epsilonDown)){
        stability <- t-1
        break
      }
      
      BPoor <- alpha[t]*BPoor
      BRich <- alpha[t]*BRich
    }
    stability
  }
  
  
  # Extract smoothed quantile(1-beta in stability definition) to 
  # get continuous stability
  smoothedQuantiles <- function(x,p){
    x <- sort(x)
    y <- unique(x)
    Ly <- length(y)
    z <- sapply(1:Ly,function(i) sum(y[i]==x))/length(x)
    zModified <- c(z[1]/2,(z[2:Ly]+z[1:(Ly-1)])/2, z[Ly]/2)
    cumul <- cumsum(zModified)
    i <- min(which(cumul>=p))
    pStar <- (p - (cumul[i] - zModified[i]))/zModified[i]
    if(i == 1){
      VaR <- y[1]
    }else if (i==(Ly+1)) {
      VaR <- y[Ly]
    }else{
      VaR <- (y[i] - y[i-1])*pStar+ y[i-1]
    }
    VaR
  }
  
}









######Example Simulation

#nb of simulations from which we will extract the 1-beta quatile from
nbSimul <- 1000 ##Low nb of simulation to keep example quick to run

# Stability Definition
beta <- .95
epsilonDown <- .1

#different scenarios (100poor+0Rich, 200poor+0Rich, 100poor+100Rich)
nbPoor <- c(100,200,100)
nbRich <- c(0,0,100)
##Other scenario details
agePoor <- 65
assetPoor <- 100000
assetRich <- 1000000
ageRich <- 65



# Comput Years of Stability (YoS)

## Compute simulation of when the number of year after which the epsilon 
## threshold is busted for 3 distinct scenario
Stability <- sapply(1:3,
                    function(x) replicate(nbSimul,
                                          StabilityCalc2Pop(nbPoor[x]
                                                            ,nbRich[x]
                                                            ,agePoor
                                                            ,ageRich
                                                            ,assetPoor
                                                            ,assetRich
                                                            ,epsilonDown = epsilonDown)))
## Compute a smoother version of the VaR (1-beta) to get the simulated number 
## of years of stability for the 3 distinct scenario
riskStability <- sapply(1:length(nbPoor),
                        function(j) smoothedQuantiles(Stability[,j],1-beta))
riskStability


