tpx <- function(t,x,A = .0001, B = .0001, c =1.089){
  exp(-A*t-(B/log(c))*c^x*(c^t-1))
}
annuity <- function(age, rate){
  years <- 0:(150-age) #probability after age 150 usually negligible
  exp(-years*rate)%*%tpx(years, age)
}
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
GainSimulation2Groups <- function(nbPoor, agePoor=65, assetPoor=100000,
                                  nbRich=0, ageRich=65, assetRich=1000000,rate = .02){
  #simul death
  survProbPoor <- c(tpx(1,agePoor:120),0) #no one survive to age 121
  NtPoor <- c(nbPoor,numeric(length(survProbPoor)))
  for (t in 1:length(survProbPoor)) {
    NtPoor[t+1] <- rbinom(1,NtPoor[t],survProbPoor[t])
  }
  
  survProbRich <- c(tpx(1,ageRich:120),0)
  NtRich <- c(nbRich,numeric(length(survProbRich)))
  for (t in 1:length(survProbRich)) {
    NtRich[t+1] <- rbinom(1,NtRich[t],survProbRich[t])
  }
  
  #preCompute the annuities
  annuitiesPoor <- sapply(agePoor:120,function(a) annuity(a,rate))
  annuitiesRich <- sapply(ageRich:120,function(a) annuity(a,rate))
  
  #compute first benefit
  B0Poor <- assetPoor/annuitiesPoor[1]
  BPoor <- B0Poor
  
  B0Rich <- assetRich/annuitiesRich[1]
  BRich <- B0Rich
  
  #create adjustment vector for those alive
  alpha <- rep(0,max(length(survProbPoor),length(survProbRich)))
  
  for (t in 1:length(survProbPoor)) {
    if((NtRich[t+1]+NtPoor[t+1])==0){#leave alpha at 0 if no survivor
      break
    }
    
    #compute adjustment
    alpha[t] <- ((NtPoor[t]*BPoor*survProbPoor[t]*annuitiesPoor[t+1]
                  +NtRich[t]*BRich*survProbRich[t]*annuitiesRich[t+1])
                 /(NtPoor[t+1]*BPoor*annuitiesPoor[t+1]
                   +NtRich[t+1]*BRich*annuitiesRich[t+1]))
    
    BPoor <- alpha[t]*BPoor
    BRich <- alpha[t]*BRich
  }
  
  alpha <- c(1,alpha)
  assetLeftPoor <- BPoor*annuitiesPoor[t]-BPoor
  assetLeftRich <- BRich*annuitiesRich[t]-BRich
  
  #nb of death per year
  DtPoor <- -diff(NtPoor)
  DtRich <- -diff(NtRich)
  
  GiPoor <- cumsum(exp(-rate*(0:(length(alpha)-1)))*cumprod(alpha))*B0Poor
  GiPoor[t] <- GiPoor[t] + assetLeftPoor*exp(-rate*(t-1)) #paid end of year of death (i.e., time t) but invested at risk free rate, so last year of discount cancel with investment rate so we only discount t-1 year
  GiPoor <- GiPoor[seq_along(DtPoor)]
  
  GiRich <- cumsum(exp(-rate*(0:(length(alpha)-1)))*cumprod(alpha))*B0Rich
  GiRich[t] <- GiRich[t] + assetLeftRich*exp(-rate*(t-1))#paid end of year of death (i.e., time t) but invested at risk free rate, so last year of discount cancel with investment rate so we only discount t-1 year
  GiRich <- GiRich[seq_along(DtRich)]
  
  
  GiAll <- c(rep(GiPoor,DtPoor) - assetPoor,
             rep(GiRich,DtRich)- assetRich)
  deathTime <- c(rep(1:length(DtPoor),DtPoor),rep(1:length(DtRich),DtRich))
  groupPoor <- c(rep(1,nbPoor), rep(0,nbRich))
  cbind(GiAll,deathTime,groupPoor)
}



##construct scenario to be shown for Poor
scenarioMat <- matrix(c(0,65,1000000),
                      dimnames=list(NULL, c("nb","age","asset")),ncol = 3, byrow = T)
ageRich <- c(60,62,65,68,70)
assetRich <- c( 200000, 500000,1000000)
nbRich <- c(10,50,100,200)
for(i in ageRich){
  for (j in assetRich) {
    for (k in nbRich) {
      scenarioMat <- rbind(scenarioMat, c(k,i,j))
    }
  }
}



nbSimul <- 1000
set.seed(1)
sdDist <- as.vector(sqrt(asymptoVariance(65, innitialAsset = 100000)))
meanDist <- 0
for (i in 1:length(scenarioMat[,1])) {
  Gains2Groups<- replicate(nbSimul,
                           GainSimulation2Groups(nbPoor = 100,
                                                 nbRich = scenarioMat[i,1],
                                                 ageRich = scenarioMat[i,2],
                                                 assetRich = scenarioMat[i,3])) #Heterogeneity
  if(scenarioMat[i,3]==100000){
    poorGains2 <- Gains2Groups[,-3,]
  }else{
    poorGains2 <- Gains2Groups[Gains2Groups[,3,1]==1,-3,]
  }
  
  sdDist <- c(sdDist,sd(poorGains2[,1,]))
  meanDist <- c(meanDist,mean(poorGains2[,1,]))
}
plot(sdDist, meanDist/100000, ylab = "% Excess EPV", xlab = "Sd of PV",
     main = "Different Pool Composition for Poor")
idOrder <- order(sdDist[-1])
scenarioMat[idOrder,]








## construct scenario to be shown for Rich ###############################################
scenarioMatRich <- matrix(c(0,65,100000),
                      dimnames=list(NULL, c("nb","age","asset")),ncol = 3, byrow = T)
agePoor <- c(60,62,65,68,70)
assetPoor <- c(100000, 200000, 500000)
nbPoor <- c(10,50,100,200)
for(i in agePoor){
  for (j in assetPoor) {
    for (k in nbPoor) {
      scenarioMat <- rbind(scenarioMat, c(k,i,j))
    }
  }
}



nbSimul <- 1000
set.seed(1)
sdDistRich <- as.vector(sqrt(asymptoVariance(65, innitialAsset = 1000000)))
meanDistRich <- 0
for (i in 1:length(scenarioMatRich[,1])) {
  Gains2Groups<- replicate(nbSimul,
                           GainSimulation2Groups(nbRich = 100,
                                                 nbPoor = scenarioMatRich[i,1],
                                                 agePoor = scenarioMatRich[i,2],
                                                 assetPoor = scenarioMatRich[i,3])) #Heterogeneity

  poorGains2 <- Gains2Groups[Gains2Groups[,3,1]==0,-3,]
  
  
  sdDistRich <- c(sdDistRich,sd(poorGains2[,1,]))
  meanDistRich <- c(meanDistRich,mean(poorGains2[,1,]))
}
plot(sdDistRich, meanDistRich/1000000, ylab = "% Excess EPV", xlab = "Sd of PV",
     main = "Different Pool Composition for Rich")

idOrder <- order(sdDistRich[-1])
scenarioMatRich[idOrder,]
