# Surival prob function with Markham law
tpx <- function(t,x,A = .0001, B = .0001, c =1.089){
  exp(-A*t-(B/log(c))*c^x*(c^t-1))
}

# AnnuityFunction
annuity <- function(age, rate){
  years <- 0:(150-age) #probability after age 150 usually negligible
  exp(-years*rate)%*%tpx(years, age)
}




# function to simulate PV of benefits from Homogeneous pool
GainSimulationHomo <- function(nb, age=65, asset=100000, rate = .02){
  #simul death
  survProb <- c(tpx(1,age:120),0) #no one survive to age 121
  Nt <- c(nb,numeric(length(survProb)))
  for (t in 1:length(survProb)) {
    Nt[t+1] <- rbinom(1,Nt[t],survProb[t])
  }
  
  #preCompute the annuities
  annuities <- sapply(age:120,function(a) annuity(a,rate))
  
  #compute first benefit
  B0 <- asset/annuities[1]
  B <- B0
  
  #create adjustment vector for those alive
  alpha <- rep(0,length(survProb))
  
  for (t in 1:length(survProb)) {
    if(Nt[t+1]==0){#leave alpha at 0 if no survivor
      break
    }
    
    #compute adjustment
    alpha[t] <- ((Nt[t]*B*survProb[t]*annuities[t+1])
                 /(Nt[t+1]*B*annuities[t+1]))
    
    B <- alpha[t]*B
  }
  alpha <- c(1,alpha)
  assetLeft <- B*annuities[t]-B
  
  Gi <- cumsum(exp(-rate*(0:(length(alpha)-1)))*cumprod(alpha))*B0
  Gi[t] <- Gi[t] + assetLeft*exp(-rate*(t-1)) #paid end of year of death (i.e., time t) but invested at risk free rate, so last year of discount cancel with investment rate so we only discount t-1 year
  
  #nb of death per year
  Dt <- -diff(Nt)
  
  GiAll <- rep(Gi[-length(Gi)],Dt)
  deathTime <- rep(1:length(Dt),Dt)
  cbind(GiAll - asset,deathTime)
}




## Example
nbSimul <- 1000
poolSize <- c(10,100,1000) #simulate 3 homogeneous pool of different size
set.seed(1)
GainsHomo1 <- replicate(nbSimul*100,GainSimulationHomo(poolSize[1]))
GainsHomo2 <- replicate(nbSimul*10,GainSimulationHomo(poolSize[2]))
GainsHomo3 <- replicate(nbSimul,GainSimulationHomo(poolSize[3]))

### plot full distribution
{
  plot(density(as.vector(GainsHomo1[,1,])),col = 1, ylim = c(0,0.00001),
       xlim = c(-100000,200000),
       main = "PDF of the Gain for homogeneous pools", xlab = "Gains")
  points(density(as.vector(GainsHomo2[,1,])), type = "l",col = 2)
  points(density(as.vector(GainsHomo3[,1,])), type = "l",col = 3)
  legend("topright", legend = c("10 people","100 people","1000 people"),
         col = c(1,2,3), lty = c(1,1,1))
}
### plot distribution of gain conditional on age of death
{
  distPerAge <- sapply(1:(120-65)
                       ,function(i) ((GainsHomo3[,1,])[which(GainsHomo3[,2,]==i)]))
  plot(density(unlist(distPerAge[6])),col = 1,xlim = c(-100000,200000),
       main = "PDF Conditional on age of death", xlab = "Gains")
  points(density(unlist(distPerAge[11])),col = 2,type = "l")
  points(density(unlist(distPerAge[21])),col = 3,type = "l")
  points(density(unlist(distPerAge[31])),col = 4,type = "l")
  
  distPerAge <- sapply(1:(120-65)
                       ,function(i) ((GainsHomo2[,1,])[which(GainsHomo2[,2,]==i)]))
  points(density(unlist(distPerAge[6])),col = 1,type = "l", lty =2)
  points(density(unlist(distPerAge[11])),col = 2,type = "l", lty =2)
  points(density(unlist(distPerAge[21])),col = 3,type = "l", lty =2)
  points(density(unlist(distPerAge[31])),col = 4,type = "l", lty =2)
  legend("top", legend = c("1000 people pool","100 people pool"),, lty = c(1,2))
  legend("topright", legend = c("70 y/o","75 y/o","85 y/o","95 y/o"),
         col = c(1,2,3,4), lty = c(1,1,1,1))
}








# function to simulate PV of benefits from Heterogeneous pool with 2 subGroup
GainSimulation2Groups <- function(nbPoor, agePoor=65, assetPoor=100000,
                                  nbRich=0, ageRich=65, assetRich=1000000,
                                  rate = .02){
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






## Example 3 Scenario with poor 65y and rich 70y
nbSimul <- 10000
set.seed(1)
Gains2Groups100_0 <- replicate(nbSimul,GainSimulation2Groups(nbPoor = 100,
                                                             nbRich = 0,
                                                             ageRich = 70))
Gains2Groups200_0 <- replicate(nbSimul,GainSimulation2Groups(nbPoor = 100,
                                                             nbRich = 100,
                                                             ageRich = 70))
Gains2Groups100_100 <- replicate(nbSimul,GainSimulation2Groups(nbPoor = 100,
                                                               nbRich = 200,
                                                               ageRich = 70))


### extract result only for poor people (col groupPoor == 1)
poorGains2Groups100_0 <- Gains2Groups100_0[Gains2Groups100_0[,3,1]==1,-3,]
poorGains2Groups200_0 <- Gains2Groups200_0[Gains2Groups200_0[,3,1]==1,-3,]
poorGains2Groups100_100 <- Gains2Groups100_100[Gains2Groups100_100[,3,1]==1,-3,]


### plot full dist from poors perspective
{
  plot(density(as.vector(poorGains2Groups100_0[,1,])),col = 1, ylim = c(0,0.00001),xlim = c(-100000,200000),
       main = "PDF of the Gain for Heterogeneous pools", xlab = "Gains")
  points(density(as.vector(poorGains2Groups200_0[,1,])), type = "l", col = 2)
  points(density(as.vector(poorGains2Groups100_100[,1,])), type = "l", col = 4)
  legend("topright", legend = c("100 Homogeneous people","200 Homogeneous people","100 poor people with 100 Rich"),
         col = c(1,2,4), lty = c(1,1,1))
}

### plot conditional dist from poors perspective
{
  distPerAgeHomo100 <- sapply(1:(120-65),function(i) ((poorGains2Groups100_0[,1,])[which(poorGains2Groups100_0[,2,]==i)]))
  plot(density(unlist(distPerAgeHomo100[6])),col = 1,xlim = c(-100000,200000),
       main = "PDF Conditional on age of death for poor people in a pool", xlab = "Gains")
  points(density(unlist(distPerAgeHomo100[11])),col = 2,type = "l")
  points(density(unlist(distPerAgeHomo100[21])),col = 3,type = "l")
  points(density(unlist(distPerAgeHomo100[31])),col = 4,type = "l")
  points(density(unlist(distPerAgeHomo100[36])),col = 6,type = "l")
  
  distPerAgeHomo200 <- sapply(1:(120-65),function(i) ((poorGains2Groups200_0[,1,])[which(poorGains2Groups200_0[,2,]==i)]))
  points(density(unlist(distPerAgeHomo200[6])),col = 1,type = "l", lty =2)
  points(density(unlist(distPerAgeHomo200[11])),col = 2,type = "l", lty =2)
  points(density(unlist(distPerAgeHomo200[21])),col = 3,type = "l", lty =2)
  points(density(unlist(distPerAgeHomo200[31])),col = 4,type = "l", lty =2)
  points(density(unlist(distPerAgeHomo200[36])),col = 6,type = "l", lty =2)
  
  distPerAge100_100 <- sapply(1:(120-65),function(i) ((poorGains2Groups100_100[,1,])[which(poorGains2Groups100_100[,2,]==i)]))
  points(density(unlist(distPerAge100_100[6])),col = 1,type = "l", lty =3)
  points(density(unlist(distPerAge100_100[11])),col = 2,type = "l", lty =3)
  points(density(unlist(distPerAge100_100[21])),col = 3,type = "l", lty =3)
  points(density(unlist(distPerAge100_100[31])),col = 4,type = "l", lty =3)
  points(density(unlist(distPerAge100_100[36])),col = 6,type = "l", lty =3)
  
  legend("topright", legend = c("100 Homogeneous people","200 Homogeneous people",
                                "100 poor people with 100 Rich"), lty = c(1,2,3))
  legend("topleft", legend = c("70 y/o","75 y/o","85 y/o","95 y/o"),
         col = c(1,2,3,4), lty = c(1,1,1,1))
}

### Create Tables of mean ans sd of dist
{
  meanHomo <- mean(unlist(distPerAgeHomo100))
  sdHomo <- sd(unlist(distPerAgeHomo100))
  for (i in 1:length(distPerAgeHomo100)) {
    meanHomo <- c(meanHomo,mean(unlist(distPerAgeHomo100[i])))
    sdHomo <- c(sdHomo,sd(unlist(distPerAgeHomo100[i])))
  }
  meanHomo <- round(meanHomo, digits = 2)
  sdHomo <- round(sdHomo, digits = 2)
  
  meanHomo200 <- mean(unlist(distPerAgeHomo200))
  sdHomo200 <- sd(unlist(distPerAgeHomo200))
  for (i in 1:length(distPerAgeHomo200)) {
    meanHomo200 <- c(meanHomo200,mean(unlist(distPerAgeHomo200[i])))
    sdHomo200 <- c(sdHomo200,sd(unlist(distPerAgeHomo200[i])))
  }
  meanHomo200 <- round(meanHomo200, digits = 2)
  sdHomo200 <- round(sdHomo200, digits = 2)
  
  mean100_100 <- mean(unlist(distPerAge100_100))
  sd100_100 <- sd(unlist(distPerAge100_100))
  for (i in 1:length(distPerAge100_100)) {
    mean100_100 <- c(mean100_100,mean(unlist(distPerAge100_100[i])))
    sd100_100 <- c(sd100_100,sd(unlist(distPerAge100_100[i])))
  }
  mean100_100 <- round(mean100_100, digits = 2)
  sd100_100 <- round(sd100_100, digits = 2)
  
  tableMean <- rbind(c(meanHomo[1],meanHomo[7],meanHomo[12],meanHomo[22],meanHomo[31],meanHomo[41]),
                     c(meanHomo200[1],meanHomo200[7],meanHomo200[12],meanHomo200[22],meanHomo200[31],meanHomo200[41]),
                     c(mean100_100[1],mean100_100[7],mean100_100[12],mean100_100[22],mean100_100[31],mean100_100[41]))
  colnames(tableMean) <- c("total", "70 y/o","75 y/o","85 y/o","95 y/o","105 y/o")
  rownames(tableMean) <- c("100 Homogeneous people","200 Homogeneous people",
                           "100 poor people with 100 Rich")
  
  tableSd <- rbind(c(sdHomo[1],sdHomo[7],sdHomo[12],sdHomo[22],sdHomo[31],sdHomo[41]),
                   c(sdHomo200[1],sdHomo200[7],sdHomo200[12],sdHomo200[22],sdHomo200[31],sdHomo200[41]),
                   c(sd100_100[1],sd100_100[7],sd100_100[12],sd100_100[22],sd100_100[31],sd100_100[41]))
  colnames(tableSd) <- c("total", "70 y/o","75 y/o","85 y/o","95 y/o","105 y/o")
  rownames(tableSd) <- c("100 Homogeneous people","200 Homogeneous people",
                         "100 poor people with 100 Rich")
  
  print(tableMean)
  print(tableSd)
  
  plot(66:120,meanHomo[-1], type = "l", ylim = c(-100000,200000),
       xlab = "age of death", ylab = "Gain")
  points(66:120,meanHomo200[-1], type = "l", col = 2)
  points(66:120,mean100_100[-1], type = "l", col = 3)
  points(66:120,sdHomo[-1], type = "l", ylim = c(0,100000),
         xlab = "age of death", ylab = "sd of gain", lty = 2)
  points(66:120,sdHomo200[-1], type = "l", col = 2, lty = 2)
  points(66:120,sd100_100[-1], type = "l", col = 3, lty = 2)
  legend("topleft", legend = c("Mean","SD"), lty = c(1,2))
  legend("bottomright", legend = c("100 Homogeneous people","200 Homogeneous people",
                                   "100 poor people with 100 Rich"),
         col = c(1,2,3), lty = c(1,1,1))
  
}

