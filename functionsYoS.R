#######Create Base functions #####
{
  # Surival prob function with Markham law
  tpx <- function(t,x,A = .0001, B = .0001, c =1.089){ #base parameter female
    exp(-A*t-(B/log(c))*c^x*(c^t-1))
  }
  
  # AnnuityFunction
  annuity <- function(age, rate=.02, sexe = "F"){
    years <- 0:(150-age) #probability after age 150 usually negligible
    v <- exp(-years*rate)
    if(sexe == "F"){
      v%*%tpx(years, age)
    }else{
      v%*%tpx(years, age, A_male, B_male)
    }
  }
  
  # Function that simulate the nb of years of stability for 1 scenario for
  StabilityCalc2Pop <- function(nb1, nb2, age1, age2, asset1,
                                asset2, epsilonDown = .1){
    rate <- .02 #fix interest
    #simul death
    survProb1 <- c(tpx(1,age1:120),0) #no one survive to age 121
    survProb2 <- c(tpx(1,age2:120),0) #no one survive to age 121
    surv1 <- c(nb1,numeric(length(survProb1)))
    for (t in 1:length(survProb1)) {
      surv1[t+1] <- rbinom(1,surv1[t],survProb1[t])
    }
    surv2 <- c(nb2,numeric(length(survProb2)))
    for (t in 1:length(survProb2)) {
      surv2[t+1] <- rbinom(1,surv2[t],survProb2[t])
    }
    
    #preCompute the annuities
    annuities1 <- sapply(age1:120,function(a) annuity(a,rate))
    annuities2<- sapply(age2:120,function(a) annuity(a,rate))
    
    #compute first benefit
    B2 <- asset2/annuities2[1]
    B1 <- asset1/annuities1[1]
    
    #create adjustment vector
    lenAlpha <- max(length(survProb2),length(survProb1))
    alpha <- rep(1,lenAlpha)
    
    for (t in 1:lenAlpha) {
      death1 <- surv1[t]-surv1[t+1]
      death2 <- surv2[t]-surv2[t+1]
      
      #if no one left in pool set stability as max value
      if((surv2[t+1]+surv1[t+1])==0){
        stability <- t-1
        break
      }
      
      #compute adjustment
      alpha[t] <- ((surv1[t]*B1*survProb1[t]*annuities1[t+1]
                    +surv2[t]*B2*survProb2[t]*annuities2[t+1])
                   /(surv1[t+1]*B1*annuities1[t+1]
                     +surv2[t+1]*B2*annuities2[t+1]))
      
      #check if benefits are still stable
      if(prod(alpha)<(1-epsilonDown)){
        stability <- t-1
        break
      }
      
      B1 <- alpha[t]*B1
      B2 <- alpha[t]*B2
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
  
  # Compute probability that cumulative adjustment is under treshold at specific time
  thetaK <- function(timeK, nb1=100,nb2=0,age1=65,age2=65, ben1=1000, ben2=1000,
                     epsilon = .1){
    kp1 <- tpx(timeK, age1)
    if(nb2==0){
      return(pbinom(nb1*kp1/(1-epsilon), nb1, kp1, lower.tail = F))
    }
    
    rfrate <- .02
    kp2 <- tpx(timeK, age2)
    ak1 <- as.vector(annuity(age1,rfrate))
    ak2 <- as.vector(annuity(age2,rfrate))
    y <- ben2/ben1
    
    c <- (nb1*kp1*ak1+y*nb2*kp2*ak2)/(1-epsilon)
    
    nk1 <- 0:nb1
    
    biggerCst <- (c-nk1*ak1)/(y*ak2)
    probC1 <- dbinom(nk1,nb1,kp1)
    probC2 <- pbinom(biggerCst, nb2, kp2, lower.tail = F)
    
    return(sum(probC1*probC2))
  }
  
  #compute theApprox for YoS presented in section 11.2
  compApproxYoS <- function(nb1=100,nb2=0,age1=65,age2=65, ben1=1000, ben2=1000,
                            epsilon = .1, beta = .95){
    prob <- 0
    for (t in 1:50) {
      prevProb <- prob
      prob <- thetaK(t, nb1,nb2,age1,age2, ben1, ben2, epsilon)
      if(prob>1-beta){
        return((t-1)*(prob-(1-beta))/(prob-prevProb)
               +t*(1-(prob-(1-beta))/(prob-prevProb)))
        # return((t-1)) 
      }
    }
  }
  
  
}









######Example Simulation #####

#nb of simulations from which we will extract the 1-beta quatile from
nbSimul <- 1000 ##Low nb of simulation to keep example quick to run

# Stability Definition
beta <- .95
epsilonDown <- .1

#different scenarios (100_C1+0_C2, 200_C1+0_C2, 100_C1+100_C2)
nb1 <- 100
nb2 <- 100
##Other scenario details
age1 <- 65
benefit1 <- 1000
benefit2 <- 10000
age2 <- c(60,65,70)



# Compute Years of Stability (YoS)

## Compute simulation of when the number of year after which the epsilon 
## threshold is busted for 3 distinct scenario
Stability <- sapply(1:3,
                    function(x) replicate(nbSimul,
                                          StabilityCalc2Pop(nb1
                                                            ,nb2
                                                            ,age1
                                                            ,age2[x]
                                                            ,benefit1*annuity(age1)
                                                            ,benefit2*annuity(age2[x])
                                                            ,epsilonDown = epsilonDown)))
## Compute a smoother version of the VaR (1-beta) to get the simulated number 
## of years of stability for the 3 distinct scenario
riskStability <- sapply(1:length(age2),
                        function(j) smoothedQuantiles(Stability[,j],1-beta))
riskStability





#Compute Approx for YoS
riskStabilityApprox <- sapply(1:3,
                        function(x) compApproxYoS(nb1,nb2,age1,age2[x]
                                                  ,benefit1,benefit2
                                                  ,epsilon = epsilonDown))
riskStabilityApprox



######Example Contour plot with Approximation #####

####Compute riskStabilityApprox (5 to 20 second to run)
{
  nb1 <- 100
  nb2 <- 100
  age1 <- 65
  benefit1 <- 1000
  
  age2 <- seq(55, 75, by = .1)
  benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
  benefit2 <- benefit1*benefitMultiplier
  
  riskStabilityApprox <- matrix(0,length(age2),length(benefitMultiplier))
  for (i in seq_along(age2)) {
    for (j in seq_along(benefit2)) {
      riskStabilityApprox[i,j] <- compApproxYoS(nb1, nb2, age1, age2[i], benefit1,
                                                benefit2[j], epsilonDown, .95)
    }
  }
}

# get better (green) and worse (red) areas
{
  allAges <- age2
  
  homoApproxYoS <- compApproxYoS(nb1+nb2, age1 = age1)
  riskSmallHomoApprox <- compApproxYoS(nb1, age1 = age1)
  
  dfbetterApproxYoS <- list(YoS = riskStabilityApprox[riskStabilityApprox>=homoApproxYoS])
  dfbetterApproxYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                    ncol = length(benefitMultiplier))[riskStabilityApprox
                                                                      >=homoApproxYoS]
  dfbetterApproxYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                                nrow = length(age2),
                                                byrow = T)[riskStabilityApprox>=homoApproxYoS]
  
  
  dfworsteApproxYoS <- list(YoS = riskStabilityApprox[riskStabilityApprox<=riskSmallHomoApprox])
  dfworsteApproxYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                    ncol = length(benefitMultiplier))[riskStabilityApprox
                                                                      <=riskSmallHomoApprox]
  dfworsteApproxYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                                nrow = length(age2),
                                                byrow = T)[riskStabilityApprox<=riskSmallHomoApprox]
  
}

# contour Plot
{
  # Prepare data in long format
  df <- melt(riskStabilityApprox)
  colnames(df) <- c("ageIndex", "benefitIndex", "YoS")
  df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
  df$age2 <- age2[df$ageIndex]
  
  # Create contour plot
  p <- plot_ly(
    data = df,
    x = ~age2,
    y = ~benefitMultiplier,
    z = ~YoS,
    type = "contour",
    showscale = FALSE,
    contours = list(
      coloring = "lines",  # or "lines", "none"
      showlabels = TRUE
    ),
    line = list(smoothing = 0),
    colorscale = list(c(0, "black"), c(1, "black")),
    reversescale = FALSE
  ) %>%
    layout(
      plot_bgcolor = "lightgrey",   # uniform background color
      paper_bgcolor = "white",  # outside background
      title = paste("Stability with Size Cohort 2:",nb2,"people"),
      xaxis = list(title = "Age Cohort 2",
                   showgrid = FALSE,
                   range = c(min(df$age2), max(df$age2))
      ),
      yaxis = list(title = "Benefit: Cohort 2/Cohort 1",
                   type = "log",
                   showgrid = FALSE,
                   range= log(c(min(df$benefitMultiplier), max(df$benefitMultiplier)),10)
      ),
      showlegend = F
    )%>%
    add_trace(
      data = dfbetterApproxYoS,
      x = ~age2,
      y = ~benefitMultiplier,
      type = "scatter",
      mode = "markers",
      marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
      name = "Preferred Region",
      inherit = FALSE
    )%>%
    add_trace(
      data = dfworsteApproxYoS,
      x = ~age2,
      y = ~benefitMultiplier,
      type = "scatter",
      mode = "markers",
      marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
      name = "No No Region",
      inherit = FALSE
    )
  
  p
}


