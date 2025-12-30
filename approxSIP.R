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
  
  # Give the approx SIP('epsilon', infinity, beta) using the approx from section 
  # 4 of paper and using a linear interpolation between the approximated stability 
  # probabilities. Does it for heterogeneous pool with 2 groups or homogeneous 
  # pool if nb2=0 (still need to put arbitrary value for age2 and benefit2)
  approxSIP2Pop <- function(nb1=100,nb2=0,age1=65,age2=65, ben1=1000, ben2=1000,
                            epsilon = .1, beta = .95,
                            m1 = 85, b1 = 10, m2 = 85, b2 = 10){
    prob <- 0
    for (k in seq(1,50,by= 1)) {
      prevProb <- prob
      prob <- pDeltaStar(k, nb1,nb2,age1,age2, ben1, ben2, epsilon, m1, b1, m2, b2)
      if(prob>1-beta){
        return((k-1)*(prob-(1-beta))/(prob-prevProb)
               +k*(1-(prob-(1-beta))/(prob-prevProb)))
      }
    }
  }
}


####### Example Simulation #####
# This is an example to create a surface of SIP for different age and benefit 
# level for group 2

## fix parameters
nb2 <- 100
nb1 <- 100 
age1 <- 65 
b1=10
b2=10
benefit1 <- 1

# axis parameters
age2 <- seq(55, 75, by = .1)
benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
benefit2 <- benefit1*benefitMultiplier

# compute the points
SIP_Approx <- matrix(0,length(age2),length(benefitMultiplier))
for (i in seq_along(age2)) {
  for (j in seq_along(benefit2)) {
    SIP_Approx[i,j] <- approxSIP2Pop(nb1, nb2, age1, age2[i], 
                                     benefit1,benefit2[j],
                                     b1 = b1, b2=b2)
  }
}


# 3D Plot
{
  library(plotly)
  winterColormap <- c(rgb(0,(0:256)/256,(1-((0:256)/256))*.5+.5))
  
  p <- plot_ly()%>%add_surface(
    x = ~benefit2,
    y = ~age2,
    z = ~SIP_Approx,
    type = "surface",
    showscale=F,
    colors = winterColormap)
  # Add grid lines in x-direction
  for (j in seq(1,length(age2), length.out = 20)) {
    p <- p %>% add_trace(
      x = benefit2,
      y = rep(age2[j], length(benefit2)),
      z = SIP_Approx[j,],
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 3),
      showlegend = FALSE
    )
  }
  
  # Add grid lines in y-direction
  for (i in seq(1,length(benefit2), length.out = 10)) {
    p <- p %>% add_trace(
      x = rep(benefit2[i], length(age2)),
      y = age2,
      z = SIP_Approx[,i],
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 3),
      showlegend = FALSE
    )
  }
  p <- p %>%
    layout(
      scene = list(
        xaxis = list(title = "Benefit 2",
                     type = "log",
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE),
        yaxis = list(title = "Age 2",
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE),
        zaxis = list(title = "approx SIP",
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE),
        aspectratio = list(x = 1, y = 2, z = 1)
      ),
      margin=list(t=0,l=0,b=0,r=0)
    )
  p
}
