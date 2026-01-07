####### Packages #####
library(plotly)
library(reshape2)
library(RColorBrewer)
library(reticulate)
library(ggplot2)
library(htmlwidgets)
library(webshot)
library(mgcv)

####### collect paths #####
args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=2){
  stop("The script requires 2 trailing arguments, 1.the path where you want the figures saved 2.the path to your python software")
}

####### Export ######
exportPath <- args[1]

reticulate::use_python(args[2])
reticulate::py_config()

# ####if need to install python
# install.packages("reticulate")
# library(reticulate)
# use_python(install_python())
# py_install(c("kaleido==0.2.1", "plotly"))

####### Export design choice ####
fontType <- 'Verdana' #'Verdana' for report 'Times New Roman' for paper
axisFont <- list(size=15, family = fontType)
titleFont <- list(size=30, family = fontType)
legendFont <- list(size=12, family = fontType)

#pixel wanted for full page plot
pixelsFullWidth <- 1240*(5.5/8) #removing margins
pixelsFullHeight <- 1754*(9.3/11) #removing margins


####### Base functions #####
{
  ######Colors
  rgbSOA <- matrix(c(2,77,124,186,191,51,119,196,213,253,206,7,210,49,56,1,1,1, 255,255,255), byrow = TRUE, ncol = 3)/255
  winterColormap <- c(rgb(0,(0:256)/256,(1-((0:256)/256))*.5+.5))
  
  # Surival prob function with Markham law
  tpx <- function(t,x,lambda = 0, b = 10, m = 85){
    exp(-lambda*t-exp((x-m)/b)*(exp(t/b)-1))
  }
  
  # AnnuityFunction
  annuity <- function(age, rate, b = 10, m = 85){
    years <- 0:(150-age) #probability after age 150 usually negligible
    v <- exp(-years*rate)
    v%*%tpx(years, age, b = b, m = m)
  }
  
  # Get variance of 1st period adjustment as proxy
  SD1Periode <- function(age1=65, B01=1000, nb1=100,
                         age2=65, B02=2000, nb2=100, rate = .02,
                         m1 = 85, b1 = 10,
                         m2 = 85, b2 = 10){
    p1 <- tpx(1,age1, m=m1, b=b1)
    p2 <- tpx(1,age2, m=m2, b=b2)
    y <- B02/B01
    
    pN1_0 <- dbinom(0,nb1,p1) 
    pN2_0 <- dbinom(0,nb2,p2) 
    
    mu <- (nb1*p1+y*nb2*p2)/((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))
    sigmaSq <- ((nb1*p1*(1-p1)+nb1^2*p1^2+2*y*nb1*p1*nb2*p2+y^2*nb2*p2*(1-p2)
                 +y^2*nb2^2*p2^2)/((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))) - mu^2
    
    varTerm1 <- ((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))*(1/mu^2+(3*sigmaSq/mu^4))
    varTerm2 <- ((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))^2*(1/mu+(sigmaSq/mu^3))^2
    
    sqrt((nb1*p1+y*nb2*p2)^2*(varTerm1-varTerm2))
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
  
  #Benefit return and SIP for different asset allocation
  meanVar <- function(nb1, nb2, age1, age2, asset1, asset2,
                      assetSplit = .5){
    totAsset <-  nb1*asset1+nb2*asset2
    
    benefit1 <- as.vector(asset1/annuity(age1,.02))
    benefit2 <- as.vector(asset2/annuity(age2,.02))
    
    newA1 <- (assetSplit*totAsset)/nb1
    newA2 <- (totAsset-newA1*nb1)/nb2
    
    newBen1 <- as.vector(newA1/annuity(age1,.02))
    newBen2 <- as.vector(newA2/annuity(age2,.02))
    
    SIP1 <- approxSIP2Pop(nb1 = nb1, age1 = age1, ben1 = benefit1, nb2 = 0)
    
    SIP2 <- approxSIP2Pop(nb1 = nb2, age1 = age2, ben1 = benefit2, nb2 = 0)
    
    SIPhat <- approxSIP2Pop(nb1 = nb1, age1 = age1, ben1 = newBen1,
                            nb2 = nb2, age2 = age2, ben2 = newBen2)
    c(SIPhat, newBen1/asset1,SIPhat, newBen2/asset2)
  }
  
  #Benefit return and SIP from SIP matrix
  meanVarMatrix <- function(nb1, nb2, age1, age2, asset1, asset2, matSIP,
                            assetSplit = .5){
    totAsset <-  nb1*asset1+nb2*asset2
    
    benefit1 <- as.vector(asset1/annuity(age1,.02))
    benefit2 <- as.vector(asset2/annuity(age2,.02))
    
    newA1 <- (assetSplit*totAsset)/nb1
    newA2 <- (totAsset-newA1*nb1)/nb2
    
    newBen1 <- as.vector(newA1/annuity(age1,.02))
    newBen2 <- as.vector(newA2/annuity(age2,.02))
    newY <- newBen2/newBen1
    
    age2Vec <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
    
    SIPhat <- matSIP[age2Vec==age2, which.min(abs(benefitMultiplier-newY))]
    c(SIPhat, newBen1/asset1,SIPhat, newBen2/asset2)
  }
  
  #Quadratic Utility
  utilityFn <- function(SIP,initialYield, alpha = 2){
    initialYield-.5*alpha*SIP^-2
  }
  
  utilityCurve <- function(utility,SIP, alpha = 2){
    utility+.5*alpha*SIP^-2
  }
  
  #find asset allocation to go on porder of red zone
  nonZeroToOptim <- function(nb1, nb2, age1, age2, asset1, asset2,
                             assetSplit = .5){
    totAsset <-  nb1*asset1+nb2*asset2
    
    benefit1 <- as.vector(asset1/annuity(age1,.02))
    benefit2 <- as.vector(asset2/annuity(age2,.02))
    
    newA1 <- (assetSplit*totAsset)/nb1
    newA2 <- (totAsset-newA1*nb1)/nb2
    
    newBen1 <- as.vector(newA1/annuity(age1,.02))
    newBen2 <- as.vector(newA2/annuity(age2,.02))
    
    SIP1 <- approxSIP2Pop(nb1 = nb1, age1 = age1, ben1 = benefit1, nb2 = 0)
    
    SIP2 <- approxSIP2Pop(nb1 = nb2, age1 = age2, ben1 = benefit2, nb2 = 0)
    
    SIPhat <- approxSIP2Pop(nb1 = nb1, age1 = age1, ben1 = newBen1,
                            nb2 = nb2, age2 = age2, ben2 = newBen2)
    
    diff1 <- SIPhat-SIP1
    diff2 <- SIPhat-SIP2
    
    minDiff <- min(diff1,diff1)
    
    ifelse(minDiff>0, -1000, minDiff)
  }
}

###############################################################################
# Section 3 SIP
###############################################################################
####### Fig 1:Empirical Dist of NB of stable years with smoothing######
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height


# adjustable parameters
nb2 <- 100     #can be anything, computed directly in code
nb1 <- 100     
nbSimul <- 10000  #can add more but <5 sec to run for 10000 simul
beta <- .95       #treshhold illustrated in plot

##preparing data
{
  set.seed(3)
  {
    age1 <- 65
    benefit1 <- 1000
    age2 <- age1
    asset1 <- as.vector(benefit1 %*% annuity(age1, rate = .02))
    asset2 <- asset1
    
    Stability <- replicate(nbSimul,StabilityCalc2Pop(nb1,
                                                     nb2,
                                                     age1,
                                                     age2,
                                                     asset1,
                                                     asset2))
    
    cumulDist <- ecdf(Stability)
    discreteVaR <- max(which(cumulDist(1:30)<=(1-beta)))
    
    p <- (1-beta)
    x <- sort(Stability)
    y <- unique(x)
    Ly <- length(y)
    z <- sapply(1:Ly,function(i) sum(y[i]==x))/length(x)
    cumul <- cumsum(z)
    i <- min(which(cumul>=p))
    pStar <- (p - (cumul[i] - z[i]))/z[i]
    if(i == 1){
      VaR <- y[1]
    }else if (i==(Ly+1)) {
      VaR <- y[Ly]
    }else{
      VaR <- (y[i] - y[i-1])*pStar+ y[i-1]
    }
    VaR
    
    #specific data manipulation for plot
    {
      # Suppose cumulDist is an ecdf object
      x_jumps <- environment(cumulDist)$x   # the sample values
      y_jumps <- cumulDist(x_jumps)
      
      # Construct horizontal segments (x[i], y[i]) to (x[i+1], y[i])
      x_segments <- c()
      y_segments <- c()
      for (i in seq_along(x_jumps)) {
        if (i < length(x_jumps)) {
          x_segments <- c(x_segments, x_jumps[i], x_jumps[i+1], NA)
          y_segments <- c(y_segments, y_jumps[i], y_jumps[i], NA)
        }
      }
    }
  }
  
  {
    fig <- plot_ly()
    fig <- fig %>%
      add_trace(x = x_segments, y = 1-y_segments, type = 'scatter', mode = 'lines',
                line = list(color = "#010101", shape = "hv"),
                name = "Non-smoothed probability")%>%
      # Add points at the jumps
      add_trace(x = x_jumps, y = 1-y_jumps,
                type = 'scatter', mode = 'markers',
                marker = list(color = "#010101", size = 6),
                name = "ECDF points",
                showlegend = F)
    
    # Continuous greenish line
    fig <- fig %>%
      add_trace(x = 1:(max(y)), y = c(rep(1,max(y)-length(cumul)),1-cumul),
                type = 'scatter', mode = 'lines',
                line = list(color = "#BABF33", shape = "linear"),
                name = "Smoothed probability")
    
    # Vertical dashed line at discreteVaR
    fig <- fig %>%
      add_trace(x = c(discreteVaR, discreteVaR),
                y = c(1-cumulDist(discreteVaR), -1),
                type = 'scatter', mode = 'lines',
                line = list(color = "#77C4D5", dash = "dash"),
                name = "SIP using non-smoothed probability"
      )
    
    # Horizontal dashed line
    fig <- fig %>%
      add_trace(x = c(0, max(Stability)), y = 1-c(p, p),
                type = 'scatter', mode = 'lines',
                line = list(color = "#010101", dash = "dash"),
                name = "95th percentile",
                showlegend = F)
    
    # Vertical dashed line at VaR
    fig <- fig %>%
      add_trace(x = c(VaR, VaR), y = c(1-0.05, -1),
                type = 'scatter', mode = 'lines',
                line = list(color = "#FDCE07", dash = "dash"),
                name = "SIP using smoothed probability"
      )
    fig <- layout(
      fig,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Horizon <i>n</i>",
                                standoff = 5),
                   range = c(9, 18),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   nticks=11,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Probability of stable adjustments"),
                   range = rev(1-c(0, 0.1)),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      legend = list(x = 0.02, y = .02,
                    xanchor = "left",
                    yanchor = "bot",
                    font = legendFont,
                    bordercolor = "black", # Set the legend border color
                    borderwidth = 1,
                    bgcolor = "rgba(255, 255, 255, 0.9)"),
      margin = list(t = 30, b=40)
    )%>%config(mathjax = 'cdn')
    fig
  }
}
save_image(fig,paste0(exportPath,"smoothedVSempericalSIP.pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"smoothedVSempericalSIP.pdf"))


####### Fig 2:SIP Homogeneous evolution with nb1 ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height


# adjustable parameters
nb1 <- seq(1,5000)     # only possibility
age1 <- c(60,65,70)   # only 60, 65 and 70 available, can select some are all of them

# not adjustable (for now)
b1=10

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # import  stability surface
  {
    riskStability <- matrix(0, length(age1), length(nb1))
    for (i in seq_along(age1)) {
      name <- paste0("simulatedData/smoothedSIP_age",age1[i],"_b",b1,".rds")
      riskStability[i,] <- readRDS(name)
    }
  }
  
  #plot
  {
    colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
    colors <- rep(colors, length.out = nrow(riskStability))  # ensure enough colors
    
    p <- plot_ly() 
    
    # Add each column of slices as a separate trace
    for (i in 1:nrow(riskStability)) {
      p <- add_trace(
        p,
        x = nb1,
        y = riskStability[i,],
        type = 'scatter',
        mode = 'lines',
        line = list(color = colors[i]),
        name = paste0("Age ", age1[i])
      )
    }
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Stable income period"),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      legend = list(x = 0.98, y = .02,
                    xanchor = "right",
                    yanchor = "bottom",
                    font = legendFont,
                    bordercolor = "black", # Set the legend border color
                    borderwidth = 1,
                    bgcolor = "rgba(255, 255, 255, 0.9)"),
      margin = list(t = 30, b=40)
    )
    p
  }
}
save_image(p,paste0(exportPath,"homoSIP_Smooth.pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"homoSIP_Smooth.pdf"))




####### Fig 3:SIP Heterogeneous Mortality evolution with nb2 ########
#dimensions as percentage of page
w <- .8    #width
h <- .3   #height


# adjustable parameters
nb2 <- seq(1,500)     # only option
age2 <- c(55,60,65,70,75)   # only 55, 60, 65, 70 and 75 available

# not adjustable (for now)
b1 = 10
b2 = 10
nb1 = 100
age1 = 65

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(age2), length(nb2))
    for (i in seq_along(age2)) {
      name <- paste0("simulatedData/smoothedSIP_y1_age",age1, age2[i],"_b",b1,b2,".rds")
      riskStability[i,] <- readRDS(name)
    }
  }
  
  #plot
  {
    colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
    colors <- rep(colors, length.out = nrow(riskStability))  # ensure enough colors
    
    p <- plot_ly() 
    
    # Add each column of slices as a separate trace
    for (i in 1:nrow(riskStability)) {
      p <- add_trace(
        p,
        x = nb2,
        y = riskStability[i,],
        type = 'scatter',
        mode = 'lines',
        line = list(color = colors[i]),
        name = paste0("Age ", age2[i])
      )
    }
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in Group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Stable income period"),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      legend = list(x = 0.02, y = .98,
                    xanchor = "left",
                    yanchor = "top",
                    font = legendFont,
                    bordercolor = "black", # Set the legend border color
                    borderwidth = 1,
                    bgcolor = "rgba(255, 255, 255, 0.9)"),
      margin = list(t = 30, b=40)
    )
    p
  }
}
save_image(p,paste0(exportPath,"SIPMortalityHeteNb.pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPMortalityHeteNb.pdf"))



####### Fig 4:SIP Heterogeneous Wealth evolution with nb2 ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- seq(1,500)     # only option
BMulti <- c(.2,.5,1,2,5) #ratio of benefit2/benefit1

# not adjustable (for now)
b1 = 10
b2 = 10
nb1 = 100
age1 = 65
age2 = 65

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(BMulti), length(nb2))
    for (i in seq_along(BMulti)) {
      name <- paste0("simulatedData/smoothedSIP_y",BMulti[i],"_age",age1, age2,"_b",b1,b2,".rds")
      riskStability[i,] <- readRDS(name)
    }
  }
  
  #plot
  {
    colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
    colors <- rep(colors, length.out = nrow(riskStability))  # ensure enough colors
    
    p <- plot_ly() 
    
    # Add each column of slices as a separate trace
    for (i in 1:nrow(riskStability)) {
      p <- add_trace(
        p,
        x = nb2,
        y = riskStability[i,],
        type = 'scatter',
        mode = 'lines',
        line = list(color = colors[i]),
        name = paste0("<i>y</i> = ", BMulti[i])
      )
    }
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in Group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Stable income period"),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      legend = list(x = 0.02, y = .98,
                    xanchor = "left",
                    yanchor = "top",
                    font = legendFont,
                    bordercolor = "black", # Set the legend border color
                    borderwidth = 1,
                    bgcolor = "rgba(255, 255, 255, 0.9)"),
      margin = list(t = 30, b=40)
    )
    p
  }
}
save_image(p,paste0(exportPath,"SIPWealthHeteNb.pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPWealthHeteNb.pdf"))



####### Fig 5:SIP 2D Plot Wealth Heterogeneity ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- 100      # 50, 100 and 200 available for all age1. 
age2 <- c(60,65,70) # Any amount of value from interval [55,75]

# not adjustable (for now)
age1 <- 65 
nb1 <- 100 
b1=10
b2=10

# Generate plot (put it in full screen before saving for better placement of legend)
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/BaseRisk",
                                   age1,"_",nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  name <- paste("simulatedData/smoothedSIP_nb",nb1,nb2,"_b",b1,b2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # extract age slice
  {
    # get age slices
    age2Vec <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    rowAge <- sapply(age2, function(x)which.min(abs(age2Vec-x)))
    slices <- riskStability[rowAge,]
    
    # get homogeneous value 
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2Vec)
    homoSIP <- riskStability[rowAge,colMultiplier]
  }
  
  #plot
  {
    colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
    colors <- rep(colors, length.out = nrow(slices))  # ensure enough colors
    
    p <- plot_ly() 
    
    # Add each column of slices as a separate trace
    for (i in 1:nrow(slices)) {
      p <- add_trace(
        p,
        x = benefitMultiplier,
        y = slices[i,],
        type = 'scatter',
        mode = 'lines',
        line = list(color = colors[i]),
        name = paste0("Age ", age2[i])
      )
      
      # Add vertical line at Benefit = some value (if desired)
      # Example if you want a vertical reference at max of ages:
      p <- add_trace(
        p,
        x = c(benefitMultiplier[which(max(slices[i,])==slices[i,])],
              benefitMultiplier[which(max(slices[i,])==slices[i,])]),
        y = c(min(slices,riskSmallHomo)-1, max(slices,homoSIP)+1),
        type = 'scatter',
        mode = 'lines',
        line = list(dash = 'dot', color = colors[i]),
        name = paste("Max stability age", age2[i]),
        showlegend=F
      )
    }
    
    # Add horizontal lines
    p <- add_trace(
      p,
      x = c(min(benefitMultiplier), max(benefitMultiplier)),
      y = c(homoSIP, homoSIP),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'gray'),
      name = paste0(nb1+nb2, " members from Group 1"),
      showlegend=T
    )
    
    p <- add_trace(
      p,
      x = c(min(benefitMultiplier), max(benefitMultiplier)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste0(nb1, " members from Group 1"),
      showlegend=T
    )
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Initial benefit of in Group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),showgrid=FALSE,
                   ticks    = "outside",
                   type = "log",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Stable income period"),
                   titlefont = axisFont,
                   range = range(slices,riskSmallHomo)*c(.95,1.05),
                   tickfont = list(size = 12),showgrid=FALSE,
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      legend = list(x = 0.02, y = .98,
                    xanchor = "left",
                    yanchor = "top",
                    font = legendFont,
                    bordercolor = "black", # Set the legend border color
                    borderwidth = 1,
                    bgcolor = "rgba(255, 255, 255, 0.9)"),
      margin = list(t = 30, b=40)
    )
    p
  }
}
save_image(p,paste0(exportPath,"SIPWealthHeteMortality.pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPWealthHeteMortality.pdf"))




####### Fig 6:SIP 2D Plot Age Heterogeneity ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- 100       # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
diffWealth <- c(.2, .5, 1, 2, 5) # Any amount of value from interval [.1,10], Benefit of Group 2 compared to group 1

# not adjustable (for now)
age1 <- 65 
nb1 <- 100 
b1=10
b2=10

# Generate plot (put it in full screen before saving for better placement of legend)
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/BaseRisk",
                                 age1,"_",nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  name <- paste("simulatedData/smoothedSIP_nb",nb1,nb2,"_b",b1,b2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # extract wealth slice
  {
    # get wealth slices
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    colMultiplier <- sapply(diffWealth, function(x)which.min(abs(benefitMultiplier-x)))
    slices <- riskStability[,colMultiplier]
    
    # get homogeneous value 
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoSIP <- riskStability[rowAge,colMultiplier]
  }
  
  #plot
  {
    colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
    colors <- rep(colors, length.out = ncol(slices))  # ensure enough colors
    p <- plot_ly() 
    
    # Add each column of slices as a separate trace
    for (i in 1:ncol(slices)) {
      p <- add_trace(
        p,
        x = age2,
        y = slices[, i],
        type = 'scatter',
        mode = 'lines',
        line = list(color = colors[i]),
        name = paste0("<i>y</i> = ", diffWealth[i])
      )
    }
    
    # Add horizontal lines
    p <- add_trace(
      p,
      x = c(min(age2), max(age2)),
      y = c(homoSIP, homoSIP),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'gray'),
      name = paste0(nb1+nb2, " members from Group 1")
    )
    
    p <- add_trace(
      p,
      x = c(min(age2), max(age2)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste0(nb1, " members from Group 1")
    )
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Age of members in Group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Stable income period"),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      legend = list(x = 0.98, y = .98,
                    xanchor = "right",
                    yanchor = "top",
                    font = legendFont,
                    bordercolor = "black", # Set the legend border color
                    borderwidth = 1,
                    bgcolor = "rgba(255, 255, 255, 0.9)"),
      margin = list(t = 30, b=40)
    )
    p
  }
}
save_image(p,paste0(exportPath,"SIPMortalityHeteWealth.pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPMortalityHeteWealth.pdf"))



####### Fig 7:SIP Contour Plot Group1's Perspective ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- 100     # 50, 100 and 200 available for all age1. 

# not adjustable (for now)
age1 <- 65 
nb1 <- 100 
b1=10
b2=10

# Generate plot
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/BaseRisk",
                                 age1,"_",nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  name <- paste("simulatedData/smoothedSIP_nb",nb1,nb2,"_b",b1,b2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # compute better and worse areas
  {
    ## get stability when group 1 and 2 are homogeneous from imported surface
    benMultiToExtract <- 1
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoSIP <- riskStability[rowAge,colMultiplier]
    
    ## extract better (green) area
    dfbetterSIP <- list(SIP = riskStability[riskStability>=homoSIP])
    dfbetterSIP$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  >=homoSIP]
    dfbetterSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability>=homoSIP]
    
    ## extract worst (red) area
    dfworsteSIP <- list(SIP = riskStability[riskStability<=riskSmallHomo])
    dfworsteSIP$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  <=riskSmallHomo]
    dfworsteSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SIP,
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
      add_trace(
        data = dfbetterSIP,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
        name = "Preferred Region",
        inherit = FALSE
      )%>%
      add_trace(
        data = dfworsteSIP,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      ) %>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",   # uniform background color
        paper_bgcolor = "white",  # outside background
        xaxis = list(title = list(text = "Age of members in Group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of Group 2",
                                  standoff = 5),
                     type = "log",
                     showgrid = FALSE,
                     range= log(c(min(df$benefitMultiplier),
                                  max(df$benefitMultiplier)),10),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b=40),
        showlegend = F
      )
    p
  }
}
save_image(p,paste0(exportPath,"SIPContour1Perspective",nb1,nb2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPContour1Perspective",nb1,nb2,".pdf"))



####### Fig 8-9-10:SIP Contour Plot Both Groups' Perspective ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

for(nb2 in c(100,50,200)){
  # Generate plot
  {
    # not adjustable
    age1 <- 65 
    nb1 <- 100 
    b1=10
    b2=10
    
    # import base stability when group 1 is on its own
    riskSmallHomo <- readRDS(paste0("simulatedData/BaseRisk", 
                                    age1,"_",nb1,".rds"))
    
    # import smoothed stability surface
    name <- paste0("simulatedData/smoothedSIP_nb",nb1,nb2,"_b",b1,b2,".rds")
    riskStability <- readRDS(name)
    
    # import smoothed stability when group 2 is on its own
    name <- paste0("simulatedData/smoothedSIP_nb",0,nb2,"_b",b1,b2,".rds")
    riskStabilitySmallHomo <- readRDS(name)
    
    # get worse areas (no better area possible for both at the same time)
    {
      age2 <- seq(55, 75, by = .1)
      benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
      
      
      dfworsteSIP <- list(SIP = riskStability[riskStability<=riskStabilitySmallHomo
                                              |riskStability<=riskSmallHomo])
      dfworsteSIP$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                  ncol = length(benefitMultiplier))[riskStability<=riskStabilitySmallHomo
                                                                    |riskStability<=riskSmallHomo]
      dfworsteSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                              nrow = length(age2),
                                              byrow = T)[riskStability<=riskStabilitySmallHomo
                                                         |riskStability<=riskSmallHomo]
      
    }
    
    # contour Plot
    {
      # Prepare data in long format
      df <- melt(riskStability)
      colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
      df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
      df$age2 <- age2[df$ageIndex]
      
      # Create contour plot
      p <- plot_ly(
        data = df,
        x = ~age2,
        y = ~benefitMultiplier,
        z = ~SIP,
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
        add_trace(
          data = dfworsteSIP,
          x = ~age2,
          y = ~benefitMultiplier,
          type = "scatter",
          mode = "markers",
          marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
          name = "No No Region",
          inherit = FALSE
        )%>%
        layout(
          font = list(family = fontType),
          plot_bgcolor = "lightgrey",   # uniform background color
          paper_bgcolor = "white",  # outside background
          xaxis = list(title = list(text = "Age of members in Group 2",
                                    standoff = 5),
                       showgrid = FALSE,
                       range = c(min(df$age2), max(df$age2)),
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE
          ),
          yaxis = list(title = list(text = "Initial benefit of Group 2",
                                    standoff = 5),
                       type = "log",
                       showgrid = FALSE,
                       range= log(c(min(df$benefitMultiplier),
                                    max(df$benefitMultiplier)),10),
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE
          ),
          margin = list(t = 50, b=40),
          showlegend = F
        )
      p
    }
  }
  save_image(p,paste0(exportPath,"SIPContour2Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPContour2Perspective",nb1,nb2,".pdf"))
  
}

####### Fig 11:SIP Contour Plot Group1's Perspective Risky Asset ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- 100     

# not adjustable (for now)
age1 <- 65 
nb1 <- 100 
b1=10
b2=10

# Generate plot
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste0("simulatedData/BaseRisk_risky_"
                                  ,age1,"_",nb1,".rds"))
  
  # import smoothed stability surface
  name <- paste0("simulatedData/smoothedSIP_risky_nb",nb1,nb2,"_b",b1,b2,".rds")
  riskStability <- readRDS(name)
  
  # compute better and worse areas
  {
    ## get stability when group 1 and 2 are homogeneous from imported surface
    benMultiToExtract <- 1
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoSIP <- riskStability[rowAge,colMultiplier]
    
    ## extract better (green) area
    dfbetterSIP <- list(SIP = riskStability[riskStability>=homoSIP])
    dfbetterSIP$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  >=homoSIP]
    dfbetterSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability>=homoSIP]
    
    ## extract worst (red) area
    dfworsteSIP <- list(SIP = riskStability[riskStability<=riskSmallHomo])
    dfworsteSIP$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  <=riskSmallHomo]
    dfworsteSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SIP,
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
      add_trace(
        data = dfbetterSIP,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
        name = "Preferred Region",
        inherit = FALSE
      )%>%
      add_trace(
        data = dfworsteSIP,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      )%>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",   # uniform background color
        paper_bgcolor = "white",  # outside background
        xaxis = list(title = list(text = "Age of members in Group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of Group 2",
                                  standoff = 5),
                     type = "log",
                     showgrid = FALSE,
                     range= log(c(min(df$benefitMultiplier),
                                  max(df$benefitMultiplier)),10),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b=40),
        showlegend = F
      )
    p
  }
}
save_image(p,paste0(exportPath,"SIPContourRiskyAsset",nb1,nb2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPContourRiskyAsset",nb1,nb2,".pdf"))

####### Fig 12:SIP Contour Plot Both Groups' Perspective Risky Asset########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- 100     # only 100

# not adjustable (for now)
age1 <- 65 
nb1 <- 100 
b1=10
b2=10

# Generate plot
{
  
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste0("simulatedData/BaseRisk_risky_", 
                                  age1,"_",nb1,".rds"))
  
  # import smoothed stability surface
  name <- paste0("simulatedData/smoothedSIP_risky_nb",nb1,nb2,"_b",b1,b2,".rds")
  riskStability <- readRDS(name)
  
  # import smoothed stability when group 2 is on its own
  name <- paste0("simulatedData/smoothedSIP_risky_nb",0,nb2,"_b",b1,b2,".rds")
  riskStabilitySmallHomo <- readRDS(name)
  
  # get worse areas (no better area possible for both at the same time)
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    
    
    dfworsteSIP <- list(SIP = riskStability[riskStability<=riskStabilitySmallHomo
                                            |riskStability<=riskSmallHomo])
    dfworsteSIP$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability<=riskStabilitySmallHomo
                                                                  |riskStability<=riskSmallHomo]
    dfworsteSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=riskStabilitySmallHomo
                                                       |riskStability<=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SIP,
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
      add_trace(
        data = dfworsteSIP,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      )%>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",   # uniform background color
        paper_bgcolor = "white",  # outside background
        xaxis = list(title = list(text = "Age of members in Group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of Group 2",
                                  standoff = 5),
                     type = "log",
                     showgrid = FALSE,
                     range= log(c(min(df$benefitMultiplier),
                                  max(df$benefitMultiplier)),10),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b=40),
        showlegend = F
      )
    p
  }
}
save_image(p,paste0(exportPath,"SIPContourRiskyAsset2Perspective",nb1,nb2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPContourRiskyAsset2Perspective",nb1,nb2,".pdf"))



####### Fig 13:SIP Different m Parameter Mortality ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

nb2 <- 100 #no other option
nb1 <- 100 #no other option

{
  #axis
  mVec <- seq(75,95, by = .1)
  b2 <- seq(6,14, by=.1)
  
  name <- paste0("simulatedData/smoothedSIP_ParrallelComputingMortalityParam",nb1,".rds")
  SIP <- readRDS(name)
  
  #slice plot
  {
    # adjustable parameters
    nb2 <- 100       # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
    age1 <- 65      # only 60, 65 and 70 available
    diff_b2 <- c(6, 8, 10, 12, 14) # Any amount of value from interval [.1,10], Benefit of Group 2 compared to group 1
    
    nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)
    
    # Generate plot (put it in full screen before saving for better placement of legend)
    {
      # extract wealth slice
      {
        # get wealth slices
        age2 <- seq(55, 75, by = .1)
        colMultiplier <- sapply(diff_b2, function(x)which.min(abs(b2-x)))
        slices <- SIP[,colMultiplier]
      }
      
      #plot w.r.t m
      {
      #   colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
      #   colors <- rep(colors, length.out = ncol(slices))  # ensure enough colors
      #   p <- plot_ly() 
      #   
      #   # Add each column of slices as a separate trace
      #   for (i in 1:ncol(slices)) {
      #     p <- add_trace(
      #       p,
      #       x = mVec,
      #       y = slices[, i],
      #       type = 'scatter',
      #       mode = 'lines',
      #       line = list(color = colors[i]),
      #       name = paste0("<i>b<sub>2</sub></i> = ", diff_b2[i])
      #     )
      #   }
      #   # Final layout
      #   p <- layout(
      #     p,
      #     font = list(family = fontType),
      #     xaxis = list(title = list(text = "Parameter m for members in Group 2",
      #                               standoff = 5),
      #                  titlefont = axisFont,
      #                  tickfont = list(size = 12),
      #                  ticks    = "outside",
      #                  ticklen  = 8,
      #                  showline = TRUE, mirror = TRUE, zeroline = FALSE),
      #     yaxis = list(title = list(text = "Stable income period"),
      #                  titlefont = axisFont,
      #                  tickfont = list(size = 12),
      #                  ticks    = "outside",
      #                  ticklen  = 8,
      #                  showline = TRUE, mirror = TRUE, zeroline = FALSE),
      #     legend = list(x = 0.02, y = .98,
      #                   xanchor = "left",
      #                   yanchor = "top",
      #                   font = legendFont,
      #                   bordercolor = "black", # Set the legend border color
      #                   borderwidth = 1,
      #                   bgcolor = "rgba(255, 255, 255, 0.9)"),
      #     margin = list(t = 30, b=40)
      #   )
      #   p
      }
      # save_image(p,paste0(exportPath,"SIPHeteMortDist_m.pdf"),
      #            width = w/2*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
      # browseURL(paste0(exportPath,"SIPHeteMortDist_m.pdf"))
      
      #plot w.r.t age
      {
        colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
        colors <- rep(colors, length.out = ncol(slices))  # ensure enough colors
        p <- plot_ly() 
        
        # Add each column of slices as a separate trace
        for (i in 1:ncol(slices)) {
          p <- add_trace(
            p,
            x = age2,
            y = rev(slices[, i]),
            type = 'scatter',
            mode = 'lines',
            line = list(color = colors[i]),
            name = paste0("<i>b<sub>2</sub></i> = ", diff_b2[i])
          )
        }
        # Final layout
        p <- layout(
          p,
          font = list(family = fontType),
          xaxis = list(title = list(text = "Age of members in Group 2",
                                    standoff = 5),
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          yaxis = list(title = list(text = "Stable income period"),
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          legend = list(x = 0.98, y = .98,
                        xanchor = "right",
                        yanchor = "top",
                        font = legendFont,
                        bordercolor = "black", # Set the legend border color
                        borderwidth = 1,
                        bgcolor = "rgba(255, 255, 255, 0.9)"),
          margin = list(t = 30, b=40)
        )
        p
      }
    }
  }
}
save_image(p,paste0(exportPath,"SIPHeteMortDist_age.pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPHeteMortDist_age.pdf"))
###############################################################################
# Section 4 Approx SIP 
###############################################################################
####### Fig 14:Approx SIP Contour Plot Group1's Perspective ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- 100     # 50, 100 and 200 to match section 3.

age1 <- 65 
nb1 <- 100 
b1=10
b2=10

{
  ##Compute riskStabilityApprox
  {
    
    benefit1 <- 1000
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    benefit2 <- benefit1*benefitMultiplier
    SIP_Approx <- matrix(0,length(age2),length(benefitMultiplier))
    for (i in seq_along(age2)) {
      for (j in seq_along(benefit2)) {
        SIP_Approx[i,j] <- approxSIP2Pop(nb1, nb2, age1, age2[i], 
                                                  benefit1,benefit2[j],
                                                  b1 = b1, b2=b2)
      }
    }
    
    
    # SD surface for smoothing purposes
    {
      sdProxy <- matrix(0, length(age2), length(benefitMultiplier))
      for (i in seq_along(age2)) {
        sdProxy[i,] <- sapply(benefitMultiplier,
                              function(x)SD1Periode(age1,1000,nb1,
                                                    age2[i],1000*x,nb2,
                                                    b1=b1, b2=b2))
      }
    }
    
    grid <- expand.grid(age = age2, benMulti = benefitMultiplier)
    
    ### Flatten matrices column-wise (assuming they match grid layout)
    grid$sdProxy <- as.vector(sdProxy)
    grid$SIP <- as.vector(SIP_Approx)
    
    ### fit model
    
    gam_fit <- gam(SIP ~ te(age, benMulti, sdProxy, k=5), data = grid)
    summary(gam_fit)
    
    grid$SIP_smooth <- predict(gam_fit, newdata = grid)
    
    SIP_ApproxSmooth <- matrix(grid$SIP_smooth, nrow = length(age2),
                                        ncol = length(benefitMultiplier))
    
    
  }
  #better and worse area
  {
    riskSmallHomoApprox <- approxSIP2Pop(nb1, 0, age1, age1, benefit1,
                                         benefit1)
    
    ## get stability when group 1 and 2 are homogeneous from imported surface
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoSIPApprox <- SIP_ApproxSmooth[rowAge,colMultiplier]
    
    ## extract better (green) area
    dfbetterSIPApprox <- list(SIP = SIP_ApproxSmooth[SIP_ApproxSmooth>=homoSIPApprox])
    dfbetterSIPApprox$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                      ncol = length(benefitMultiplier))[SIP_ApproxSmooth
                                                                        >=homoSIPApprox]
    dfbetterSIPApprox$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                                  nrow = length(age2),
                                                  byrow = T)[SIP_ApproxSmooth>=homoSIPApprox]
    
    ## extract worst (red) area
    dfworsteSIPApprox <- list(SIP = SIP_ApproxSmooth[SIP_ApproxSmooth<=riskSmallHomoApprox])
    dfworsteSIPApprox$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                      ncol = length(benefitMultiplier))[SIP_ApproxSmooth
                                                                        <=riskSmallHomoApprox]
    dfworsteSIPApprox$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                                  nrow = length(age2),
                                                  byrow = T)[SIP_ApproxSmooth<=riskSmallHomoApprox]
    
  }
  ## approx surface contour
  {
    # Prepare data in long format
    df <- melt(SIP_ApproxSmooth)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SIP,
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
      add_trace(
        data = dfbetterSIPApprox,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
        name = "Preferred Region",
        inherit = FALSE
      )%>%
      add_trace(
        data = dfworsteSIPApprox,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      ) %>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",   # uniform background color
        paper_bgcolor = "white",  # outside background
        xaxis = list(title = list(text = "Age of members in Group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of Group 2",
                                  standoff = 5),
                     type = "log",
                     showgrid = FALSE,
                     range= log(c(min(df$benefitMultiplier),
                                  max(df$benefitMultiplier)),10),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b=40),
        showlegend = F
      )
    p
  }
}
save_image(p,paste0(exportPath,"ApproxSIPContour1Perspective",nb1,nb2,"_b",b1,b2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"ApproxSIPContour1Perspective",nb1,nb2,"_b",b1,b2,".pdf"))

####### Fig 15:Approx SIP Contour Plot Both Groups' Perspective ########
#dimensions as percentage of page
w <- .8    #width
h <- .3  #height

# adjustable parameters
nb2 <- 100     # 50, 100 and 200 to match section 3.

age1 <- 65 
nb1 <- 100 
b1=10
b2=10


{
  ####Compute riskStabilityApprox
  {
    
    benefit1 <- 1000
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    benefit2 <- benefit1*benefitMultiplier
    SIP_Approx <- matrix(0,length(age2),length(benefitMultiplier))
    for (i in seq_along(age2)) {
      for (j in seq_along(benefit2)) {
        SIP_Approx[i,j] <- approxSIP2Pop(nb1, nb2, age1, age2[i],
                                         benefit1,benefit2[j], 
                                         b1=b1, b2=b2)
      }
    }
    
    
    
    SIP_SmallHomoApprox <- matrix(0,length(age2),length(benefitMultiplier))
    for (i in seq_along(age2)) {
        SIP_SmallHomoApprox[i,] <- approxSIP2Pop(0, nb2, age1, age2[i], 
                                                  benefit1,benefit1, 
                                                  b1=b1, b2=b2)
    }
    
    
    # SD surface for smoothing purposes
    {
      
      sdProxy <- matrix(0, length(age2), length(benefitMultiplier))
      for (i in seq_along(age2)) {
        sdProxy[i,] <- sapply(benefitMultiplier,
                              function(x)SD1Periode(age1,1000,nb1,
                                                    age2[i],1000*x,nb2,
                                                    b1=b1, b2=b2))
      }
      sdProxySmallHomo <- matrix(0,length(age2),length(benefitMultiplier))
      for (i in seq_along(age2)) {
        sdProxySmallHomo[i,] <- SD1Periode(age1, 1000, 0,
                                           age2[i],1000,nb2,
                                           b1=b1, b2=b2)
      }
      
    }
    
    grid <- expand.grid(age = age2, benMulti = benefitMultiplier)
    
    ### Flatten matrices column-wise (assuming they match grid layout)
    grid$sdProxy <- as.vector(sdProxy)
    grid$SIP <- as.vector(SIP_Approx)
    
    ### fit model
    
    gam_fit <- gam(SIP ~ te(age, benMulti, sdProxy, k=5), data = grid)
    summary(gam_fit)
    
    grid$SIP_smooth <- predict(gam_fit, newdata = grid)
    
    SIP_ApproxSmooth <- matrix(grid$SIP_smooth, nrow = length(age2),
                                        ncol = length(benefitMultiplier))
    
    
    ### add SmallHomo (assuming they match grid layout)
    grid$sdProxySmallHomo <- as.vector(sdProxySmallHomo)
    grid$SIPsmallHomo <- as.vector(SIP_SmallHomoApprox)
    
    ### fit model
    
    gam_fit <- gam(SIPsmallHomo ~ te(age, k=3), data = grid)
    summary(gam_fit)
    
    grid$SIP_smooth <- predict(gam_fit, newdata = grid)
    
    SIP_SmallHomoApproxSmooth <- matrix(grid$SIP_smooth, nrow = length(age2),
                                                 ncol = length(benefitMultiplier))
    
  }
  
  
  #worse area
  {
    riskSmallHomoApprox <- approxSIP2Pop(nb1, 0, age1, age1, benefit1,
                                         benefit1)
    
    dfworsteSIPApprox <- list(SIP = SIP_ApproxSmooth[SIP_ApproxSmooth<=SIP_SmallHomoApproxSmooth
                                                     |SIP_ApproxSmooth<=riskSmallHomoApprox])
    dfworsteSIPApprox$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                      ncol = length(benefitMultiplier))[SIP_ApproxSmooth<=SIP_SmallHomoApproxSmooth
                                                                        |SIP_ApproxSmooth<=riskSmallHomoApprox]
    dfworsteSIPApprox$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                                  nrow = length(age2),
                                                  byrow = T)[SIP_ApproxSmooth<=SIP_SmallHomoApproxSmooth
                                                             |SIP_ApproxSmooth<=riskSmallHomoApprox]
    
  }
  ## approx surface contour
  {
    # Prepare data in long format
    df <- melt(SIP_ApproxSmooth)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SIP,
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
      add_trace(
        data = dfworsteSIPApprox,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      )%>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",   # uniform background color
        paper_bgcolor = "white",  # outside background
        xaxis = list(title = list(text = "Age of members in Group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of Group 2",
                                  standoff = 5),
                     type = "log",
                     showgrid = FALSE,
                     range= log(c(min(df$benefitMultiplier),
                                  max(df$benefitMultiplier)),10),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b=40),
        showlegend = F
      )
    p
  }
}
save_image(p,paste0(exportPath,"ApproxSIPContour2Perspective",nb1,nb2,"_b",b1,b2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"ApproxSIPContour2Perspective",nb1,nb2,"_b",b1,b2,".pdf"))



###############################################################################
# Section 5 Mitigation
###############################################################################
####### Fig 16:Example pool post mitigation SIP only #######
#dimensions as percentage of page
w <- .7    #width
h <- .25  #height

# adjustable parameters
nb2 <- 100        # 50, 100 and 200 available for all

age1 <- 65        #only 60, 65 and 70 available
benefit1 <- 100
age2 <- 70
benefit2 <- 500
b1 = 10
b2 = 10


{
  asset1 <- as.vector(benefit1*annuity(age1,.02))
  asset2 <- as.vector(benefit2*annuity(age2,.02))
  nb1 <- 100      # Don't change
  
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste0("simulatedData/BaseRisk", 
                                  age1,"_",nb1,".rds"))
  
  # import smoothed stability surface
  name <- paste0("simulatedData/smoothedSIP_nb",nb1,nb2,"_b",b1,b2,".rds")
  riskStability <- readRDS(name)
  
  # import smoothed stability when group 2 is on its own
  name <- paste0("simulatedData/smoothedSIP_nb",0,nb2,"_b",b1,b2,".rds")
  riskStabilitySmallHomo <- readRDS(name)
  
  # get worse areas (no better area possible for both at the same time)
  {
    age2Vec <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    
    
    dfworsteSIP <- list(SIP = riskStability[riskStability<=riskStabilitySmallHomo
                                            |riskStability<=riskSmallHomo])
    dfworsteSIP$age2Vec <-  matrix(rep(age2Vec, length(benefitMultiplier)),
                                   ncol = length(benefitMultiplier))[riskStability<=riskStabilitySmallHomo
                                                                     |riskStability<=riskSmallHomo]
    dfworsteSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2Vec)),
                                            nrow = length(age2Vec),
                                            byrow = T)[riskStability<=riskStabilitySmallHomo
                                                       |riskStability<=riskSmallHomo]
    
  }
  
  {
    totAsset <- (nb1*asset1+nb2*asset2)
    currentSplit <- (nb1*asset1)/totAsset
    currentY <- benefit2/benefit1
    
    newY <- benefitMultiplier[which.max(riskStability[age2Vec ==age2,])]
  }
  
  {#table info
    pps1natural <- asset1/benefit1
    pps2natural <- asset2/benefit2
    naturalSIP <- riskStability[age2Vec ==age2,
                                which.min(abs(benefitMultiplier-currentY))]
    
    optimizedSIP <- riskStability[age2Vec == age2,
                                  which.min(abs(benefitMultiplier-newY))]
    
    SIP1benchmark <- riskSmallHomo
    SIP2benchmark <- riskStabilitySmallHomo[age2Vec ==age2,
                                            which.min(abs(benefitMultiplier-currentY))]
    
    newAssetRatio <- as.vector(newY*annuity(age2,.02)/annuity(age1,.02))
    newAsset1 <- totAsset/(nb1+nb2*newAssetRatio)
    newAsset2 <- (totAsset-nb1*newAsset1)/nb2
    newBen1 <- as.vector(newAsset1/annuity(age1, .02))
    newBen2 <- as.vector(newAsset2/annuity(age2, .02))
    pps1optim <- asset1/newBen1
    pps2optim <- asset2/newBen2
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2Vec <- age2Vec[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2Vec,
      y = ~benefitMultiplier,
      z = ~SIP,
      type = "contour",
      showscale = FALSE,
      contours = list(
        coloring = "lines",
        showlabels = TRUE
      ),
      line = list(smoothing = 0),
      colorscale = list(c(0, "black"), c(1, "black")),
      reversescale = FALSE
    ) %>%
      add_trace(
        data = dfworsteSIP,
        x = ~age2Vec,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      ) %>%
      add_trace(
        x = c(age2),
        y = c(currentY),
        type = "scatter",
        mode = "markers",
        marker = list(color = "blue", size = 7, symbol = "circle"),
        name = "Current point",
        inherit = FALSE
      )%>%
      add_trace(
        x = c(age2),
        y = c(newY),
        type = "scatter",
        mode = "markers",
        marker = list(color = "green", size = 7, symbol = "circle"),
        name = "Current point",
        inherit = FALSE
      ) %>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",
        paper_bgcolor = "white",
        xaxis = list(
          title = list(text = "Age of members in Group 2", standoff = 5),
          showgrid = FALSE,
          range = c(min(df$age2Vec), max(df$age2Vec)),
          titlefont = axisFont,
          tickfont = list(size = 12),
          ticks = "outside",
          ticklen = 8,
          showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(
          title = list(text = "Initial benefit of Group 2", standoff = 5),
          type = "log",
          showgrid = FALSE,
          range = log10(c(min(df$benefitMultiplier), max(df$benefitMultiplier))),
          titlefont = axisFont,
          tickfont = list(size = 12),
          ticks = "outside",
          ticklen = 8,
          showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b = 40),
        showlegend = FALSE,
        annotations = list(
          list(
            # ⚠️ Use log10(y) for both head and tail because y-axis is logarithmic
            x = age2, y = log10(newY),
            ax = age2, ay = log10(currentY*1.05),
            xref = "x", yref = "y",
            axref = "x", ayref = "y",
            showarrow = TRUE,
            arrowhead = 3,
            arrowsize = 1.5,
            arrowwidth = 2,
            arrowcolor = "black",
            standoff = 2,
            layer = "above"   # ensures it’s drawn over contour lines
          )
        )
      )
    p
  }
}
save_image(p,paste0(exportPath,"SIPContourPostMitigate",nb1,nb2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPContourPostMitigate",nb1,nb2,".pdf"))


####### Fig 17:Initial Benefit and SIP Preferences #######
# dimensions as percentage of page
w <- .7    #width
h <- .25  #height

# not adjustable parameters
nb1 <- 100
age1 <- 65
benefit1 <- 100
nb2 <- 100
age2 <- 70
benefit2 <- 500

alpha <- 2 #risk aversion level

b1=10 #mortality param
b2=10

# contour plot Fig 17
{
  
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste0("simulatedData/BaseRisk", 
                                  age1,"_",nb1,".rds"))
  
  # import smoothed stability surface
  name <- paste0("simulatedData/smoothedSIP_nb",nb1,nb2,"_b",b1,b2,".rds")
  riskStability <- readRDS(name)
  
  # import smoothed stability when group 2 is on its own
  name <- paste0("simulatedData/smoothedSIP_nb",0,nb2,"_b",b1,b2,".rds")
  riskStabilitySmallHomo <- readRDS(name)[,1]
  
  # import smoothed stability when group 2 is on its own
  name <- paste0("simulatedData/smoothedSIP_nb",0,nb2,"_b",b1,b2,".rds")
  riskStabilitySmallHomo <- readRDS(name)
  
  asset1 <- as.vector(benefit1*annuity(age1,.02))
  asset2 <- as.vector(benefit2*annuity(age2,.02))
  percBenefit1 <- benefit1/asset1
  percBenefit2 <- benefit2/asset2
  
  
  age2Vec <- seq(55, 75, by = .1)
  benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
  
  SIP1 <- riskSmallHomo
  SIP2 <- riskStabilitySmallHomo[age2==age2Vec,1]
  
  
  utile1 <- utilityFn(SIP1,percBenefit1, alpha)
  SIPIndif1 <- seq(1,30,.5)
  indifferenceCurve1 <- utilityCurve(utile1,SIPIndif1, alpha)
  
  utile2 <- utilityFn(SIP2,percBenefit2, alpha)
  SIPIndif2 <- seq(1,30,.5)
  indifferenceCurve2 <- utilityCurve(utile2,SIPIndif2, alpha)
  
  baseSplit <- (asset1*nb1)/(asset1*nb1+asset2*nb2)
  baseCase <- meanVarMatrix(nb1 = nb1, age1 = age1,
                            asset1 = asset1,
                            nb2 = nb2, age2 = age2,
                            asset2 = asset2,
                            matSIP = riskStability,
                            assetSplit = baseSplit)
  
  splits <- seq(0.01,.99,0.001)
  traceVec <- sapply(splits, function(a)  meanVarMatrix(nb1 = nb1, age1 = age1,
                                                        asset1 = asset1,
                                                        nb2 = nb2, age2 = age2,
                                                        asset2 = asset2,
                                                        matSIP = riskStability,
                                                        assetSplit = a))
  # Create grid
  grid1 <- list(x = traceVec[1,], y = traceVec[2,], splits = splits)
  #fit model
  gam_fit1 <- gam(x ~ te(y, k=20), data = grid1)
  summary(gam_fit1)
  grid1$smoothyX <- predict(gam_fit1, newdata = grid1)
  # Create grid
  grid2 <- list(x = traceVec[3,], y = traceVec[4,], splits = splits)
  #fit model
  gam_fit2 <- gam(x ~ te(y, k=20), data = grid2)
  summary(gam_fit2)
  grid2$smoothyX <- predict(gam_fit2, newdata = grid2)
  
  #get worst points
  {
    #extend indifference curve (check for all benefit level of SIP curve)
    indifferenceCurve1Extended <- utilityCurve(utile1,grid1$smoothyX,alpha)
    indifferenceCurve2Extended <- utilityCurve(utile2,grid2$smoothyX,alpha)
    
    #compare for all benefit level
    diff1 <- grid1$y - indifferenceCurve1Extended
    diff2 <- grid2$y - indifferenceCurve2Extended
    
    noGoodIndex1 <- which(diff1<=0)
    noGoodIndex2 <- which(diff2<=0)
    
  }
  
  goodRebalanceId <- seq_along(splits)[!seq_along(splits)%in%c(noGoodIndex1,
                                                               noGoodIndex2)]
  minY <- (grid2$y[min(goodRebalanceId)]*asset2)/(grid1$y[min(goodRebalanceId)]*asset1)
  maxY <- (grid2$y[max(goodRebalanceId)]*asset2)/(grid1$y[max(goodRebalanceId)]*asset1)
  
  currentY <- benefit2/benefit1
  
  # get worse areas (no better area possible for both at the same time)
  {
    
    dfworsteSIP <- list(SIP = riskStability[riskStability<=riskStabilitySmallHomo
                                            |riskStability<=riskSmallHomo])
    dfworsteSIP$age2Vec <-  matrix(rep(age2Vec, length(benefitMultiplier)),
                                   ncol = length(benefitMultiplier))[riskStability<=riskStabilitySmallHomo
                                                                     |riskStability<=riskSmallHomo]
    dfworsteSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2Vec)),
                                            nrow = length(age2Vec),
                                            byrow = T)[riskStability<=riskStabilitySmallHomo
                                                       |riskStability<=riskSmallHomo]
    
  }
  
  #table info
  {
    totAsset <- (nb1*asset1+nb2*asset2)
    
    
    benchmarkObjFn1 <- utile1
    benchmarkObjFn2 <- utile2
    SIP1benchmark <- riskSmallHomo
    SIP2benchmark <- riskStabilitySmallHomo[age2Vec ==age2,
                                            which.min(abs(benefitMultiplier-currentY))]
    
    pps1natural <- asset1/benefit1
    pps2natural <- asset2/benefit2
    naturalSIP <- riskStability[age2Vec ==age2,
                                which.min(abs(benefitMultiplier-currentY))]
    naturalObjFn1 <- utilityFn(naturalSIP,1/pps1natural)
    naturalObjFn2 <- utilityFn(naturalSIP,1/pps2natural)
    
    
    lowLimSIP <- riskStability[age2Vec == age2,
                               which.min(abs(benefitMultiplier-minY))]
    lowLimAssetRatio <- as.vector(minY*annuity(age2,.02)/annuity(age1,.02))
    lowLimAsset1 <- totAsset/(nb1+nb2*lowLimAssetRatio)
    lowLimAsset2 <- (totAsset-nb1*lowLimAsset1)/nb2
    lowLimBen1 <- as.vector(lowLimAsset1/annuity(age1, .02))
    lowLimBen2 <- as.vector(lowLimAsset2/annuity(age2, .02))
    pps1lowLim <- asset1/lowLimBen1
    pps2lowLim <- asset2/lowLimBen2
    lowLimObjFn1 <- utilityFn(lowLimSIP,1/pps1lowLim)
    lowLimObjFn2 <- utilityFn(lowLimSIP,1/pps2lowLim)
    
    upLimSIP <- riskStability[age2Vec == age2,
                              which.min(abs(benefitMultiplier-maxY))]
    upLimAssetRatio <- as.vector(maxY*annuity(age2,.02)/annuity(age1,.02))
    upLimAsset1 <- totAsset/(nb1+nb2*upLimAssetRatio)
    upLimAsset2 <- (totAsset-nb1*upLimAsset1)/nb2
    upLimBen1 <- as.vector(upLimAsset1/annuity(age1, .02))
    upLimBen2 <- as.vector(upLimAsset2/annuity(age2, .02))
    pps1upLim <- asset1/upLimBen1
    pps2upLim <- asset2/upLimBen2
    upLimObjFn1 <- utilityFn(upLimSIP,1/pps1upLim)
    upLimObjFn2 <- utilityFn(upLimSIP,1/pps2upLim)
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2Vec <- age2Vec[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2Vec,
      y = ~benefitMultiplier,
      z = ~SIP,
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
      add_trace(
        data = dfworsteSIP,
        x = ~age2Vec,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      ) %>%
      add_trace(
        x = c(age2),
        y = c(currentY),
        type = "scatter",
        mode = "markers",
        marker = list(color = "blue", size = 5, symbol = "circle"),
        name = "Current point",
        inherit = FALSE
      ) %>%
      add_trace(
        x = c(age2,age2),
        y = c(minY,maxY),
        type = "scatter",
        mode = "lines",
        line = list(color = "green", width = 5),
        name = "good reallocation",
        inherit = FALSE
      )%>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",   # uniform background color
        paper_bgcolor = "white",  # outside background
        xaxis = list(title = list(text = "Age of members in Group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2Vec), max(df$age2Vec)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of Group 2",
                                  standoff = 5),
                     type = "log",
                     showgrid = FALSE,
                     range= log(c(min(df$benefitMultiplier),
                                  max(df$benefitMultiplier)),10),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b=40),
        showlegend = F
      )
    p
  }
}
save_image(p,paste0(exportPath,"SIPContourPostUtilMitigate",nb1,nb2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPContourPostUtilMitigate",nb1,nb2,".pdf"))

####### Fig 18:(w/ risky asset) Initial Benefit and SIP Preferences #######
# dimensions as percentage of page
w <- .7    #width
h <- .25  #height

# not adjustable parameters
nb1 <- 100
age1 <- 65
benefit1 <- 100
nb2 <- 100
age2 <- 70
benefit2 <- 500

alpha <- 2 #risk aversion level

b1=10 #mortality param
b2=10
# contour plot Fig 18
{
  
  asset1 <- as.vector(benefit1*annuity(age1,.02))
  asset2 <- as.vector(benefit2*annuity(age2,.02))
  
  percBenefit1 <- benefit1/asset1
  percBenefit2 <- benefit2/asset2
  
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste0("simulatedData/BaseRisk_risky_", 
                                  age1,"_",nb1,".rds"))
  
  # import smoothed stability surface
  name <- paste0("simulatedData/smoothedSIP_risky_nb",nb1,nb2,"_b",b1,b2,".rds")
  riskStability <- readRDS(name)
  
  # import smoothed stability when group 2 is on its own
  name <- paste0("simulatedData/smoothedSIP_risky_nb",0,nb2,"_b",b1,b2,".rds")
  riskStabilitySmallHomo <- readRDS(name)
  
  
  age2Vec <- seq(55, 75, by = .1)
  benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
  
  SIP1 <- riskSmallHomo
  SIP2 <- riskStabilitySmallHomo[age2==age2Vec,1]
  
  
  utile1 <- utilityFn(SIP1,percBenefit1, alpha)
  SIPIndif1 <- seq(1,30,.5)
  indifferenceCurve1 <- utilityCurve(utile1,SIPIndif1, alpha)
  
  utile2 <- utilityFn(SIP2,percBenefit2, alpha)
  SIPIndif2 <- seq(1,30,.5)
  indifferenceCurve2 <- utilityCurve(utile2,SIPIndif2, alpha)
  
  baseSplit <- (asset1*nb1)/(asset1*nb1+asset2*nb2)
  baseCase <- meanVarMatrix(nb1 = nb1, age1 = age1,
                            asset1 = asset1,
                            nb2 = nb2, age2 = age2,
                            asset2 = asset2,
                            matSIP = riskStability,
                            assetSplit = baseSplit)
  #background indifference
  splits <- seq(0.01,.99,0.001)
  traceVec <- sapply(splits, function(a)  meanVarMatrix(nb1 = nb1, age1 = age1,
                                                        asset1 = asset1,
                                                        nb2 = nb2, age2 = age2,
                                                        asset2 = asset2,
                                                        matSIP = riskStability,
                                                        assetSplit = a))
  # Create grid
  grid1 <- list(x = traceVec[1,], y = traceVec[2,], splits = splits)
  #fit model
  gam_fit1 <- gam(x ~ te(y, k=20), data = grid1)
  summary(gam_fit1)
  grid1$smoothyX <- predict(gam_fit1, newdata = grid1)
  # Create grid
  grid2 <- list(x = traceVec[3,], y = traceVec[4,], splits = splits)
  #fit model
  gam_fit2 <- gam(x ~ te(y, k=20), data = grid2)
  summary(gam_fit2)
  grid2$smoothyX <- predict(gam_fit2, newdata = grid2)
  
  #get worst points
  {
    #extend indifference curve (check for all benefit level of SIP curve)
    indifferenceCurve1Extended <- utilityCurve(utile1,grid1$smoothyX,alpha)
    indifferenceCurve2Extended <- utilityCurve(utile2,grid2$smoothyX,alpha)
    
    #compare for all benefit level
    diff1 <- grid1$y - indifferenceCurve1Extended
    diff2 <- grid2$y - indifferenceCurve2Extended
    
    noGoodIndex1 <- which(diff1<=0)
    noGoodIndex2 <- which(diff2<=0)
    
  }
  
  goodRebalanceId <- seq_along(splits)[!seq_along(splits)%in%c(noGoodIndex1,
                                                               noGoodIndex2)]
  minY <- (grid2$y[min(goodRebalanceId)]*asset2)/(grid1$y[min(goodRebalanceId)]*asset1)
  maxY <- (grid2$y[max(goodRebalanceId)]*asset2)/(grid1$y[max(goodRebalanceId)]*asset1)
  
  currentY <- benefit2/benefit1
  
  # get worse areas (no better area possible for both at the same time)
  {
    age2Vec <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    
    
    dfworsteSIP <- list(SIP = riskStability[riskStability<=riskStabilitySmallHomo
                                            |riskStability<=riskSmallHomo])
    dfworsteSIP$age2Vec <-  matrix(rep(age2Vec, length(benefitMultiplier)),
                                   ncol = length(benefitMultiplier))[riskStability<=riskStabilitySmallHomo
                                                                     |riskStability<=riskSmallHomo]
    dfworsteSIP$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2Vec)),
                                            nrow = length(age2Vec),
                                            byrow = T)[riskStability<=riskStabilitySmallHomo
                                                       |riskStability<=riskSmallHomo]
    
  }
  
  #table info
  {
    totAsset <- (nb1*asset1+nb2*asset2)
    
    
    benchmarkObjFn1 <- utile1
    benchmarkObjFn2 <- utile2
    SIP1benchmark <- riskSmallHomo
    SIP2benchmark <- riskStabilitySmallHomo[age2Vec ==age2,
                                            which.min(abs(benefitMultiplier-currentY))]
    
    pps1natural <- asset1/benefit1
    pps2natural <- asset2/benefit2
    naturalSIP <- riskStability[age2Vec ==age2,
                                which.min(abs(benefitMultiplier-currentY))]
    naturalObjFn1 <- utilityFn(naturalSIP,1/pps1natural)
    naturalObjFn2 <- utilityFn(naturalSIP,1/pps2natural)
    
    
    lowLimSIP <- riskStability[age2Vec == age2,
                               which.min(abs(benefitMultiplier-minY))]
    lowLimAssetRatio <- as.vector(minY*annuity(age2,.02)/annuity(age1,.02))
    lowLimAsset1 <- totAsset/(nb1+nb2*lowLimAssetRatio)
    lowLimAsset2 <- (totAsset-nb1*lowLimAsset1)/nb2
    lowLimBen1 <- as.vector(lowLimAsset1/annuity(age1, .02))
    lowLimBen2 <- as.vector(lowLimAsset2/annuity(age2, .02))
    pps1lowLim <- asset1/lowLimBen1
    pps2lowLim <- asset2/lowLimBen2
    lowLimObjFn1 <- utilityFn(lowLimSIP,1/pps1lowLim)
    lowLimObjFn2 <- utilityFn(lowLimSIP,1/pps2lowLim)
    
    upLimSIP <- riskStability[age2Vec == age2,
                              which.min(abs(benefitMultiplier-maxY))]
    upLimAssetRatio <- as.vector(maxY*annuity(age2,.02)/annuity(age1,.02))
    upLimAsset1 <- totAsset/(nb1+nb2*upLimAssetRatio)
    upLimAsset2 <- (totAsset-nb1*upLimAsset1)/nb2
    upLimBen1 <- as.vector(upLimAsset1/annuity(age1, .02))
    upLimBen2 <- as.vector(upLimAsset2/annuity(age2, .02))
    pps1upLim <- asset1/upLimBen1
    pps2upLim <- asset2/upLimBen2
    upLimObjFn1 <- utilityFn(upLimSIP,1/pps1upLim)
    upLimObjFn2 <- utilityFn(upLimSIP,1/pps2upLim)
  
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SIP")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2Vec <- age2Vec[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2Vec,
      y = ~benefitMultiplier,
      z = ~SIP,
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
      add_trace(
        data = dfworsteSIP,
        x = ~age2Vec,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(255, 99, 71, 0.2)", size = 6, symbol = "circle"),
        name = "No No Region",
        inherit = FALSE
      ) %>%
      add_trace(
        x = c(age2),
        y = c(currentY),
        type = "scatter",
        mode = "markers",
        marker = list(color = "blue", size = 5, symbol = "circle"),
        name = "Current point",
        inherit = FALSE
      ) %>%
      add_trace(
        x = c(age2,age2),
        y = c(minY,maxY),
        type = "scatter",
        mode = "lines",
        line = list(color = "green", width = 5),
        name = "good reallocation",
        inherit = FALSE
      )%>%
      layout(
        font = list(family = fontType),
        plot_bgcolor = "lightgrey",   # uniform background color
        paper_bgcolor = "white",  # outside background
        xaxis = list(title = list(text = "Age of members in Group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2Vec), max(df$age2Vec)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of Group 2",
                                  standoff = 5),
                     type = "log",
                     showgrid = FALSE,
                     range= log(c(min(df$benefitMultiplier),
                                  max(df$benefitMultiplier)),10),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        margin = list(t = 50, b=40),
        showlegend = F
      )
    p
  }
}
save_image(p,paste0(exportPath,"SIPContourRiskyPostUtilReallocation",nb1,nb2,".pdf"),
           width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
browseURL(paste0(exportPath,"SIPContourRiskyPostUtilReallocation",nb1,nb2,".pdf"))


