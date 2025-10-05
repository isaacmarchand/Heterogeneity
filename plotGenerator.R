library(plotly)
library(reshape2)
library(RColorBrewer)
library(reticulate)
library(ggplot2)
library(htmlwidgets)
library(webshot)
reticulate::use_python("/opt/anaconda3/bin/python3")

# ####if need to install python
# install.packages("reticulate")
# library(reticulate)
# use_python(install_python())
# py_install(c("kaleido==0.2.1", "plotly"))

###### Design and export choice ####
exportPath <- "/Users/macbook/Library/Mobile\ Documents/com~apple~CloudDocs/School/SFU/Research/Coding/Plots/October10th/"
fontType <- 'Times New Roman' #'Verdana' for report 'Times New Roman' for paper
axisFont <- list(size=18, family = fontType)
titleFont <- list(size=30, family = fontType)
legendFont <- list(size=12, family = fontType)

pixelsFullWidth <- 1240*(5.5/8) #removing margins
pixelsFullHeight <- 1754*(9.3/11) #removing margins
###### Base functions to run first #####
{
  ######Colors
  rgbSOA <- matrix(c(2,77,124,186,191,51,119,196,213,253,206,7,210,49,56,1,1,1, 255,255,255), byrow = TRUE, ncol = 3)/255
  winterColormap <- c(rgb(0,(0:256)/256,(1-((0:256)/256))*.5+.5))
  
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
  
  # Get variance of 1st period adjustment as proxy
  SD1Periode <- function(age1=65, B01=10000, nb1=100,
                         age2=65, B02=20000, nb2=100, rate = .02){
    p1 <- tpx(1,age1)
    p2 <- tpx(1,age2)
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
}

###############################################################################
# Plots for Section 5 SIP
###############################################################################

####### SIP Contour Plot Group1's Perspective Risky Asset ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height

# adjustable parameters
nb2 <- 100     # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available

nb1 <- 100      # Don't recommend to change, but some scenario available at nb1 = (10 and 500)

# Generate plot
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste0("simulatedData/SIPBase_riskyAsset",
                                 nb1,".rds"))
  
  # import smoothed stability surface
  diffAge <- ifelse(age1!=65,paste(age1,"_",sep=""),"")
  name <- paste("simulatedData/controlledSIP_ParrallelComputingRiskyAsset",
                diffAge,nb1,nb2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # compute better and worse areas
  {
    ## get stability when group 1 and 2 are homogeneous from imported surface
    benMultiToExtract <- 1
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoYoS <- riskStability[rowAge,colMultiplier]
    
    ## extract better (green) area
    dfbetterYoS <- list(YoS = riskStability[riskStability>=homoYoS])
    dfbetterYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  >=homoYoS]
    dfbetterYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability>=homoYoS]
    
    ## extract worst (red) area
    dfworsteYoS <- list(YoS = riskStability[riskStability<=riskSmallHomo])
    dfworsteYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  <=riskSmallHomo]
    dfworsteYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
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
      add_trace(
        data = dfbetterYoS,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
        name = "Preferred Region",
        inherit = FALSE
      )%>%
      add_trace(
        data = dfworsteYoS,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  
  save_image(p,paste0(exportPath,"SIPContourRiskyAsset",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPContourRiskyAsset",nb1,nb2,".pdf"))
  
  # # other export option
  # saveWidget(p, paste0(exportPath,"SIPContourRiskyAsset.html"), selfcontained = TRUE)
  # webshot(paste0(exportPath,"SIPContourRiskyAsset.html"),
  #         file   = paste0(exportPath,"SIPContourRiskyAsset.jpeg"),
  #         vwidth = 720, vheight = 480,  # canvas size
  #         zoom   = 3)
}
####### SIP Contour Plot Group1's Perspective ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height

# adjustable parameters
nb2 <- 100     # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available

nb1 <- 100      # Don't recommend to change, but some scenario available at nb1 = (10 and 500)

# Generate plot
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/BaseRisk",
                                 age1,"_",nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  diffAge <- ifelse(age1!=65,paste(age1,"_",sep=""),"")
  name <- paste("simulatedData/controlledYoS_ParrallelComputingBenefits",
                diffAge,nb1,nb2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # compute better and worse areas
  {
    ## get stability when group 1 and 2 are homogeneous from imported surface
    benMultiToExtract <- 1
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoYoS <- riskStability[rowAge,colMultiplier]
    
    ## extract better (green) area
    dfbetterYoS <- list(YoS = riskStability[riskStability>=homoYoS])
    dfbetterYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  >=homoYoS]
    dfbetterYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability>=homoYoS]
    
    ## extract worst (red) area
    dfworsteYoS <- list(YoS = riskStability[riskStability<=riskSmallHomo])
    dfworsteYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  <=riskSmallHomo]
    dfworsteYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
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
      add_trace(
        data = dfbetterYoS,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
        name = "Preferred Region",
        inherit = FALSE
      )%>%
      add_trace(
        data = dfworsteYoS,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  
  save_image(p,paste0(exportPath,"SIPContour1Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPContour1Perspective",nb1,nb2,".pdf"))
}



####### SIP Contour Plot Both Groups' Perspective ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height

# adjustable parameters
nb2 <- 100       # 50, 100 and 200 available for all
age1 <- 65      #only 60, 65 and 70 available

# Generate plot
{
  nb1 <- 100      # Don't change
  
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/BaseRisk",
                                 age1,"_",nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  diffAge <- ifelse(age1!=65,paste(age1,"_",sep=""),"")
  name <- paste("simulatedData/controlledYoS_ParrallelComputingBenefits",
                diffAge,nb1,nb2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # import smoothed stability when group 2 is on its own
  name <- paste("simulatedData/controlledYoS_ParrallelComputingBenefits",0,
                nb2,".rds", sep = "")
  riskStabilitySmallHomo <- readRDS(name)
  
  # get worse areas (no better area possible for both at the same time)
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
    
    
    dfworsteYoS <- list(YoS = riskStability[riskStability<=riskStabilitySmallHomo
                                            |riskStability<=riskSmallHomo])
    dfworsteYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability<=riskStabilitySmallHomo
                                                                  |riskStability<=riskSmallHomo]
    dfworsteYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=riskStabilitySmallHomo
                                                       |riskStability<=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
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
      add_trace(
        data = dfworsteYoS,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  }
  save_image(p,paste0(exportPath,"SIPContour2Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPContour2Perspective",nb1,nb2,".pdf"))
}



####### SIP Contour Plot Both Groups' Perspective Risky Asset########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height

# adjustable parameters
nb2 <- 100       # 50, 100 and 200 available for all
age1 <- 65      #only 60, 65 and 70 available

# Generate plot
{
  nb1 <- 100      # Don't change
  
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/SIPBase_riskyAsset",
                                 nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  diffAge <- ifelse(age1!=65,paste(age1,"_",sep=""),"")
  name <- paste("simulatedData/controlledSIP_ParrallelComputingRiskyAsset",
                diffAge,nb1,nb2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # import smoothed stability when group 2 is on its own
  name <- paste("simulatedData/controlledSIP_ParrallelComputingRiskyAsset",0,
                nb2,".rds", sep = "")
  riskStabilitySmallHomo <- readRDS(name)
  
  # get worse areas (no better area possible for both at the same time)
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
    
    
    dfworsteYoS <- list(YoS = riskStability[riskStability<=riskStabilitySmallHomo
                                            |riskStability<=riskSmallHomo])
    dfworsteYoS$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability<=riskStabilitySmallHomo
                                                                  |riskStability<=riskSmallHomo]
    dfworsteYoS$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=riskStabilitySmallHomo
                                                       |riskStability<=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
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
      add_trace(
        data = dfworsteYoS,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  save_image(p,paste0(exportPath,"SIPContourRiskyAsset2Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPContourRiskyAsset2Perspective",nb1,nb2,".pdf"))
  
  # # other export option
  # saveWidget(p, paste0(exportPath,"SIPContourRiskyAsset2Perspective.html"), selfcontained = TRUE)
  # webshot(paste0(exportPath,"SIPContourRiskyAsset2Perspective.html"),
  #         file   = paste0(exportPath,"SIPContourRiskyAsset2Perspective.jpeg"),
  #         vwidth = 720, vheight = 480,  # canvas size
  #         zoom   = 3)
}



####### SIP 2D Plot Age Heterogeneity ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height

# adjustable parameters
nb2 <- 100       # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available
diffWealth <- c(.2, .5, 1, 2, 5) # Any amount of value from interval [.1,10], Benefit of group 2 compared to group 1

nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/BaseRisk",
                                 age1,"_",nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  diffAge <- ifelse(age1!=65,paste(age1,"_",sep=""),"")
  name <- paste("simulatedData/controlledYoS_ParrallelComputingBenefits",
                diffAge,nb1,nb2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # extract wealth slice
  {
    # get wealth slices
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
    colMultiplier <- sapply(diffWealth, function(x)which.min(abs(benefitMultiplier-x)))
    slices <- riskStability[,colMultiplier]
    
    # get homogeneous value 
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoYoS <- riskStability[rowAge,colMultiplier]
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
      y = c(homoYoS, homoYoS),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'gray'),
      name = paste(nb1+nb2, " members from group 1")
    )
    
    p <- add_trace(
      p,
      x = c(min(age2), max(age2)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste(nb1, " members from group 1")
    )
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in group 2",
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
  save_image(p,paste0(exportPath,"SIPMortalityHeteWealth.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPMortalityHeteWealth.pdf"))
}



####### SIP 2D Plot Wealth Heterogeneity ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb2 <- 100      # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available
age2 <- c(60,65,70) # Any amount of value from interval [55,75]

nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  # import base stability when group 1 is on its own
  riskSmallHomo <- readRDS(paste("simulatedData/BaseRisk",
                                 age1,"_",nb1,".rds", sep = ""))
  
  # import smoothed stability surface
  diffAge <- ifelse(age1!=65,paste(age1,"_",sep=""),"")
  name <- paste("simulatedData/controlledYoS_ParrallelComputingBenefits",
                diffAge,nb1,nb2,".rds", sep = "")
  riskStability <- readRDS(name)
  
  # extract age slice
  {
    # get age slices
    age2Vec <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 100))
    rowAge <- sapply(age2, function(x)which.min(abs(age2Vec-x)))
    slices <- riskStability[rowAge,]
    
    # get homogeneous value 
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2Vec)
    homoYoS <- riskStability[rowAge,colMultiplier]
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
        y = c(min(slices,riskSmallHomo)-1, max(slices,homoYoS)+1),
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
      y = c(homoYoS, homoYoS),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'gray'),
      name = paste(nb1+nb2, " members from group 1"),
      showlegend=T
    )
    
    p <- add_trace(
      p,
      x = c(min(benefitMultiplier), max(benefitMultiplier)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste(nb1, " members from group 1"),
      showlegend=T
    )
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in group 2",
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
  save_image(p,paste0(exportPath,"SIPWealthHeteMortality.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPWealthHeteMortality.pdf"))
}




####### SIP Homogeneous evolution with nb1 ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 <- seq(1,5000)     # only possibility
age1 <- c(60,65,70)   # only 60, 65 and 70 available, can select some are all of them

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # import  stability surface
  {
    riskStability <- matrix(0, length(age1), length(nb1))
    for (i in seq_along(age1)) {
      name <- paste("simulatedData/controlledYoS_ParrallelComputingHomo",age1[i],".rds", sep = "")
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
        name = paste0("Age:", age1[i])
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
  save_image(p,paste0(exportPath,"homoSIP_Smooth.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"homoSIP_Smooth.pdf"))
}



####### SIP Heterogeneous Mortality evolution with nb2 ########
#dimensions as percentage of page
w <- 1    #width
h <- .4   #height


# adjustable parameters
nb1 = 100
age1 = 65
nb2 <- seq(0,200)     # only option
age2 <- c(55,60,65,70,75)   # only 55, 60, 65, 70 and 75 available

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(age2), length(nb2))
    for (i in seq_along(age2)) {
      name <- paste("simulatedData/controlledYoS_ParrallelComputingHeteAge",age2[i],".rds", sep = "")
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
        name = paste0("Age:", age2[i])
      )
    }
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in group 2",
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
  save_image(p,paste0(exportPath,"SIPMortalityHeteNb.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPMortalityHeteNb.pdf"))
}



####### SIP Heterogeneous Wealth evolution with nb2 ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = 100
age1 = 65
nb2 <- seq(0,200)     
BMulti <- c(.2,.5,1,2,5) #ratio of benefit2/benefit1

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(BMulti), length(nb2))
    for (i in seq_along(BMulti)) {
      name <- paste("simulatedData/controlledYoS_ParrallelComputingHeteWealth",BMulti[i],".rds", sep = "")
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
      xaxis = list(title = list(text = "Number of members in group 2",
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
  save_image(p,paste0(exportPath,"SIPWealthHeteNb.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPWealthHeteNb.pdf"))
}




####### SIP Heterogeneous Wealth evolution with TOTAL participants ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = seq(0,100)
BMulti <- c(.1,.2,1,5,10) #ratio of benefit2/benefit1

age = 65
nb2 <- rev(nb1)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # Import Data
  {
    riskStability <- matrix(0, length(BMulti), length(nb2))
    for (i in seq_along(BMulti)) {
      name <- paste("simulatedData/ControlledSIP_ParrallelComputingHeteWealth",BMulti[i],"Tot.rds", sep = "")
      riskStability[i,] <- readRDS(name)
    }
  }
  
  #plot
  {
    colors <- rgb(rgbSOA[,1],rgbSOA[,2],rgbSOA[,3]) #can use rgb code instead
    colors <- rep(colors, length.out = nrow(riskStability))  # ensure enough colors
    
    p <- plot_ly() 
    
    # Add each column of slices as a separate trace
    for (i in seq_along(BMulti)) {
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
      xaxis = list(title = list(text = "Number of members in group 2",
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
  save_image(p,paste0(exportPath,"SIPWealthHeteTot.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SIPWealthHeteTot.pdf"))
  
  #other export option
  # saveWidget(p, paste0(exportPath,"SIPWealthHeteTot.html"), selfcontained = TRUE)
  # webshot(paste0(exportPath,"SIPWealthHeteTot.html"),
  #         file   = paste0(exportPath,"SIPWealthHeteTot.jpeg"),
  #         vwidth = 720, vheight = 480,  # canvas size
  #         zoom   = 3)
}
####### SIP Heterogeneous Mortality evolution with TOTAL participants ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = seq(0,100)
ages <- c(55,60,65,70,75) 
nb2 <- rev(nb1)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  # Import Data
  {
    riskStability <- matrix(0, length(BMulti), length(nb2))
    for (i in seq_along(BMulti)) {
      name <- paste("simulatedData/ControlledSIP_ParrallelComputingHeteMortality",ages[i],"Tot.rds", sep = "")
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
        name = paste0("Age ", ages[i])
      )
    }
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"SDMortalityHeteTot.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDMortalityHeteTot.pdf"))
  
  #other export option
  # saveWidget(p, paste0(exportPath,"SDMortalityHeteTot.html"), selfcontained = TRUE)
  # webshot(paste0(exportPath,"SDMortalityHeteTot.html"),
  #         file   = paste0(exportPath,"SDMortalityHeteTot.jpeg"),
  #         vwidth = 720, vheight = 480,  # canvas size
  #         zoom   = 3)
}



#########################################
# Section 5: other figures then results
#########################################
####### Empirical Dist of NB of stable years with smoothing######
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb2 <- 50     #can be anything, computed directly in code
nb1 <- 100     
nbSimul <- 10000  #can add more but <5 sec to run for 10000 simul
beta <- .95       #treshhold illustrated in plot

##preparing data
{
  set.seed(1)
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
    discreteVaR <- min(which(cumulDist(1:30)>(1-beta)))
    
    p <- (1-beta)
    x <- sort(Stability)
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
      add_trace(x = x_segments, y = y_segments, type = 'scatter', mode = 'lines',
                line = list(color = "#010101", shape = "hv"),
                name = "ECDF")%>%
      # Add points at the jumps
      add_trace(x = x_jumps, y = y_jumps,
                type = 'scatter', mode = 'markers',
                marker = list(color = "#010101", size = 6),
                name = "ECDF points",
                showlegend = F)
    
    # Vertical dashed line at discreteVaR
    fig <- fig %>%
      add_trace(x = c(discreteVaR, discreteVaR),
                y = c(cumulDist(discreteVaR), -1),
                type = 'scatter', mode = 'lines',
                line = list(color = "#77C4D5", dash = "dash"),
                name = "5th percentile ECDF"
      )
    
    # Continuous greenish line
    fig <- fig %>%
      add_trace(x = c(y, 34), y = cumul,
                type = 'scatter', mode = 'lines',
                line = list(color = "#BABF33", shape = "linear"),
                name = "Smoothed ECDF")
    
    # Horizontal dashed line
    fig <- fig %>%
      add_trace(x = c(0, max(Stability)), y = c(p, p),
                type = 'scatter', mode = 'lines',
                line = list(color = "#010101", dash = "dash"),
                name = "5th percentile",
                showlegend = F)
    
    # Vertical dashed line at VaR
    fig <- fig %>%
      add_trace(x = c(VaR, VaR), y = c(0.05, -1),
                type = 'scatter', mode = 'lines',
                line = list(color = "#FDCE07", dash = "dash"),
                name = "5th percentile smoothed ECDF"
                )
    fig <- layout(
      fig,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Years",
                                standoff = 5),
                   range = c(4, 14),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Distribution of the number of years"),
                   range = c(0, 0.1),
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
    fig
  } 
  save_image(fig,paste0(exportPath,"smoothedVSempericalSIP.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"smoothedVSempericalSIP.pdf"))
}

###############################################################################
# Plots for Section 4 SD 
# (Adjustable parameter constraint are only suggestion to match section 5,
#  any values are possible since data is generated in the script)
###############################################################################
####### SD Contour Plot Group1's Perspective ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb2 <- 100      # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available

nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)

# Generate plot
{
  # compute base stability when group 1 is on its own
  riskSmallHomo <- SD1Periode(age1, nb1 = nb1, nb2 = 0)
  
  # compute stability surface
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    riskStability <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      riskStability[i,] <- sapply(benefitMultiplier,function(x)SD1Periode(age1,
                                                                          1000,
                                                                          nb1,
                                                                          age2[i],
                                                                          1000*x,
                                                                          nb2))
    }
  }
  
  # compute better and worse areas
  {
    ## get stability when group 1 and 2 are homogeneous from imported surface
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoSD <- riskStability[rowAge,colMultiplier]
    
    ## extract better (green) area
    dfbetterSD <- list(SD = riskStability[riskStability<=homoSD])
    dfbetterSD$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  <=homoSD]
    dfbetterSD$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability<=homoSD]
    
    ## extract worst (red) area
    dfworsteSD <- list(SD = riskStability[riskStability>=riskSmallHomo])
    dfworsteSD$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability
                                                                  >=riskSmallHomo]
    dfworsteSD$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability>=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SD")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SD,
      type = "contour",
      showscale = FALSE,
      contours = list(
        coloring = "lines",  # or "lines", "none"
        showlabels = TRUE
      ),
      line = list(smoothing = 0),
      colorscale = list(c(0, "black"), c(1, "black")),
      reversescale = FALSE
    )%>%
      add_trace(
        data = dfbetterSD,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
        name = "Preferred Region",
        inherit = FALSE
      )%>%
      add_trace(
        data = dfworsteSD,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                              standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  save_image(p,paste0(exportPath,"SDContour1Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDContour1Perspective",nb1,nb2,".pdf"))
}



####### SD Contour Plot Both Groups' Perspective ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb2 <- 100       # 50, 100 and 200 available for all
age1 <- 65      #only 60, 65 and 70 available

# Generate plot
{
  nb1 <- 100      # Don't change, you can but won't match plot of section 5
  
  # compute base stability when group 1 is on its own
  riskSmallHomo <- SD1Periode(age1, nb1 = nb1, nb2 = 0)
  
  # compute stability surface
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    riskStability <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      riskStability[i,] <- sapply(benefitMultiplier,function(x)SD1Periode(age1,
                                                                          1000,
                                                                          nb1,
                                                                          age2[i],
                                                                          1000*x,
                                                                          nb2))
    }
  }
  
  # compute stability when group 2 is on its own
  {
    riskStabilitySmallHomo <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      riskStabilitySmallHomo[i,] <- sapply(benefitMultiplier,
                                           function(x)SD1Periode(age1,1000,0,
                                                                 age2[i],1000*x,
                                                                 nb2))
    }
  }
  
  # get worse areas (no better area possible for both at the same time)
  {
    dfworsteSD <- list(SD = riskStability[riskStability>=riskStabilitySmallHomo
                                            |riskStability>=riskSmallHomo])
    dfworsteSD$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                                ncol = length(benefitMultiplier))[riskStability>=riskStabilitySmallHomo
                                                                  |riskStability>=riskSmallHomo]
    dfworsteSD$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                            nrow = length(age2),
                                            byrow = T)[riskStability>=riskStabilitySmallHomo
                                                       |riskStability>=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(riskStability)
    colnames(df) <- c("ageIndex", "benefitIndex", "SD")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SD,
      type = "contour",
      showscale = FALSE,
      contours = list(
        coloring = "lines",  # or "lines", "none"
        showlabels = TRUE
      ),
      line = list(smoothing = 0),
      colorscale = list(c(0, "black"), c(1, "black")),
      reversescale = FALSE
    )%>%
      add_trace(
        data = dfworsteSD,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  save_image(p,paste0(exportPath,"SDContour2Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDContour2Perspective",nb1,nb2,".pdf"))
}



####### SD 2D Plot Age Heterogeneity ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb2 <- 100      # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available
diffWealth <- c(.2, .5, 1, 2, 5) # Any amount of value from interval [.1,10], Benefit of group 2 compared to group 1

nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  # compute base stability when group 1 is on its own
  riskSmallHomo <- SD1Periode(age1, nb1 = nb1, nb2 = 0)
  
  # compute smoothed stability surface
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    riskStability <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      riskStability[i,] <- sapply(benefitMultiplier,function(x)SD1Periode(age1,
                                                                          1000,
                                                                          nb1,
                                                                          age2[i],
                                                                          1000*x,
                                                                          nb2))
    }
  }
  
  # extract wealth slice
  {
    # get wealth slices
    colMultiplier <- sapply(diffWealth, function(x)which.min(abs(benefitMultiplier-x)))
    slices <- riskStability[,colMultiplier]
    
    # get homogeneous value 
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoSD <- riskStability[rowAge,colMultiplier]
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
      y = c(homoSD, homoSD),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'gray'),
      name = paste0(nb1+nb2, " members from group 1")
    )
    
    p <- add_trace(
      p,
      x = c(min(age2), max(age2)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste0(nb1, " members from group 1")
    )
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"SDMortalityHeteWealth.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDMortalityHeteWealth.pdf"))
}



####### SD 2D Plot Wealth Heterogeneity ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb2 <- 100     # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available
age2 <- c(60,65,70) # Any amount of value from interval [55,75]

nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  # compute base stability when group 1 is on its own
  riskSmallHomo <- SD1Periode(age1, nb1 = nb1, nb2 = 0)
  
  # compute stability surface
  {
    age2Vec <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    riskStability <- matrix(0, length(age2Vec), length(benefitMultiplier))
    for (i in seq_along(age2Vec)) {
      riskStability[i,] <- sapply(benefitMultiplier,function(x)SD1Periode(age1,
                                                                          1000,
                                                                          nb1,
                                                                          age2Vec[i],
                                                                          1000*x,
                                                                          nb2))
    }
  }
  
  # extract age slice
  {
    # get age slices
    age2Vec <- seq(55, 75, by = .1)
    rowAge <- sapply(age2, function(x)which.min(abs(age2Vec-x)))
    slices <- riskStability[rowAge,]
    
    # get homogeneous value 
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2Vec)
    homoSD <- riskStability[rowAge,colMultiplier]
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
        x = c(benefitMultiplier[which(min(slices[i,])==slices[i,])],
              benefitMultiplier[which(min(slices[i,])==slices[i,])]),
        y = c(min(slices,riskSmallHomo)-1, max(slices,homoSD)+1),
        type = 'scatter',
        mode = 'lines',
        line = list(dash = 'dot', color = colors[i]),
        name = paste("Max SD age", age2[i]),
        showlegend=F
      )
    }
    
    # Add horizontal lines
    p <- add_trace(
      p,
      x = c(min(benefitMultiplier), max(benefitMultiplier)),
      y = c(homoSD, homoSD),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'gray'),
      name = paste0(nb1+nb2, " members from group 1"),
      showlegend=T
    )
    
    p <- add_trace(
      p,
      x = c(min(benefitMultiplier), max(benefitMultiplier)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste0(nb1, " members from group 1"),
      showlegend=T
    )
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),showgrid=FALSE,
                   ticks    = "outside",
                   type = "log",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"SDWealthHeteMortality.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDWealthHeteMortality.pdf"))
}



####### SD Homogeneous evolution with nb1 ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 <- seq(1,200)     # Can't contain 0
age1 <- c(60,65,70)   # only 60, 65 and 70 available in section 5

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(age1), length(nb1))
    for (i in seq_along(age1)) {
      riskStability[i,] <- sapply(nb1,function(x)SD1Periode(age1[i], 1000,
                                                            x, nb2 = 0))
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
                   nticks = 4,
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"homoSD_Smooth.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"homoSD_Smooth.pdf"))
}


####### SD Heterogeneous Mortality evolution with nb2 ########
#dimensions as percentage of page
w <- 1   #width
h <- .4  #height


# adjustable parameters
nb1 = 100
age1 = 65
nb2 <- seq(0,500)     
age2 <- c(55, 60,65,70,75)   # only 55, 60, 65, 70 and 75 available in section 5

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(age2), length(nb2))
    for (i in seq_along(age2)) {
      riskStability[i,] <- sapply(nb2,function(x) SD1Periode(age1, 1000, nb1,
                                                             age2 = age2[i],
                                                             B02 = 1000,
                                                             nb2 = x))
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
      xaxis = list(title = list(text = "Number of members in group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"SDMortalityHeteNb.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDMortalityHeteNb.pdf"))
}



####### SD Heterogeneous Wealth evolution with nb2 ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = 100
age1 = 65
nb2 <- seq(0,500)     
BMulti <- c(.2,.5,1,2,5) #ratio of benefit2/benefit1

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(BMulti), length(nb2))
    for (i in seq_along(BMulti)) {
      riskStability[i,] <- sapply(nb2,function(x) SD1Periode(age1, 1000, nb1,
                                                             age2 = age1,
                                                             B02 = 1000*BMulti[i],
                                                             nb2 = x))
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
      xaxis = list(title = list(text = "Number of members in group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"SDWealthHeteNb.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDWealthHeteNb.pdf"))
}




####### SD Heterogeneous Wealth evolution with TOTAL participants ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = seq(0,100)
BMulti <- c(.2,.5,1,2,5) #ratio of benefit2/benefit1

age = 65
nb2 <- rev(nb1)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(BMulti), length(nb1))
    for (i in seq_along(nb1)) {
      riskStability[,i] <- sapply(BMulti,function(x) SD1Periode(age, 1000, nb1[i],
                                                             age2 = age,
                                                             B02 = 1000*x,
                                                             nb2 = nb2[i]))
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
      xaxis = list(title = list(text = "Number of members in group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"SDWealthHeteTot.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDWealthHeteTot.pdf"))
}


####### SD Heterogeneous Mortality evolution with TOTAL participants ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = seq(0,100)
ages <- c(55,60,65,70,75) 
nb2 <- rev(nb1)

# Generate plot (put it in full screen before saving for better placement of legend)
{
  
  # compute  stability surface
  {
    riskStability <- matrix(0, length(BMulti), length(nb1))
    for (i in seq_along(nb1)) {
      riskStability[,i] <- sapply(ages,function(x) SD1Periode(ages[3], 1000, nb1[i],
                                                              age2 = x,
                                                              B02 = 1000,
                                                              nb2 = nb2[i]))
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
        name = paste0("Age ", ages[i])
      )
    }
    
    # Final layout
    p <- layout(
      p,
      font = list(family = fontType),
      xaxis = list(title = list(text = "Number of members in group 2",
                                standoff = 5),
                   titlefont = axisFont,
                   tickfont = list(size = 12),
                   ticks    = "outside",
                   ticklen  = 8,
                   showline = TRUE, mirror = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = "Standard deviation"),
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
  save_image(p,paste0(exportPath,"SDMortalityHeteTot.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SDMortalityHeteTot.pdf"))
}

####### SD 3D Homogeneous Nb People to Age ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
age1 <- seq(55,85)
nb1 <- seq(100,500,25) 

nb2 <- 0 # to stay homogeneous
# Generate plot
{
  # compute stability surface
  {
    riskStability <- matrix(0, length(age1), length(nb1))
    for (i in seq_along(age1)) {
      riskStability[i,] <- sapply(nb1,function(x)SD1Periode(age1[i],
                                                            1000,
                                                            x,
                                                            nb2 = nb2))
    }
  }
  
  # 3D Plot
  {
    p <- plot_ly()%>%add_surface(
      x = ~age1,
      y = ~nb1,
      z = ~t(riskStability),
      type = "surface",
      showscale=F,
      colors = winterColormap)
    # Add grid lines in x-direction
    for (j in seq(1,length(nb1))) {
      p <- p %>% add_trace(
        x = age1,
        y = rep(nb1[j], length(age1)),
        z = riskStability[,j ],
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
    }
    
    # Add grid lines in y-direction
    for (i in seq(1,length(age1))) {
      p <- p %>% add_trace(
        x = rep(age1[i], length(nb1)),
        y = nb1,
        z = riskStability[i,],
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
    }
    p <- p %>%
      layout(
        font = list(family = fontType),
        scene = list(
          xaxis = list(title = "Age",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          yaxis = list(title = "Number of members",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          zaxis = list(title = "Standard deviation",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          aspectratio = list(x = 1, y = 2, z = 1),
          camera = list(
            # Position the camera to look at the plot from a new angle
            eye = list(x = -2.2*.9, y = 1*.9, z = .7*.9),
            # shift image center
            center = list(x=1*.05,y=1.9*.03,z=-.25))
        ),
        margin=list(t=0,l=0,b=0,r=0)
      )
    p
  }
  save_image(p,paste0(exportPath,"SD3dMortalityNB1.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SD3dMortalityNB1.pdf"))
}



####### SD 3D Heterogeneous Mortality evolution with nb2 ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = 100
age1 = 65
age2 <- seq(55,85)
nb2 <- seq(0,500,1) 

# Generate plot
{
  # compute stability surface
  {
    riskStability <- matrix(0, length(age2), length(nb2))
    for (i in seq_along(age2)) {
      riskStability[i,] <- sapply(nb2,function(x)SD1Periode(age1,
                                                            1000,
                                                            nb1,
                                                            nb2 = x,
                                                            age2 = age2[i],
                                                            B02 = 1000))
    }
  }
  
  # 3D Plot
  {
    p <- plot_ly()%>%add_surface(
      x = ~age2,
      y = ~nb2,
      z = ~t(riskStability),
      type = "surface",
      showscale=F,
      colors = winterColormap)
    # Add grid lines in x-direction
    for (j in seq(1,length(nb2), by= 20)) {
      p <- p %>% add_trace(
        x = age2,
        y = rep(nb2[j], length(age2)),
        z = riskStability[,j ],
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
    }
    
    # Add grid lines in y-direction
    for (i in seq(1,length(age2))) {
      p <- p %>% add_trace(
        x = rep(age2[i], length(nb2)),
        y = nb2,
        z = riskStability[i,],
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
    }
    
    #thick line on homogeneous case
    p <- p %>% add_trace(
      x = rep(age1, length(nb2)),
      y = nb2,
      z = riskStability[which(age1==age2),],
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 10),
      showlegend = FALSE
    )
    p <- p %>%
      layout(
        font = list(family = fontType),
        scene = list(
          xaxis = list(title = "Age",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          yaxis = list(title = "Number of members in group 2",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          zaxis = list(title = "Standard deviation",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          aspectratio = list(x = 1, y = 2, z = 1),
          camera = list(
            # Position the camera to look at the plot from a new angle
            eye = list(x = -2.2*.9, y = 1*.9, z = .7*.9),
            # shift image center
            center = list(x=1*.05,y=1.9*.03,z=-.25))
        ),
        margin=list(t=0,l=0,b=0,r=0)
      )
    p
  }  
  save_image(p,paste0(exportPath,"SD3dMortalityNB2.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SD3dMortalityNB2.pdf"))
}

####### SD 3D Heterogeneous Wealth evolution with nb2 ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height


# adjustable parameters
nb1 = 100
age1 = 65
age2 <- 65
y <- exp(seq(log(1/10),log(10), length.out = 21))
nb2 <- seq(0,500,1) 

# Generate plot
{
  # compute stability surface
  {
    riskStability <- matrix(0, length(y), length(nb2))
    for (i in seq_along(y)) {
      riskStability[i,] <- sapply(nb2,function(x)SD1Periode(age1,
                                                            1000,
                                                            nb1,
                                                            nb2 = x,
                                                            age2 = age2,
                                                            B02 = 1000*y[i]))
    }
  }
  
  # 3D Plot
  {
    p <- plot_ly()%>%add_surface(
      x = ~y,
      y = ~nb2,
      z = ~t(riskStability),
      type = "surface",
      showscale=F,
      colors = winterColormap)
    # Add grid lines in x-direction
    for (j in seq(1,length(nb2), by = 20)) {
      p <- p %>% add_trace(
        x = y,
        y = rep(nb2[j], length(y)),
        z = riskStability[,j ],
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
    }
    
    # Add grid lines in y-direction
    for (i in seq(1,length(y))) {
      p <- p %>% add_trace(
        x = rep(y[i], length(nb2)),
        y = nb2,
        z = riskStability[i,],
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
    }
    
    #thick line
    p <- p %>% add_trace(
      x = rep(1, length(nb2)),
      y = nb2,
      z = riskStability[which.min(abs(y-1)),],
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 10),
      showlegend = FALSE
    )
    p <- p %>%
      layout(
        font = list(family = fontType),
        scene = list(
          xaxis = list(title = "Initial benefit of group 2",
                       type = "log",
                       showgrid = F,
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          yaxis = list(title = "Number of members in group 2",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          zaxis = list(title = "Standard deviation",
                       titlefont = axisFont,
                       titlefont = axisFont,
                       tickfont = list(size = 12),
                       ticks    = "outside",
                       ticklen  = 8,
                       showline = TRUE, mirror = TRUE, zeroline = FALSE),
          aspectratio = list(x = 1, y = 2, z = 1),
          camera = list(
            # Position the camera to look at the plot from a new angle
            eye = list(x = 2.2*.85, y = 1.5*.85, z = .7*.85),
            # shift image center
            center = list(x=0,y=0,z=-.3))
        ),
        margin=list(t=0,l=0,b=0,r=0)
      )
    p
  }   
  save_image(p,paste0(exportPath,"SD3dWealthNB2.pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"SD3dWealthNB2.pdf"))
}


#########################################
# Section 4 with implied nb1 scale instead of SD
#########################################
####### Implied Nb1 Contour Plot Group1's Perspective ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height

# adjustable parameters
nb2 <- 100      # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      # only 60, 65 and 70 available

nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)

# Generate plot
{
  # compute base stability when group 1 is on its own
  riskSmallHomo <- SD1Periode(age1, nb1 = nb1, nb2 = 0)
  
  # compute SD surface
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    riskStability <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      riskStability[i,] <- sapply(benefitMultiplier,function(x)SD1Periode(age1,
                                                                          1000,
                                                                          nb1,
                                                                          age2[i],
                                                                          1000*x,
                                                                          nb2))
    }
  }
  
  # compute implied Nb surface
  {
    homoSdGrid <- sapply(1:1000, function(x) SD1Periode(age1 = age1,1000,
                                                        nb1 = x, nb2 = 0))
    impliedNb1 <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      impliedNb1[i,] <- sapply(seq_along(riskStability[i,]),
                               function(j) which.min(abs(riskStability[i,j]-homoSdGrid)))
    }
  }
  
  # compute better and worse areas
  {
    ## get stability when group 1 and 2 are homogeneous from imported surface
    benMultiToExtract <- 1
    colMultiplier <- which.min(abs(benefitMultiplier-benMultiToExtract)) 
    rowAge <- which(age1==age2)
    homoSD <- riskStability[rowAge,colMultiplier]
    
    ## extract better (green) area
    dfbetterSD <- list(SD = riskStability[riskStability<=homoSD])
    dfbetterSD$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                               ncol = length(benefitMultiplier))[riskStability
                                                                 <=homoSD]
    dfbetterSD$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                           nrow = length(age2),
                                           byrow = T)[riskStability<=homoSD]
    
    ## extract worst (red) area
    dfworsteSD <- list(SD = riskStability[riskStability>=riskSmallHomo])
    dfworsteSD$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                               ncol = length(benefitMultiplier))[riskStability
                                                                 >=riskSmallHomo]
    dfworsteSD$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                           nrow = length(age2),
                                           byrow = T)[riskStability>=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(impliedNb1)
    colnames(df) <- c("ageIndex", "benefitIndex", "SD")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SD,
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
        data = dfbetterSD,
        x = ~age2,
        y = ~benefitMultiplier,
        type = "scatter",
        mode = "markers",
        marker = list(color = "rgba(144, 238, 144, 0.3)", size = 6, symbol = "circle"),
        name = "Preferred Region",
        inherit = FALSE
      )%>%
      add_trace(
        data = dfworsteSD,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  save_image(p,paste0(exportPath,"ImpliedSDContour1Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"ImpliedSDContour1Perspective",nb1,nb2,".pdf"))
}



####### Implied Nb1 Contour Plot Both Groups' Perspective ########
#dimensions as percentage of page
w <- 1    #width
h <- .4  #height

# adjustable parameters
nb2 <- 100       # 50, 100 and 200 available for all
age1 <- 65      #only 60, 65 and 70 available

# Generate plot
{
  nb1 <- 100      # Don't change, you can but won't match plot of section 5
  
  # compute base stability when group 1 is on its own
  riskSmallHomo <- SD1Periode(age1, nb1 = nb1, nb2 = 0)
  
  # compute stability surface
  {
    age2 <- seq(55, 75, by = .1)
    benefitMultiplier <- exp(seq(log(1/10),log(10), length.out = 101))
    riskStability <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      riskStability[i,] <- sapply(benefitMultiplier,function(x)SD1Periode(age1,
                                                                          1000,
                                                                          nb1,
                                                                          age2[i],
                                                                          1000*x,
                                                                          nb2))
    }
  }
  
  # compute implied Nb surface
  {
    homoSdGrid <- sapply(1:1000, function(x) SD1Periode(age1 = age1,1000,
                                                        nb1 = x, nb2 = 0))
    impliedNb1 <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      impliedNb1[i,] <- sapply(seq_along(riskStability[i,]),
                               function(j) which.min(abs(riskStability[i,j]-homoSdGrid)))
    }
  }
  
  # compute stability when group 2 is on its own
  {
    riskStabilitySmallHomo <- matrix(0, length(age2), length(benefitMultiplier))
    for (i in seq_along(age2)) {
      riskStabilitySmallHomo[i,] <- sapply(benefitMultiplier,function(x)SD1Periode(age1,
                                                                                   1000,
                                                                                   0,
                                                                                   age2[i],
                                                                                   1000*x,
                                                                                   nb2))
    }
  }
  
  # get worse areas (no better area possible for both at the same time)
  {
    dfworsteSD <- list(SD = riskStability[riskStability>=riskStabilitySmallHomo
                                          |riskStability>=riskSmallHomo])
    dfworsteSD$age2 <-  matrix(rep(age2, length(benefitMultiplier)),
                               ncol = length(benefitMultiplier))[riskStability>=riskStabilitySmallHomo
                                                                 |riskStability>=riskSmallHomo]
    dfworsteSD$benefitMultiplier <- matrix(rep(benefitMultiplier, length(age2)),
                                           nrow = length(age2),
                                           byrow = T)[riskStability>=riskStabilitySmallHomo
                                                      |riskStability>=riskSmallHomo]
    
  }
  
  # contour Plot
  {
    # Prepare data in long format
    df <- melt(impliedNb1)
    colnames(df) <- c("ageIndex", "benefitIndex", "SD")
    df$benefitMultiplier <- benefitMultiplier[df$benefitIndex]
    df$age2 <- age2[df$ageIndex]
    
    # Create contour plot
    p <- plot_ly(
      data = df,
      x = ~age2,
      y = ~benefitMultiplier,
      z = ~SD,
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
        data = dfworsteSD,
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
        xaxis = list(title = list(text = "Age of members in group 2",
                                  standoff = 5),
                     showgrid = FALSE,
                     range = c(min(df$age2), max(df$age2)),
                     titlefont = axisFont,
                     tickfont = list(size = 12),
                     ticks    = "outside",
                     ticklen  = 8,
                     showline = TRUE, mirror = TRUE, zeroline = FALSE
        ),
        yaxis = list(title = list(text = "Initial benefit of group 2",
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
  save_image(p,paste0(exportPath,"ImpliedSDContour2Perspective",nb1,nb2,".pdf"),
             width = w*pixelsFullWidth, height = h*pixelsFullHeight, scale = 1)
  browseURL(paste0(exportPath,"ImpliedSDContour2Perspective",nb1,nb2,".pdf"))
}


