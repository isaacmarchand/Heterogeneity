library(plotly)
library(reshape2)
###### Base functions #####
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
  
  # Get variance of 1st period adjustment as proxy
  Variance1Periode <- function(age1=65, asset1=100000, nb1=100,
                               age2=65, asset2=200000, nb2=100, rate = .02){
    p1 <- tpx(1,age1)
    p2 <- tpx(1,age2)
    B01 <- as.vector(asset1/annuity(age1, .02))
    B02 <- as.vector(asset2/annuity(age2, .02))
    y <- B02/B01
    
    pN1_0 <- dbinom(0,nb1,p1) 
    pN2_0 <- dbinom(0,nb2,p2) 
    
    mu <- (nb1*p1+y*nb2*p2)/((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))
    sigmaSq <- ((nb1*p1*(1-p1)+nb1^2*p1^2+2*y*nb1*p1*nb2*p2+y^2*nb2*p2*(1-p2)
                 +y^2*nb2^2*p2^2)/((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))) - mu^2
    
    varTerm1 <- ((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))*(1/mu^2+(3*sigmaSq/mu^4))
    varTerm2 <- ((1-pN1_0)+(1-pN2_0)-(1-pN1_0)*(1-pN2_0))^2*(1/mu+(sigmaSq/mu^3))^2
    
    (nb1*p1+y*nb2*p2)^2*(varTerm1-varTerm2)
  }

}

###############################################################################
# Plots for Section 5 SIP (or HSAFE)
###############################################################################

####### Contour Plot Group1's Perspective ########

# adjustable parameters
nb2 <- 100      # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
age1 <- 65      #only 60, 65 and 70 available

nb1 <- 100      # Don't change, but some scenario available at nb1 = (10 and 500)

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
      )
    
    p
  }
  
}



####### Contour Plot Both Groups' Perspective ########

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
        data = dfworsteYoS,
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
}



####### 2D Plot Age Heterogeneity ########

# adjustable parameters
nb2 <- 100      # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
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
    colors <- c( "red", "green", "blue", "cyan", "magenta", "gold","black") #can use rgb code instead
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
        name = paste0("2 has ", diffWealth[i], " times banefits of 1")
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
      name = paste(nb1+nb2, " Homogeneous from Group 1", sep = "")
    )
    
    p <- add_trace(
      p,
      x = c(min(age2), max(age2)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste(nb1, " Homogeneous from Group 1", sep = "")
    )
    
    # Add vertical line at age2 = some value (if desired)
    # Example if you want a vertical reference at 65:
    # p <- add_trace(
    #   p,
    #   x = c(65,65),
    #   y = c(min(slices), max(slices)),
    #   type = 'scatter',
    #   mode = 'lines',
    #   line = list(dash = 'dot', color = 'gray'),
    #   name = 'Reference Age'
    # )
    
    # Final layout
    p <- layout(
      p,
      xaxis = list(title = "age of Group 2"),
      yaxis = list(title = "Nb of Stable Years"),
      legend = list(x = 0.8, y = 1),
      title = paste("Age Groupe 1: ", age1, " / Nb Groupe 1: ", nb1, " / Nb Groupe 2: ", nb2, sep = "")
    )
    
    p
  }
  
}



####### 2D Plot Wealth Heterogeneity ########

# adjustable parameters
nb2 <- 10      # 50, 100 and 200 available for all age1. Also, 5, 10, 500, 1000 available for age1=65
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
    colors <- c( "red", "green", "blue", "cyan", "magenta", "gold","black") #can use rgb code instead
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
        name = paste0("Age 2:", age2[i])
      )
      
      # Add vertical line at Benefit = some value (if desired)
      # Example if you want a vertical reference at max of ages:
      p <- add_trace(
        p,
        x = c(benefitMultiplier[which(max(slices[i,])==slices[i,])],
              benefitMultiplier[which(max(slices[i,])==slices[i,])]),
        y = c(min(slices,riskSmallHomo), max(slices,homoYoS)),
        type = 'scatter',
        mode = 'lines',
        line = list(dash = 'dot', color = colors[i]),
        name = paste("Max stability age", age2[i])
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
      name = paste(nb1+nb2, " Homogeneous from Group 1", sep = "")
    )
    
    p <- add_trace(
      p,
      x = c(min(benefitMultiplier), max(benefitMultiplier)),
      y = c(riskSmallHomo, riskSmallHomo),
      type = 'scatter',
      mode = 'lines',
      line = list(dash = 'dash', color = 'black'),
      name = paste(nb1, " Homogeneous from Group 1", sep = "")
    )
    
    # Final layout
    p <- layout(
      p,
      xaxis = list(title = "Benefit Mutiplier of Group 2",
                   type = "log"),
      yaxis = list(title = "Nb of Stable Years"),
      legend = list(x = 0.8, y = 1),
      title = paste("Age Groupe 1: ", age1, " / Nb Groupe 1: ", nb1, " / Nb Groupe 2: ", nb2, sep = "")
    )
    
    p
  }
  
}




###############################################################################
# Plots for Section 4 SD
###############################################################################