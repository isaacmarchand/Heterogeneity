# Heterogeneity
Code used for research about heterogeneity in pension pool

## functionsYoS.R
Contain the functions used to estimate the number of years of stability with simulations and a quick example. 

It also contain the function to approximate the number of stable years based on the approximation presented in section 11.2, it is much quicker to compute but overestimate the number of stable years.

Some code to create a surface plot is also presented

## functionSimulatedPVDist.R
Contain the function used to simulate an emperical distribution of the PV of the benefits. It has a function for homogeneous pool and one for heterogeneous pool with 2 subgroup. It also contains example code to plot unconditional distribution, ditribution conditional on the age of death and to plot the the mean and sd of the present value per age.

## functionAsymptoticPVDist.R
Contain a function that compute the variance of the asymptotic distribution of the PV of the benefits for an homogeneous pool.
