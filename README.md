# Heterogeneity

Code used for research about heterogeneity in pension pool

## plotGenerator.R

Contain all the plot in the paper. All the sections a separated by comment lines.

The first few section should be run first to make sure the code generating the plots can work properly. Those sections contains the following.

-   The top of the script contains the packages used for the figures generation and export.
-   It is followed ny the design and export variable for the formatting of the plots at export.
-   Then all the functions used for some of the data generated in the plot script are in one bracket.

The plots are separated by comments line and are ordered as presented in the paper. At the to of each plot code, some dimension formatting options are available

All the data contained in the plot comes from the pre simulated data contained in the 'simulatedData' folder.

## functionsYoS.R

Contain the functions used to estimate the number of years of stability with simulations and a quick example.

It also contain the function to approximate the number of stable years based on the approximation presented in section 11.2, it is much quicker to compute but overestimate the number of stable years.

Some code to create a surface plot is also presented
