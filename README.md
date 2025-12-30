# Heterogeneity

Code used for article ...

## figureGenerator.R

Contain all the figures in the paper.

### Generating all the figures from the paper

The script can be run as a whole and all Figures should be exported in a specified directory. This export directory should be specified on line 12, right after the list of package required. A working version of Python3 is also required and the location of your python should be specified on line 15. If you don't already have python installed, the commented lines 16-20 can be used to installe it directly from R.

### Generating select figures

The script can also be run manually to generate specific figure with or without exporting them. All the sections a separated by comment lines.

The first few section should be run first to make sure the code generating the figures can work properly. Those sections contains the following.

-   The 'Package' section contains the packages used for the figures generation and export.
-   The 'Export' section contains the Export path and python3 path
-   The 'Export design choice' section contains some general formatting for all figures and the resolution wanted for a full page figure as a baseline.
-   The 'Base functions' section contains all the functions used for the data generated directly in the figure script.

The Following sections separated by comments line are all the figures from the paper in order. All Figures section have the same structure:

-   Dimension of the figure as a percentage of a page for the width 'w' and height 'h'. (for the export only)
-   Some adjustable parameters (for some of the figures)
-   In bracket, all the code required to generate the figure
-   After the brackets, some code to export the figures if needed

Almost all the data contained in the figures of section 3 and 5 comes from the pre-simulated data contained in the 'simulatedData' folder. The data for the figures of section 4 is generated directly in the script.

## simulatedData

Folder containing all the data required to generate the figures from the paper. This data was generated using the functions in 'SIP.R'. Those SIP were pre calculated as they are estimated using simulations which can take some time to compute.

## SIP.R

Contains a function to estimate the SIP using simulations with the control variate estimate introduced in Section 4.2 of the paper. It works for Homogeneous pools and Heterogeneous pool with 2 groups. It only works for a SIP with no upper bound on stability $\epsilon_2 = \infty$.

The first section 'Create functions' contains all the function used to estimate the SIP, with SIP2Pop_CV() being the function used to estimate the SIP for a given pool. You can give to SIP2Pop_CV()
-   (mandatory) The number of member from group 1 and 2 (nb1 and nb2)
-   (mandatory) The age of the members from group 1 and 2 (age1 and age2)
-   (mandatory) The initial benefit of the members from group 1 and 2 (benefit1 and benefit2)
-   (optional) The lower stability bound $\epsilon_1$ (epsilon)
-   (optional) The stability confidence $\beta$ (beta)
-   (optional) The number of simulations used to estimate the SIP (nbSimul) (recommend 10,000)
-   (optional) The Gompertz law modal and dispersion parameters for group 1 and 2 (m1, b1, m2, b2)

It also contains a quick example of how the function SIP2Pop_CV() can be used to compute the SIP for different level of wealth heterogeneity through benefit ratio 'y'.

## SIP_risky.R

Contains a function to estimate the SIP using simulations (CV estimate doesn't allow for risky investment for now). It works for Homogeneous pools and Heterogeneous pool with 2 groups. It only works for a SIP with no upper bound on stability $\epsilon_2 = \infty$.

The first section 'Create functions' contains all the function used to estimate the SIP, with SIP2Pop_risky() being the function used to estimate the SIP for a given pool and risky asset characteristics. You can give to SIP2Pop_risky()
-   (mandatory) The number of member from group 1 and 2 (nb1 and nb2)
-   (mandatory) The age of the members from group 1 and 2 (age1 and age2)
-   (mandatory) The initial benefit of the members from group 1 and 2 (benefit1 and benefit2)
-   (optional) The lower stability bound $\epsilon_1$ (epsilon)
-   (optional) The stability confidence $\beta$ (beta)
-   (optional) The number of simulations used to estimate the SIP (nbSimul) (recommend 100,000)
-   (optional) The characteristic of the balanced portfolio, always (50 risky/50 risk-free), average return and volatility of risky asset and return of risk-free asset (avgRisky, volRisky, and rfrate)
-   (optional) The Gompertz law modal and dispersion parameters for group 1 and 2 (m1, b1, m2, b2)

It also contains a quick example of how the function SIP2Pop_risky() can be used to compute the SIP for different level of wealth heterogeneity through benefit ratio 'y'.

## approxSIP.R

Contains a function to approximate the SIP using the approximation presented in section 4.1 of paper. It works for Homogeneous pools and Heterogeneous pool with 2 groups. It only works for an approximate SIP with no upper bound on stability $\epsilon_2 = \infty$.

The first section 'Create functions' contains all the function used to approximate the SIP, with approxSIP2Pop() being the function used to approximate the SIP for a given pool. You can give to approxSIP2Pop()
-   (optional) The number of member from group 1 and 2 (nb1 and nb2)
-   (optional) The age of the members from group 1 and 2 (age1 and age2)
-   (optional) The initial benefit of the members from group 1 and 2 (ben1 and ben2)
-   (optional) The lower stability bound $\epsilon_1$ (epsilon)
-   (optional) The stability confidence $\beta$ (beta)
-   (optional) The Gompertz law modal and dispersion parameters for group 1 and 2 (m1, b1, m2, b2)

It also contains a quick example of how the function approxSIP2Pop() can be used to compute a surface of approximate SIP for varying wealth and age of the second group, it is then presented in a 3d plot.
