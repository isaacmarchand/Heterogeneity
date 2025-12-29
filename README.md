# Heterogeneity

Code used for research about heterogeneity in pension pool

## figureGenerator.R

Contain all the figures in the paper.

####Generating all the figures from the paper

The script can be run as a whole and all Figures should be exported in a specified directory. This export directory should be specified on line 12, right after the list of package required. A working version of Python3 is also required and the location of your python should be specified on line 15. If you don't already have python installed, the commented lines 16-20 can be used to installe it directly from R. 

####Generating select figures

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
-   After hte brackets, some code to export the figures if needed

Almost all the data contained in the figures of section 3 and 5 comes from the pre-simulated data contained in the 'simulatedData' folder. The data for the figures of section 4 is generated directly in the script.

## simulatedData

Folder containing all the data required to generate the figures from the paper. This data was generated using the functions in 'SIP.R'. Those SIP were pre calculated as they are estimated using simulations which can take some time to compute.

## SIP.R

## approxSIP.R

