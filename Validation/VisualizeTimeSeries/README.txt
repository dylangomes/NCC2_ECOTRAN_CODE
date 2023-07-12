MatlabOutputToCSV.m pulls in any .mat ECOTRAN output files, and converts them to csv, based on the region selected of interest (see code annotation).

PlotStability.R reads in a 150 year average CUTI run "NCC2_09032022_001_1_0003_28-Mar-2023_SR1.csv" and creates the the 150-year average CUTI stability plot "Equilibrium_150Year.pdf"

TimeSeriesPlotting.R reads in a simulation run with real CUTI values (not averaged) "NCC2_09032022_001_1_0004_28-Mar-2023_SR1.csv" and produces validation plots found in "Figures/"

Fratio_PP_check.m is code that was used to look at calculated F-ratios of primary production and absolute primary production values. Used in tuning of detritus recycling parameters.

QB.csv are consumption to biomass ratios of functional groups - used in calculations.