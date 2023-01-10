CUTI_daily.nc is actual CUTI values.
CUTI_MEAN.nc is mean CUTI values (mean over course of entire timeseries) repeated for every day. The primary productivity and biomass values were about 10 times higher than expected so CUTI_MEAN_Div10.nc is a file with all these values divided by 10.

CUTI_Visualization.pdf shows 1988 - 2021 CUTI data with an average (red) line on top of it. 

ECOTRANdynamic_NCC2_10072022_batch.m matlab script uses an input ERD_CUTI OR ERD_CONST (constant). CUTI calls the real CUTI file "CUTI_daily.nc" within the script "f_ECOTRANphysics_NCC2_upwelling_09042022.m" ERD_CONST needs to be changed within "f_ECOTRANphysics_NCC2_upwelling_09042022" to the constant file that you'd like to use.