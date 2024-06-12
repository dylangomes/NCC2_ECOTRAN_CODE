function GetEnergyBudget(setWD, Model_name)


% step 1b: define food web model to use -----------------------------------
readFile            	= [setWD Model_name];

% step 1c: load ECOPATH (EwE) model from Aydin VisualBasic file (.csv format)
dat                  	= f_readEwEcsv_10pp_07072021(readFile);	% use for models with up to 10 primary producers

% step 1d: aggregate model results & prep EwEResult for analysis ----------
[EwEResult, PEDIGREE] 	= f_AggregateBiologicalModel_02052021(dat);



% *************************************************************************
% STEP 2: ECOTRAN conversion-----------------------------------------------
MonteCarloStore         = [];
[ECOTRAN]            	= ECOTRANheart_09032021(EwEResult, MonteCarloStore);
% *************************************************************************

writematrix(ECOTRAN.EnergyBudget,strcat(setWD,regexprep(Model_name,".csv",""),"_EnergyBudget.csv"))

