% Calculate production rates

% *************************************************************************
% STEP 1: 
% step 1a: define food web model to use -----------------------------------

% ReadFile_directory      = 'C:/Users/dgome/Documents/NCC_Comparison/post_MHW_Web/';
% BiologicalModel_name	= 'NCC2_NCC1Matched_10042022.csv';

ReadFile_directory      = 'C:/Users/dgome/Documents/NCC_Comparison/pre_MHW_Web/';
BiologicalModel_name	= 'NCC_11242020_JimsOld_updatedKrill_10052022.csv';
% *************************************************************************

% calls:
%       f_readEwEcsv_10pp_08042020                	read in ECOPATH (EwE) model from VisualBasic .csv file and store as variable 'dat'; (use for VisualBasic food web files allowing up to 10 primary producers)
%       f_AggregateBiologicalModel_02052021        	prepare EwE model for use by ECOTRAN; also, aggregate functional groups here if wanted
%           f_calcEE_12292020                          calculate Ecotrophic Efficiency
%           f_VarianceDivision_12132018                calculate the variance of one term divided by another term
%           f_VarianceMultiplication_12132018          calculate the variance of two products
%
%       ECOTRANheart_02062021                       generate an ECOTRAN (E2E) model; the heart of ECOTRAN
%           f_ECOfunction_02062021                     returns a single ECOTRAN model for 1 "type" EwE model or 1 MonteCarlo EwE model
%               f_RedistributeCannibalism_11202019       remove cannibalism terms on diagonal of matrix EwE_diet (f_ECOfunction_12292020 makes adjustments to feces & metabolism terms to account for cannibalism)
%               f_calcEE_12292020                        calculate Ecotrophic Efficiency 
%               f_calcPredationBudget_12102019      	 for each producer, p, and consumer, c: ((b_pc * Q_c) / M2_p)) = the fraction of total predation on each producer, p, going to each consumer, c; (2D matrix: num_grps X num_grps; consumers X producers)
%
%       f_E2Epedigree_08042020                      calculate uncertainty terms for all elements of the ECOTRAN EnergyBudget (A_cp) as Coefficients of Variation (CV)
%           f_VarianceDivision_12132018                calculate the variance of one term divided by another term
%           f_VarianceMultiplication_12132018          calculate the variance of two products
%       f_E2E_MonteCarlo_08042020                   generate a set of randomly generated ECOTRAN models, multiple alternate versions of the EnergyBudget (A_cp); model 1 of the stack is the "type" model defined by the parameters of the VisualBasic .csv file
%
%       f_WebProductivity_03272019                  calculate production rates of each functional group (vertical vector: num_grps X 1)
%
%       f_ScenarioGenerator_10212020                modify EnergyBudget, FFF-->ConsumptionBudget, FFF-->fates, & (QQQ DiscardFraction)
%           f_WebProductivity_03272019                 calculate production rates of each functional group (vertical vector: num_grps X 1)
%       f_CompileScenarioResults_09102020           simple function to compile scenario model results relative to base model conditions
%       p_PlotScenarioResults_07022020
%           myboxplot_5                             modified from the Jorn_Diedrichsen_Toolbox
%           box utility function                    from the Jorn_Diedrichsen_Toolbox
%           round2                                  Round to a specified number of decimals. (From Matlab community user group)

% ***************************************************

readFile            	= [ReadFile_directory BiologicalModel_name];

% step 1b: load ECOPATH (EwE) model from Aydin VisualBasic file (.csv format)
dat                  	= f_readEwEcsv_10pp_07072021(readFile);	% use for models with up to 10 primary producers (use for 10 pp version of Aydin .xlsm, regardless of number of actual pp in model)

% step 1c: aggregate model results & prep EwEResult for analysis ----------
[EwEResult, PEDIGREE] 	= f_AggregateBiologicalModel_02052021(dat);
% *************************************************************************


% *************************************************************************
% STEP 2: load pre-generated Monte Carlo EwE models------------------------
% step 2a: Identify Monte Carlo method ------------------------------------
%          method 1: PREFERRED METHOD --> start with the one original ECOTRAN base model and generate a set of Monte Carlo E2E models of the ECOTRAN production matrix using a predefined CV (or defualt pedigree CV = 0.5)
%          method 2: use pre-generated Monte Carlo EwE models, convert to ECOTRAN E2E models, and run scenarios directly on EACH of these individually
%                    (each matrix being its own base model and scenario model)
MonteCarlo_method           = 1;    % SSS; choose 1 or 2; NOTE: method 2 is no longer supported
use_MonteCarlo              = 'y';  % SSS; use MonteCarlo uncertainty estimates?


% step 2b: if method 2 is chosen, load MonteCarloSet ----------------------
if strcmp(use_MonteCarlo, 'y') && MonteCarlo_method == 2
    
    % define MonteCarloSet  = name of Monte Carlo results with file directory (loads 'MonteCarloStore')
    MonteCarloSet_FileDirectory = '';   % NOTE laptop MAC directory format
    MonteCarloSet               = '';   % sampling interval = +/1 1 STD sampling (1 std in denominator)
    MonteCarloSet               = [MonteCarloSet_FileDirectory MonteCarloSet];
 %   
    load(MonteCarloSet, 'MonteCarloStore')
else
    MonteCarloStore          = [];  % need this line for MonteCarlo_method = 1
end
% *************************************************************************





% *************************************************************************
% STEP 3: ECOTRAN conversion-----------------------------------------------
% step 3a: perform ECOTRAN conversion -------------------------------------
[ECOTRAN]                           = ECOTRANheart_09032021(EwEResult, MonteCarloStore);
% *************************************************************************





% *************************************************************************
% STEP 4: Generate E2E Monte Carlo models based on ECOTRAN EnergyBudget----
%         Start with the one original ECOTRAN base model and generate a set of 
%           Monte Carlo models from the ECOTRAN EnergyBudget & ConsumptionBudget
%           matrices using predefined CV values
%         NOTE: this is Monte Carlo method 1

% step 4a: generate a set of Monte Carlo models ---------------------------
if strcmp(use_MonteCarlo, 'y') && MonteCarlo_method == 1
    
    num_MC              = 1000;       % SSS set this value
    ECOTRAN.num_MC      = num_MC;
    
    PEDIGREE.ee_eggs_CV                               = 0.1; % SSS egg pedigree W.R.T. production budget for all groups; (CV); (scaler)
    PEDIGREE.BacterialMTBLSM_CV                       = 0.1; % SSS pedigree for implicit bacterial metabolism of terminal detritus (CV)
	PEDIGREE.Oxidation_NH4_CV                         = 0.1;	% fraction of NH4 produced oxidized directly back to NO3 abiologically; QQQ scaler?? (vertical vector: num_NH4 X 1)??
    PEDIGREE.NutrientUptake_CV                        = 0.1; % SSS pedigree for nutrient uptake by primary producers (CV)
    ECOTRAN_PEDIGREE                               	  = f_E2Epedigree_08042020(ECOTRAN, PEDIGREE); % NEW!!!
    
%             % SSS use for standardized pedigree
%             %     overwrite the pedigree values from the ECOPATH (EwE) model from VisualBasic file (.csv format)
%             [rows, clms]                                = size(ECOTRAN_PEDIGREE.EnergyBudget_CV);
%             ECOTRAN_PEDIGREE.EnergyBudget_CV            = 0.75 * ones(rows, clms);
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV       = zeros(7, clms);
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV(1, :) = 0.25; % feces
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV(2, :) = 0.25; % metabolism
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV(3, :) = 0.25; % eggs
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV(4, :) = 0.5;	% predation
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV(5, :) = 0.5;  % senescence
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV(6, :) = 0.5;  % ba
%             ECOTRAN_PEDIGREE.ConsumptionBudget_CV(7, :) = 0.5;  % em
% 	          ECOTRAN_PEDIGREE.DiscardFraction_CV         = 0.05 * ECOTRAN_PEDIGREE.DiscardFraction_CV; % QQQ reduce DiscardFraction_CV
    
    MonteCarloConditions.num_MC               	    = num_MC;	% number of random Monte Carlo models to generate
    MonteCarloConditions.DistributionType         	= 'normal';	% SSS distribution to draw random values from ('normal' or 'uniform'); (NOTE: code not fully proofed for uniform)
    
    ECOTRAN_MC                                      = f_E2E_MonteCarlo_08042020(MonteCarloConditions, ECOTRAN, ECOTRAN_PEDIGREE);
    EnergyBudget_MC                                 = ECOTRAN_MC.EnergyBudget_MC;
    ConsumptionBudget_MC                            = ECOTRAN_MC.ConsumptionBudget_MC;
    DiscardFraction_MC                           	= ECOTRAN_MC.DiscardFraction_MC;
	ECOTRAN_MC.num_MC                               = num_MC;

elseif strcmp(use_MonteCarlo, 'n')
    
	num_MC                          = 1;       % only the "type" model is used
    EnergyBudget_MC                 = ECOTRAN.EnergyBudget;
    ConsumptionBudget_MC            = ECOTRAN.ConsumptionBudget;
    DiscardFraction_MC              = ECOTRAN.DiscardFraction;
    
    ECOTRAN_MC.num_MC               = num_MC;
	ECOTRAN_MC.EnergyBudget_MC      = ECOTRAN.EnergyBudget; % (3D matrix: num_grps (consumers) X num_grps (producers) X num_MC (1))
    ECOTRAN_MC.ConsumptionBudget_MC	= ECOTRAN.ConsumptionBudget; % (3D matrix: 7 X num_grps (producers) X num_MC (1))
    ECOTRAN_MC.DiscardFraction_MC	= ECOTRAN.DiscardFraction;

end % (choose MonteCarlo_method)
% -------------------------------------------------------------------------


% step 4b: read in ECOTRAN structure variables ----------------------------
%          (so that no changes are made to original values)
GroupType                           = ECOTRAN.GroupType;
label                            	= ECOTRAN.label;
% EnergyBudget_MC                 	  = ECOTRAN.EnergyBudget;
biomass                          	= ECOTRAN.biomass;              % (vertical vector: num_grps X 1); note inclusion of separately constructed regional models
pb                               	= ECOTRAN.pb;                   % (vertical vector: num_grps X 1)
qb                               	= ECOTRAN.qb;                   % (vertical vector: num_grps X 1)
fate_feces                       	= ECOTRAN.fate_feces;
fate_metabolism                  	= ECOTRAN.fate_metabolism;
fate_eggs                        	= ECOTRAN.fate_eggs;
fate_senescence                  	= ECOTRAN.fate_senescence;
ProductionLossScaler             	= ECOTRAN.ProductionLossScaler;	% (vertical vector: num_grps X 1)
RetentionScaler                 	= ECOTRAN.RetentionScaler;      % sensitivity to advection & mixing (0 = more advection <--> less advection =1); (vertical vector: num_grps X 1)
FunctionalResponseParams         	= ECOTRAN.FunctionalResponseParams;
num_grps                            = ECOTRAN.num_grps;             % number of model groups
num_MC                              = ECOTRAN.num_MC;               % number of Monte Carlo models
% TransferEfficiency                  = ECOTRAN.TransferEfficiency;	  % gets redefined manually below
% -------------------------------------------------------------------------


% step 4c: find detritus, nutrients, ba & em ------------------------------
%          by row address in EnergyBudget_MC
looky_NO3                       	= find(GroupType        == ECOTRAN.GroupTypeDef_NO3);
looky_plgcNH4                       = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4);
looky_bnthNH4                       = find(GroupType        == ECOTRAN.GroupTypeDef_bnthNH4);
looky_NH4                           = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4 | GroupType == ECOTRAN.GroupTypeDef_bnthNH4);
looky_nutrients                     = find(floor(GroupType)	== ECOTRAN.GroupTypeDef_ANYNitroNutr);	% row addresses of nutrients
looky_ANYPrimaryProducer        	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYPrimaryProd);
looky_macroalgae                  	= find(GroupType        == ECOTRAN.GroupTypeDef_Macrophytes);
looky_phytoplankton                 = looky_ANYPrimaryProducer(~ismember(looky_ANYPrimaryProducer, looky_macroalgae));
looky_ANYconsumer                 	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYConsumer);
looky_micrograzers                  = find(GroupType        == ECOTRAN.GroupTypeDef_micrograzers);
looky_bacteria                      = find(floor(GroupType) == ECOTRAN.GroupTypeDef_bacteria);
looky_eggs                          = find(GroupType        == ECOTRAN.GroupTypeDef_eggs);
looky_ANYdetritus                   = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYDetritus);
looky_terminalPLGCdetritus          = find(GroupType        == ECOTRAN.GroupTypeDef_terminalPlgcDetr);
looky_terminalBNTHdetritus          = find(GroupType        == ECOTRAN.GroupTypeDef_terminalBnthDetr);
looky_eggsANDdetritus               = sort([looky_ANYdetritus; looky_eggs]);
looky_livingANDdetritus             = sort([looky_ANYPrimaryProducer; looky_ANYconsumer; looky_bacteria; looky_eggsANDdetritus]);
looky_terminalANYdetritus           = find(GroupType == ECOTRAN.GroupTypeDef_terminalPlgcDetr | GroupType == ECOTRAN.GroupTypeDef_terminalBnthDetr);
looky_fleets                        = find(floor(GroupType) == ECOTRAN.GroupTypeDef_fleet);
looky_livingANDfleets               = [looky_ANYPrimaryProducer; looky_ANYconsumer; looky_bacteria; looky_fleets]; % includes primary producers & bacteria
looky_NONnutrients                  = sort([looky_livingANDdetritus; looky_fleets]);	% addresses of all groups EXCEPT nutrients (needed to append nutrients)
looky_nonNO3                     	= 1:num_grps;
looky_nonNO3(looky_NO3)          	= [];

num_nutrients                   	= length(looky_nutrients);
num_NH4                          	= length(looky_NH4);
num_ANYPrimaryProd                	= length(looky_ANYPrimaryProducer);
num_phytoplankton                	= length(looky_phytoplankton);
num_macroalgae                   	= length(looky_macroalgae);
num_ANYconsumers                	= length(looky_ANYconsumer);
num_livingANDfleets                 = length(looky_livingANDfleets);
num_eggs                         	= length(looky_eggs);
num_ANYdetritus                   	= length(looky_ANYdetritus);
% -------------------------------------------------------------------------







% *************************************************************************
% STEP 5: define TransferEfficiency terms----------------------------------
%         SET MANUALLY: TransferEfficiency = 1 for all groups because we are
%                       working with the consumption budget matrix that tracks the fate of ALL consumption (not just predation)
%                       system losses (groups where TransferEfficiency < 1 include benthic detritus, fishery, (others?)
%          NOTE: terminal benthic detritus column in EnergyBudget sums to 1 and has NO OtherMortality, 
%                TE is the only means of removing terminal benthic detritus from system
TransferEfficiency      = ones(num_MC, num_grps); % re-initialize as rows of horizontal vector of ones
TransferEfficiency(:, [looky_terminalBNTHdetritus]) = 0.1; % NOTE: I use 0.1 TE for terminal benthic detritus as a standard default (JRuz 9/27/2018)
% *************************************************************************





% *************************************************************************
% STEP 6: define production loss fractions---------------------------------
%         (FFF in future, might use vector read in from .csv file; but for now, this is all handled here)
% LX                   = 0; % NOTE: set to 0 to ignore physical losses; set to 1 for the default relative upwelling index (1 is used for the average yearly relative upwelling index)
LX_base                       = ones(1, num_grps);
ProductionLossScaler_base     = zeros(1, num_grps); % initialize
PhysicalLossFraction_base     = ProductionLossScaler_base .* LX_base;
% *************************************************************************


SmlPhyto_rowclm                 = 4;
LargePhyto_rowclm               = 5;

% *************************************************************************
% STEP 7: initialize InputProductionVector---------------------------------
%         NOTE: here you can run a scenario with an alternate driver
InputProductionVector_base                          = zeros(1, num_grps);   % initialize InputProductionVector_base
InputProductionVector_scenario                      = zeros(1, num_grps);   % initialize InputProductionVector_scenario

phytoplankton_base                                  = pb(SmlPhyto_rowclm)   .* biomass(SmlPhyto_rowclm);	% (t/km2/y); (vertical vector: num_grps X 1)

InputProductionVector_base(SmlPhyto_rowclm)    	= phytoplankton_base;	% (t/km2/y); (vertical vector: num_grps X 1)

InputProductionVector_scenario(SmlPhyto_rowclm)	= phytoplankton_base;	% base for now, can change these for a scenario; (t/km2/y); (vertical vector: num_grps X 1)



% step 7b: evaluate type of driver (nutrient vs living group) -------------
% NOTE: THIS IS NEEDED IF DRIVEN BY ANYTHING OTHER THAN NO3
TestInput_base                      = InputProductionVector_base;
TestInput_scenario                  = InputProductionVector_scenario;
TestInput_base(looky_nutrients)     = [];
TestInput_scenario(looky_nutrients) = [];
if max(TestInput_base > 0 | TestInput_scenario > 0) 
    disp('-->WARNING: non-NO3 driver. Deactivating nutrient recycling; set Nutrient columns in EnergyBudget to 0');
    EnergyBudget_MC(:, looky_nutrients, :)	= 0; % shut off recycling; % deactivate nutrient recycling; set nutrient columns in EnergyBudget to 0 so that recycled nutrients don't flow into any other group
    fate_metabolism                         = ECOTRAN_MC.fate_metabolism; % (3D matrix: num_nutrients X num_grps X num_MC)
	fate_metabolism(:, looky_nutrients, :)	= 0; % deactivate nutrient recycling; set nutrient columns in EnergyBudget to 0 so that recycled nutrients don't flow into any other group
    fate_predation                          = ECOTRAN_MC.fate_predation; % (3D matrix: num_livingANDfleets X num_grps X num_MC)
    fate_predation(:, looky_nutrients, :)	= 0; % deactivate nutrient recycling; set nutrient columns in EnergyBudget to 0 so that recycled nutrients don't flow into any other group
    % NOTE: this automatic switch assumes that both the base and the scenario are using non-nutrient drivers
    % Will need to add a more sophisticated switch to allow a nutrient driver base and a non-nutrient driver scenario (or vice-versa)
    % FFF in future, can have nutrient and non-nutrient drives for different subregions
else
    disp('-->WARNING: NO3 driver. Be sure to keep nutrient recycling ACTIVATED'); 
end
% *************************************************************************





% *************************************************************************
% STEP 8: calculate base productivity--------------------------------------
Production_base	= zeros(num_grps, num_MC); % initialize; (2D matrix: num_grps X num_MC)
for MonteCarlo_loop = 1:num_MC    
    current_EnergyBudget                = EnergyBudget_MC(:, :, MonteCarlo_loop);
    current_TransferEfficiency          = TransferEfficiency(MonteCarlo_loop, :);
    Production_base(:, MonteCarlo_loop)	= f_WebProductivity_03272019(current_TransferEfficiency, current_EnergyBudget, InputProductionVector_base, PhysicalLossFraction_base); % production (actually CONSUMPTION) rates as "initial" or "mean" conditions; (mmole N/m3/d); (2D matrix: num_grps X num_MC)
end % (MonteCarlo_loop)
% *************************************************************************

% if you want all MC production values, change from Production_base(:,1) to Production_base(:,:)
csvwrite(strcat(ReadFile_directory,regexprep(BiologicalModel_name,".csv",""),"_Production.csv"),Production_base(:,1))
