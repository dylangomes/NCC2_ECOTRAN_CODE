% ECOTRAN_FootprintReach_NCC_01102022
%
% Calculate FOOTPRINT & REACH metrics upon individual groups and at the ecosystem level
%
% The FOOTPRINT metrics are the fraction of each PRODUCER group's production
%       flowing to CONSUMER = TraceGroup. The FOOTPRINT is the fraction of 
%       each PRODUCER's production required to support a particular (TraceGroup)
%       CONSUMER group. The code calculates FOOTPRINT for all functional groups as TraceGroup.
%
% The REACH metrics are the fraction of each CONSUMER group's production ultimately 
%       originating from PRODUCER = TraceGroup. The REACH is the fraction of 
%       a particular (TraceGroup) PRODUCER group's production going to support 
%       each CONSUMER group. The code calculates REACH for all functional groups as TraceGroup.
%
% Calculate the FOOTPRINT and REACH traces for each trophic linkage for 
%       one (1) specific functional group of interest = TraceGroup.
% FOOTPRINT_trace is the fraction of each trophic link ultimately contributing
%       to production of CONSUMER = TraceGroup
%       It is the relative contribution of each linkage to the production of 
%       CONSUMER = TraceGroup.
% REACH_trace is the fraction of each trophic link ultimately originating
%       from PRODUCER = TraceGroup.
%       It is the fraction of energy within each linkge ultimately 
%       originating from PRODUCER = TraceGroup.
%
% The user defines the TraceGroup
%
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
%       f_FootprintReach_05182022_b                   calculate FOOTPRINT & REACH metrics and network trace variables
%           f_WebProductivity_03272019              calculate consumption rate (q) for each group given a group driver (e.g. phytoplankton production) input production rate
%           f_DietTrace_08032020                    use DIET matrix to calculate the "REACH" metrics of each functional group, and calculate the FOOTPRINT & REACH traces for each trophic linkage for one (1) specific functional group of interest = TraceGroup.
%
%       p_WebPlotter_08042020
%
% returns:
%       DIET_MC         	- diet of each consumer
%                       	- (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%
%       FOOTPRINT_array_MC	- "FOOTPRINT" of each TraceGroup = consumer (rows) upon each producer (columns)
%                        	- the fraction of each producer group's production flowing to each consumer = TraceGroup
%                         	- (3D matrix: num_grps (consumer) X num_grps (producer) X num_MC)
%                         	- NOTE: use for web plotting footprint box colors relative to TraceGroup
%                         	- NOTE: does NOT include import diet.
%                                   There is no footprint calculated for producer = import
%                                   because the value of import production can be arbitrary.
%
%       REACH_array_MC      - "REACH" of TraceGroup = producer (rows) to each consumer (columns)
%                           - the fraction of each consumer group's production ultimately originating from producer = TraceGroup                 
%                           - (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%                           - NOTE: use to plot REACH box colors in food web diagram
%                           - NOTE: formerly called "TraceFraction_upward"
%                           - NOTE: DOES include import diet
%
%       FOOTPRINT_trace_MC	- Fraction of each trophic link ultimately contributing to production of CONSUMER = TraceGroup
%                           - the relative contribution of each linkage to production of CONSUMER = TraceGroup
%                           - (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%                           - NOTE: use to plot FOOTPRINT arrow colors in food web diagram
%                           - NOTE: formerly called "DietTrace_downward"
%                           - NOTE: DOES include import diet
%
%       REACH_trace_MC      - Fraction of each trophic link ultimately originating from PRODUCER = TraceGroup
%                           - It is the fraction of energy within each linkge ultimately originating from PRODUCER = TraceGroup.
%                           - (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%                           - NOTE: use to plot REACH arrow colors in food web diagram
%                           - NOTE: formerly called "DietTrace_upward"
%                           - NOTE: DOES include import diet
%
%       FOOTPRINT_trace_MC	- SYSTEM-LEVEL footprint = ratio of total TraceGroup footprint on all PRODUCERS over total production of all CONSUMER groups
%                           - total consumer production of ecosystem excludes microzooplankton
%                           - (vertical vector: num_MC X 1)
%
%       REACH_system_MC     - SYSTEM-LEVEL reach = ratio of total TraceGroup production going to all consumers over total production of all CONSUMER groups
%                           - total consumer production of ecosystem excludes microzooplankton
%                           - (vertical vector: num_MC X 1)
%
%       And, the means & standard deviations of each variable across all Monte Carlo models 
%           (NOTE: row or layer 1 = mean; row or layer 2 = standard deviation)
%               DIET_recalculated	(3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
%               FOOTPRINT_array     (3D matrix: num_grps (consumers) X num_grps (producers) X 2)
%               REACH_array         (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
%               FOOTPRINT_trace     (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
%               REACH_trace         (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
%               FOOTPRINT_system	(vertical vector: 2 X 1)
%               REACH_system        (vertical vector: 2 X 1)
%
% revision date: 1-9-2023


% *************************************************************************
% STEP 1: load & aggregate EwE results-------------------------------------
fname_FootprintReach     = 'ECOTRAN_FootprintReach_NCC_01102022'; % save name of this m-file to keep in saved model results

% step 1a: define food web model to use -----------------------------------
ReadFile_directory      = 'C:/Users/dgome/Documents/NCC_Comparison/';
BiologicalModel_name	= {'pre_MHW_Web/NCC_11242020_JimsOld_updatedKrill_10052022.csv',...
                             'post_MHW_Web/NCC2_NCC1Matched_10042022.csv'};
z = []; % initiate object to save to
num_MC  = 1000;  % set value for number of MC models to produce

% select functional group
FG = {'phytoplankton','copepods','jellyfish','pyrosomes',...
    'E. pacifica','T. spinifera',...
    'sardine','anchovy','herring',...
    'hake','rockfish','jack mackerel','dogfish',...
    'dungeness',...
    'pinnipeds','murre'...
    };


% *************************************
% stop editing for general use (make all selections above)
% *************************************


for m = 1:size(BiologicalModel_name,2)
% step 1b: load ECOPATH (EwE) model from Aydin VisualBasic file (.csv format)
readFile            	= strcat(ReadFile_directory, BiologicalModel_name(m));
dat                  	= f_readEwEcsv_10pp_07072021(char(readFile));	% use for models with up to 10 primary producers (use for 10 pp version of Aydin .xlsm, regardless of number of actual pp in model)
display(['Processing ecosystem model ' num2str(m) ' of ' num2str(size(BiologicalModel_name,2)) ' : ' char(BiologicalModel_name(m))])
display(['Number of MC models: ' num2str(num_MC)])

% step 1c: aggregate model results & prep EwEResult for analysis ----------
[EwEResult, PEDIGREE] 	= f_AggregateBiologicalModel_02052021(dat);
% *************************************************************************

% *************************************************************************
% STEP 2: Monte Carlo EwE models------------------------
%      start with the one original ECOTRAN base model and generate a set of Monte Carlo E2E models of the ECOTRAN production matrix using a predefined CV (or defualt pedigree CV = 0.5)
use_MonteCarlo              = 'y';  % use MonteCarlo uncertainty estimates?
MonteCarloStore          = [];
% *************************************************************************

% STEP 3: ECOTRAN conversion-----------------------------------------------
[ECOTRAN]                           = ECOTRANheart_09032021(EwEResult, MonteCarloStore);

% *************************************************************************
% STEP 4: Generate E2E Monte Carlo models based on ECOTRAN EnergyBudget----
%         Start with the one original ECOTRAN base model and generate a set of Monte Carlo models from the ECOTRAN EnergyBudget & ConsumptionBudget matrices using predefined CV values
% step 4a: generate a set of Monte Carlo models ---------------------------
if strcmp(use_MonteCarlo, 'y')
    
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
% *************************************************************************



% *************************************************************************
% STEP 5: define TransferEfficiency terms----------------------------------
%         SET MANUALLY: TransferEfficiency = 1 for all groups because we are
%                       working with the consumption budget matrix that tracks the fate of ALL consumption (not just predation)
%                       system losses (groups where TransferEfficiency < 1 include benthic detritus, fishery, (others?)
%          NOTE: terminal benthic detritus column in EnergyBudget sums to 1 and has NO OtherMortality, 
%                TE is the only means of removing terminal benthic detritus from system
TransferEfficiency                  = ones(num_MC, num_grps); % re-initialize as rows of horizontal vector of ones
TransferEfficiency(:, [looky_terminalBNTHdetritus]) = 0.1; % NOTE: I use 0.1 TE for terminal benthic detritus as a standard default (JRuz 9/27/2018)
ECOTRAN_MC.TransferEfficiency_MC	= TransferEfficiency;
% *************************************************************************


% *************************************************************************
% STEP 6: define production loss fractions---------------------------------
%         (FFF in future, might use vector read in from .csv file; but for now, this is all handled here)
% LX                   = 0; % NOTE: set to 0 to ignore physical losses; set to 1 for the default relative upwelling index (1 is used for the average yearly relative upwelling index)
LX                          = ones(1, num_grps);
ProductionLossScaler        = zeros(1, num_grps); % initialize
PhysicalLossFraction        = ProductionLossScaler .* LX;
% *************************************************************************



% STEP 7: calculate FOOTPRINT & REACH metrics------------------------------
% The FOOTPRINT metrics are the fraction of each PRODUCER group's production
%       flowing to CONSUMER = TraceGroup. Code calculates FOOTPRINT for
%       all functional groups as TraceGroup.
%
% The REACH metrics are the fraction of each CONSUMER group's production ultimately 
%       originating from PRODUCER = TraceGroup. Code calculates REACH for
%       all functional groups as TraceGroup.
%
% Calculate the FOOTPRINT and REACH traces for each trophic linkage for 
%       one (1) specific functional group of interest = TraceGroup.
% FOOTPRINT_trace is the fraction of each trophic link ultimately contributing
%       to production of CONSUMER = TraceGroup
%       It is the relative contribution of each linkage to the production of 
%       CONSUMER = TraceGroup.
% REACH_trace is the fraction of each trophic link ultimately originating
%       from PRODUCER = TraceGroup.
%       It is the fraction of energy within each linkge ultimately 
%       originating from PRODUCER = TraceGroup.
%

% step 7a: define TraceGroup & calculate metrics --------------------------


for j = 1:size(FG,2)
   display(['Functional Group: ' num2str(j) ' of ' num2str(size(FG,2)) ' : ' char(FG(j))])
TraceGroup              = find(contains(ECOTRAN.label,FG(j)));

for t = 1:size(TraceGroup,1)
[DIET_MC, FOOTPRINT_array_MC, REACH_array_MC, FOOTPRINT_trace_MC, REACH_trace_MC, FOOTPRINT_system_MC, REACH_system_MC] = ...
    f_FootprintReach_05182022_b(ECOTRAN, ECOTRAN_MC, PhysicalLossFraction, TraceGroup(t));
% -------------------------------------------------------------------------


% STEP 7b: pack-up results-------------------------------------------------
%          take averages % standard deviations of Monte Carlo models
DIET_recalculated(:, :, 1)	= mean(DIET_MC, 3);              % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
DIET_recalculated(:, :, 2)	= std(DIET_MC, 0, 3);            % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)

FOOTPRINT_array(:, :, 1)	= mean(FOOTPRINT_array_MC, 3);   % (3D matrix: num_grps (consumers) X num_grps (producers) X 2)
FOOTPRINT_array(:, :, 2)	= std(FOOTPRINT_array_MC, 0, 3); % (3D matrix: num_grps (consumers) X num_grps (producers) X 2)

REACH_array(:, :, 1)        = mean(REACH_array_MC, 3);       % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
REACH_array(:, :, 2)        = std(REACH_array_MC, 0, 3);     % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)

FOOTPRINT_trace(:, :, 1)	= mean(FOOTPRINT_trace_MC, 3);	 % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
FOOTPRINT_trace(:, :, 2)	= std(FOOTPRINT_trace_MC, 0, 3); % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)

REACH_trace(:, :, 1)        = mean(REACH_trace_MC, 3);	     % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
REACH_trace(:, :, 2)        = std(REACH_trace_MC, 0, 3);	 % (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)

FOOTPRINT_system(1)         = mean(FOOTPRINT_system_MC);     % (vertical vector: 2 X 1)
FOOTPRINT_system(2)         = std(FOOTPRINT_system_MC);      % (vertical vector: 2 X 1)

REACH_system(1)             = mean(REACH_system_MC);         % (vertical vector: 2 X 1)
REACH_system(2)             = std(REACH_system_MC);          % (vertical vector: 2 X 1)
% *************************************************************************





% % *************************************************************************
% % STEP 8: plot food webs---------------------------------------------------
% % step 8a: box x-position -------------------------------------------------
% %           NOTE: The order of x-pos is the same as the ECOPATH order
% x_position = [random('Uniform',0,10,num_grps+1,1)];
% % -------------------------------------------------------------------------
% 
% 
% % step 8b: y-position (= TL) ----------------------------------------------
% %           NOTE: uses EwE trophic levels, but you must add fisheries & import prey
% TL_position             = ECOTRAN.TL;
% TL_position(num_grps+1)	= 1.5; % Import diet
% % -------------------------------------------------------------------------
% 
% 
% % step 8c: box labels -----------------------------------------------------
% %           NOTE: you must add fisheries & import prey
% % GRPlabels               = ([1:num_grps+1]');
% GRPlabels               = regexprep(ECOTRAN.label,'.*: ()','$1');
% GRPlabels{(num_grps+1)}	= {'import'; 'diet'};
% 
% Edit              = find(contains(ECOTRAN.label,'amphipods isopods'));
% x_position(Edit) = 2;
% % overwrite the ECOTRAN.label designations if wanted
% % GRPlabels{NO3_rowclm}                       = {};
% % GRPlabels{pelagicNH4_rowclm}             	= {};
% % GRPlabels{benthicNH4_rowclm}            	= {};
% % GRPlabels{SmlPhyto_rowclm}                  = {'small'; 'phytoplankton'};
% % -------------------------------------------------------------------------


% % step 8d: plot web -------------------------------------------------------
% grid_switch = 'off';
% 
% % p_WebPlotter_08042020(ECOTRAN, FOOTPRINT_array, FOOTPRINT_trace, REACH_array, REACH_trace, x_position, TL_position, GRPlabels, TraceGroup, grid_switch)
% % p_WebPlotter_08042020(ECOTRAN, FOOTPRINT_array, REACH_array, x_position, TL_position, GRPlabels, TraceGroup, grid_switch)
% 
% % save
% % NAMES = repmat({FG(j)},num_MC,1) 


for i = 1:num_MC
z = [z; FG(j) TraceGroup(t) FOOTPRINT_system_MC(i) REACH_system_MC(i) char(BiologicalModel_name(m))];
end
    end
        end
            end

writecell(z, strcat(ReadFile_directory,"Footprint_Reach_",date,".csv"))
% print -dtiff -r300 'testtoss_TEST8b.tif'
% *************************************************************************

% end m-file***************************************************************