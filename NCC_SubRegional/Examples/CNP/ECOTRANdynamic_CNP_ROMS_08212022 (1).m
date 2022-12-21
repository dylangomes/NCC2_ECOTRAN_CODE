% ECOTRANdynamic_CNP_ROMS_08212022
% run a dynamic model over time
% by Jim Ruzicka
%
% calls:
%       f_readEwEcsv_10pp_08042020                	read in ECOPATH (EwE) model from VisualBasic .csv file and store as variable 'dat'; (use for VisualBasic food web files with 10 primary producers)
%       f_AggregateBiologicalModel_02052021        	prepare EwE model for use by ECOTRAN; also, aggregate functional groups here if wanted
%           f_calcEE_12292020                          calculate Ecotrophic Efficiency
%           f_VarianceDivision_12132018                calculate the variance of one term divided by another term
%           f_VarianceMultiplication_12132018          calculate the variance of two products
%
%       ECOTRANheart_09032021                       generate an ECOTRAN (E2E) model; the heart of ECOTRAN
%           f_ECOfunction_09032021                     returns a single ECOTRAN model for 1 "type" EwE model or 1 MonteCarlo EwE model
%               f_RedistributeCannibalism_11202019       remove cannibalism terms on diagonal of matrix EwE_diet
%               f_calcEE_12292020                        calculate Ecotrophic Efficiency 
%               f_calcPredationBudget_12102019      	 for each producer, p, and consumer, c: ((b_pc * Q_c) / M2_p)) = the fraction of total predation on each producer, p, going to each consumer, c; (2D matrix: num_grps X num_grps; consumers X producers)
%
%       f_E2Epedigree_08042020                      calculate uncertainty terms for all elements of the ECOTRAN EnergyBudget (A_cp) as Coefficients of Variation (CV)
%           f_VarianceDivision_12132018                calculate the variance of one term divided by another term
%           f_VarianceMultiplication_12132018          calculate the variance of two products
%       f_E2E_MonteCarlo_08042020                   generate a set of randomly generated ECOTRAN models, multiple alternate versions of the EnergyBudget (A_cp); model 1 of the stack is the "type" model defined by the parameters of the VisualBasic .csv file
%
%       f_OrdinalDate                               calculate ordinal dates; day of year with January 1 of any year = 1
%
%       alternate physical models:
%           f_ECOTRANphysics_upwelling_05242022      GoMexOcn system: define the physical geometry of the model and prepare physical advection, mixing, and sinking rate time-series
%               f_read_ERDdata_05122021	calculate daily median upwelling intensities from ERD products (http://http://www.pfeg.noaa.gov/products/PFELData/upwell/6_hourly/upwell45N125W [may be a defunct link])
%               f_LightIntensity_12112020        	solar light intensity for given location and time
%               f_CompactFluxTimeSeries_11182019	compact flux time series when arranged as 3D matrix (time X source box X destiny box)
%               f_EvaluateFluxBalance_11192019      examine for flux time-series for imbalances IN & OUT of invidual boxes and IN & OUT of the overall domain
%               f_UnCompactFluxTimeSeries_12112019	UnCompact a flux time series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole
%                   f_CalcNetFlux_12112019          calculate net flux into and net flux out of each model box and across outer domain boundaries
%
%       f_ECOTRANmigration_NCC_05122021        area of overlap of neighboring model sub-domains
%
%       f_CompactFluxTimeSeries_11182019            compact physical flux 3D matrix
%
%       f_DVMsinusoid_02192021                      Calculate DVM flux rates between all model domain boxes for each functional group at each time point; Uses a daily sinusoidal migration pattern
%
%       f_StaticProductionTimeseries_09042017       (deactivated, not yet converted to ECOTRAN II) calculate production rates for all model groups under a given driver (NO3) input rate; use for ocean boundary conditions
%           f_WebProductivityWLoss                    calculate production rates of all groups under a given driver (e.g., NO3 or primary production); also accounts for defined rates of group production export when running static scenarios
%
%       f_FunctionalResponse_MonteCarlo_05132021    prepare array of vulnerability terms and allows for random generation of functional response terms within a predefined uncertainty level
%
%       f_SeasonalPB_PrimaryProducers               (deactivated)
%       f_MichaelisMenten_05152016                  (optional) phytoplankton uptake rate of NO3 & NH4 (mmole N/m3/d) & p/b @ t as calculated from Michaelis-Menten uptake kinetics
%
%       f_InitialProductionRates_11202019          	calculate initial or mean production conditions
%           f_WebProductivity_03272019                calculate production rates of all groups under a given driver (e.g., NO3 or primary production); also accounts for defined rates of group production export when running static scenarios
%
%       2 ODE solver options:
%       f_PrepMexODE_10272020                       prepare ECOTRAN variables for using C++ ODE solver mex function. Pack parameters & drivers along proper dimensions
%           f_unspoolMATRIX_04282020                  linearize multidimenional matrices up to 4-D for use in C++
%       mex_ECOTRANode_03092021_D                   solve the ecosystem ODE in C++
%
%       f_ECOTRANode_11202019                       Ordinary Differential Equation (use this for getting soln at each time-point; default is for reflective boundary, but has built-in options for defined boundary conditions)
%           f_PhysicalFlux_intraODE_09092019          calculate biomass fluxes into and out of each box and across domain boundaries (mmoles N/m3/d)	(3D matrix: 1 X num_grps X num_boxes DESTINATION)
%           f_MichaelisMenten_05152016                (optional) Michaelis-Menton uptake kinetics for phytoplankton
%
% returns:
%       store_T                     time-series of day umbers
%       store_ProductionRates       time-series of production rates for each functional group; (t/km2/d) 
%
% NOTE: code cannot currently accomodate seasonal changes in flows to non-terminal detritus pools (e.g., fisheries and detritus columns of EnergyBudget)
%
% revision date: 7-20-2022
%       5/23/2022 adding ThorntonLessem max consumption rate scaling
%       5/24/2022 updated f_ECOTRANphysics_upwelling_05242022 to correct NH-Line temperature time-series error
%       5/27/2022 added loop to combine ROMS output years
%       6/4/2022 new C++ solver with ThrontonLessem temperature
%       6/9/2022 QQQ need to add temperature timeseries to make Thornton-Lessem work in C++ solver
%       7/18/2022 JR corrected Q10 expense accounting in ConsumptionBudget
%       8/21/2022 JR corrected use of pb & qb in solver, added MiningPlume effects model


%% *************************************************************************
% STEP 1: load & aggregate EwE results-------------------------------------
% step 1a: set operating conditions ---------------------------------------
fname_ECOTRANdynamic	= 'ECOTRANdynamic_CNP_ROMS_08212022'; % save name of this m-file to keep in saved model results

switch_ODEsolver            = 'CppSolver';                          % SOLVER 1: C++
% switch_ODEsolver            = 'MatlabSolver';                       % SOLVER 1: Matlab

% switch_SubModel             = 'identical_SubModel';                 % OPTION 1: use the same food web across the shelf
switch_SubModel             = 'independent_SubModel';               % OPTION 2: use independently defined webs for each shelf zone

% switch_INITIALproduction	= 'INITIALproduction_MichaelisMenton';	% METHOD 1: use for driving initial model conditions with primary production defined by Michaelis-Menton uptake
% switch_INITIALproduction	= 'INITIALproduction_pb';               % METHOD 2: use for driving initial model conditions with primary production defined by p = [(p/b) * b]
switch_INITIALproduction	= 'INITIALproduction_SubModel';         % METHOD 3: use for driving initial model conditions with values loaded along with regional sub-model definitions 
% switch_INITIALproduction	= 'INITIALproduction_nutrients';      	% METHOD 4: use for driving initial model conditions with mean annual nutrient input rates

% switch_FunctionalResponse	= 'Linear'; % linear functional response
switch_FunctionalResponse	= 'NonLinear_default'; % NonLinear_default functional response
% switch_FunctionalResponse	= 'NonLinear_alt'; % NonLinear_default functional response

% switch_MonteCarlo           = 'MonteCarlo_build';	% generate (and optionally save) a stack of MonteCarlo food webs
% switch_MonteCarlo           = 'MonteCarlo_load';	% load a saved stack of MonteCarlo food webs
switch_MonteCarlo           = 'MonteCarlo_TypeModel';	% use NO MonteCarlo food webs

switch_PhysicalModel        = '2D_physics';
% switch_PhysicalModel        = '3D_ROMS';

switch_physiology           = 'physiology_ON';
% switch_physiology           = 'physiology_OFF';

% switch_ThorntonLessem_physiology = 'ThorntonLessem_ON';
switch_ThorntonLessem_physiology = 'ThorntonLessem_OFF';

% switch_MiningPlume_physiology = 'MiningPlume_ON';
switch_MiningPlume_physiology = 'MiningPlume_OFF';

switch_physiology_scenario  = 'up25perc'; % 'dwn10perc', 'base', 'up10perc', 'up25perc', 'dwn25perc', 'up40perc' 'range'
% -------------------------------------------------------------------------


% step 1b: define food web model to use -----------------------------------
% ReadFile_directory      = '\Users\gails\Desktop\CNP-ECOTRAN_11242021\3_ECOTRAN_code_files\'; % Jake's directory
ReadFile_directory      = '/Users/jimsebi/Documents/11_FoodWeb_models/4_ChoyCNP/Choy_Rebuild/'; % Jim's directory
BiologicalModel_name	= 'CNP_11182021_P.csv';

readFile            	= [ReadFile_directory BiologicalModel_name];

% step 1c: load ECOPATH (EwE) model from Aydin VisualBasic file (.csv format)
dat                  	= f_readEwEcsv_10pp_07072021(readFile);	% use for models with up to 10 primary producers

% step 1d: aggregate model results & prep EwEResult for analysis ----------
[EwEResult, PEDIGREE] 	= f_AggregateBiologicalModel_02052021(dat);

% step 1e: define filename and directory for saving results ---------------
% SaveFile_directory      = '\Users\gails\Desktop\'; % QQQ activate for Jake's directory
SaveFile_directory      = '/Users/jimsebi/Desktop/CNP_Hollings_Runs/'; % QQQ activate for Jim's directory
SaveFile_label          = 'CNP_testTOSS';
% *************************************************************************





%% ************************************************************************
% STEP 2: ECOTRAN conversion-----------------------------------------------
MonteCarloStore         = [];
[ECOTRAN]            	= ECOTRANheart_09032021(EwEResult, MonteCarloStore);
% *************************************************************************





%% *************************************************************************
% STEP 3: Generate E2E Monte Carlo models based on ECOTRAN EnergyBudget----
%         Start with the one original ECOTRAN base model and generate a set of 
%           Monte Carlo models from the ECOTRAN EnergyBudget & ConsumptionBudget
%           matrices using predefined CV values

switch switch_MonteCarlo

	case 'MonteCarlo_build' % generate (and save) a stack of Monte Carlo models
    
        num_MC              = 3;        % SSS set this value
        
        disp(['MonteCarlo: building stack of ' num2str(num_MC) ' food webs'])
        
        ECOTRAN.num_MC      = num_MC;

        PEDIGREE.ee_eggs_CV                               = 0.01; % SSS egg pedigree W.R.T. production budget for all groups; (CV); (scaler)
        PEDIGREE.BacterialMTBLSM_CV                       = 0.01; % SSS pedigree for implicit bacterial metabolism of terminal detritus (CV)
        PEDIGREE.Oxidation_NH4_CV                         = 0.01;	% fraction of NH4 produced oxidized directly back to NO3 abiologically; QQQ scaler?? (vertical vector: num_NH4 X 1)??
        PEDIGREE.NutrientUptake_CV                        = 0.01; % SSS pedigree for nutrient uptake by primary producers (CV)
        ECOTRAN_PEDIGREE                               	  = f_E2Epedigree_08042020(ECOTRAN, PEDIGREE); % NEW!!!

                % SSS use for standardized pedigree
                %     overwrite the pedigree values from the ECOPATH (EwE) model from VisualBasic file (.csv format)
                [rows, clms]                                = size(ECOTRAN_PEDIGREE.EnergyBudget_CV);
                ECOTRAN_PEDIGREE.EnergyBudget_CV            = 0.001 * ones(rows, clms); % changes how important predators are relative to eachother
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV       = zeros(7, clms);
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(1, :) = 0.05; % feces
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(2, :) = 0.05; % metabolism
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(3, :) = 0.05; % eggs
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(4, :) = 0.05;	% predation
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(5, :) = 0.05;  % senescence
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(6, :) = 0.05;  % ba
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(7, :) = 0.05;  % em
                ECOTRAN_PEDIGREE.DiscardFraction_CV         = 0.05 * ECOTRAN_PEDIGREE.DiscardFraction_CV; % QQQ reduce DiscardFraction_CV

        MonteCarloConditions.num_MC               	    = num_MC;	% number of random Monte Carlo models to generate
        MonteCarloConditions.DistributionType         	= 'normal';	% SSS distribution to draw random values from ('normal' or 'uniform'); (NOTE: code not fully proofed for uniform)

        ECOTRAN_MC                                      = f_E2E_MonteCarlo_08042020(MonteCarloConditions, ECOTRAN, ECOTRAN_PEDIGREE);
        EnergyBudget_MC                                 = ECOTRAN_MC.EnergyBudget_MC;
        ConsumptionBudget_MC                            = ECOTRAN_MC.ConsumptionBudget_MC;
        DiscardFraction_MC                              = ECOTRAN_MC.DiscardFraction_MC;
        ECOTRAN_MC.num_MC                               = num_MC;

%         % activate to save this stack of Monte Carlo food webs
%         filename_MC      = [SaveFile_directory 'MonteCarlo_NCC_stack_' date '.mat'];
%         disp(['MonteCarlo: SAVING stack: ' filename_MC])
%         save(filename_MC, 'ECOTRAN_MC')

    % end (case build_MonteCarlo) -----------------------------------------

    
    case 'MonteCarlo_load' % load a set of Monte Carlo models
        filename_MC = [SaveFile_directory 'MonteCarlo_NCC_stack_27-Jun-2022.mat']; % SSS be sure to give correct saved file name here
        disp(['MonteCarlo: LOADING stack: ' filename_MC])
        load(filename_MC, 'ECOTRAN_MC')
        
        num_MC                            	= ECOTRAN_MC.num_MC;
        EnergyBudget_MC                     = ECOTRAN_MC.EnergyBudget_MC;
        ConsumptionBudget_MC                = ECOTRAN_MC.ConsumptionBudget_MC;
        DiscardFraction_MC                  = ECOTRAN_MC.DiscardFraction_MC;
    % end (case load_MonteCarlo) ------------------------------------------
        
    
    case 'MonteCarlo_TypeModel'
        
        disp('MonteCarlo: using the defining TypeModel')
    
        num_MC                          = 1;       % only the "type" model is used
        EnergyBudget_MC                 = ECOTRAN.EnergyBudget;
        ConsumptionBudget_MC            = ECOTRAN.ConsumptionBudget;
        DiscardFraction_MC              = ECOTRAN.DiscardFraction;

        ECOTRAN_MC.num_MC               = num_MC;
        ECOTRAN_MC.EnergyBudget_MC      = ECOTRAN.EnergyBudget; % (3D matrix: num_grps (consumers) X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.ConsumptionBudget_MC	= ECOTRAN.ConsumptionBudget; % (3D matrix: 7 X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.DiscardFraction_MC	= ECOTRAN.DiscardFraction;
        
        ECOTRAN_MC.fate_metabolism      = ECOTRAN.fate_metabolism;	% (3D matrix: num_nutrients X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_eggs            = ECOTRAN.fate_eggs;        % (3D matrix: num_eggs X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_feces           = ECOTRAN.fate_feces;       % (3D matrix: num_ANYdetritus X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_senescence      = ECOTRAN.fate_senescence;  % (3D matrix: num_ANYdetritus X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_predation       = ECOTRAN.fate_predation;   % (3D matrix: num_livingANDfleets X num_grps (producers) X num_MC (1))
	% end (case MonteCarlo_TypeModel) -------------------------------------

end % (switch_MonteCarlo)
% *************************************************************************





%% ************************************************************************
% STEP 4: prep and pack ECOTRAN model parameters---------------------------

% step 4a: read in ECOTRAN structure variables ----------------------------
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


%% step 4b: find detritus, nutrients, ba & em -----------------------------
%           row addresses in EnergyBudget_MC
%           NOTE: ba = biomass accumulation term, em = emigration term

looky_NO3                       	= find(GroupType        == ECOTRAN.GroupTypeDef_NO3);
looky_plgcNH4                       = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4);
looky_bnthNH4                       = find(GroupType        == ECOTRAN.GroupTypeDef_bnthNH4);
looky_NH4                           = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4 | GroupType == ECOTRAN.GroupTypeDef_bnthNH4);
looky_nutrients                     = find(floor(GroupType)	== ECOTRAN.GroupTypeDef_ANYNitroNutr);	% row addresses of nutrients
looky_ANYPrimaryProducer        	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYPrimaryProd);
looky_macroalgae                  	= find(GroupType        == ECOTRAN.GroupTypeDef_Macrophytes);
looky_phytoplankton                 = looky_ANYPrimaryProducer(~ismember(looky_ANYPrimaryProducer, looky_macroalgae));
looky_ANYconsumer                 	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYConsumer);
% looky_PLGCbacteria                  = find(GroupType        == ECOTRAN.GroupTypeDef_plgcBacteria);
% looky_BNTHbacteria                  = find(GroupType        == ECOTRAN.GroupTypeDef_bnthBacteria);
% looky_ANYbacteria                   = find(GroupType        == ECOTRAN.GroupTypeDef_plgcBacteria | GroupType == ECOTRAN.GroupTypeDef_bnthBacteria);
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
num_NO3                             = length(looky_NO3);
num_NH4                          	= length(looky_NH4);
num_plgcNH4                         = length(looky_plgcNH4); % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
num_bnthNH4                         = length(looky_bnthNH4); % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
num_ANYPrimaryProd                	= length(looky_ANYPrimaryProducer);
num_phytoplankton                	= length(looky_phytoplankton);
num_macroalgae                   	= length(looky_macroalgae);
num_ANYconsumers                	= length(looky_ANYconsumer);
num_fleets                          = length(looky_fleets);
num_predators                       = num_ANYconsumers + num_fleets;
num_livingANDfleets                 = length(looky_livingANDfleets);
num_eggs                         	= length(looky_eggs);
num_ANYdetritus                   	= length(looky_ANYdetritus);
% -------------------------------------------------------------------------


% step 4c: pack variables for ODE solver ----------------------------------
ODEinput.looky_nutrients         	= looky_nutrients;
ODEinput.looky_NO3                	= looky_NO3;
ODEinput.looky_NH4                  = looky_NH4;
ODEinput.looky_plgcNH4            	= looky_plgcNH4;
ODEinput.looky_bnthNH4            	= looky_bnthNH4;
ODEinput.looky_ANYPrimaryProducer	= looky_ANYPrimaryProducer;
ODEinput.looky_phytoplankton      	= looky_phytoplankton;
ODEinput.looky_macroalgae       	= looky_macroalgae;
ODEinput.looky_fleets             	= looky_fleets;
ODEinput.looky_ANYconsumer       	= looky_ANYconsumer;
ODEinput.looky_livingANDfleets      = looky_livingANDfleets;
ODEinput.looky_eggs             	= looky_eggs;
ODEinput.looky_terminalPLGCdetritus	= looky_terminalPLGCdetritus;
ODEinput.looky_terminalBNTHdetritus	= looky_terminalBNTHdetritus;
ODEinput.looky_ANYdetritus      	= looky_ANYdetritus;
ODEinput.looky_nonNO3             	= looky_nonNO3;

ODEinput.num_grps                   = num_grps;
ODEinput.num_nutrients              = num_nutrients;
ODEinput.num_NO3                    = num_NO3;
ODEinput.num_NH4                    = num_NH4;
ODEinput.num_plgcNH4              	= num_plgcNH4; % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
ODEinput.num_bnthNH4              	= num_bnthNH4; % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
ODEinput.num_ANYPrimaryProd     	= num_ANYPrimaryProd;
ODEinput.num_phytoplankton       	= num_phytoplankton;
ODEinput.num_macroalgae          	= num_macroalgae;
ODEinput.num_ANYconsumers           = num_ANYconsumers;
ODEinput.num_predators           	= num_predators;
ODEinput.num_livingANDfleets        = num_livingANDfleets;
ODEinput.num_eggs                 	= num_eggs;
ODEinput.num_ANYdetritus            = num_ANYdetritus;
% -------------------------------------------------------------------------


%% step 4e: matrix addresses of ECOTRAN generic model groups ---------------
%          NOTE: these must be entered manually for each new model, for easier plotting

%          CNP_11182021_P
rc_NO3                          = 1;
rc_plgcNH4                      = 2;
rc_bnthNH4                      = 3;
rc_diatoms                      = 4;
rc_diazotrophs                  = 5;
rc_OtherLrgPhytoplankton        = 6;
rc_SmlPhytoplankton             = 7;
rc_microzooplankton             = 8;
rc_copepods                     = 9;
rc_MesoCrustaceans              = 10;
rc_pteropods                    = 11;
rc_chaetognaths                 = 12;
rc_OtherMesozooplankton         = 13;
rc_euphausiids                  = 14;
rc_MicronektonCrustaceans       = 15;
rc_decapods                     = 16;
rc_FilterJellies                = 17;
rc_PredatorJellies              = 18;
rc_EpipelagicSquid              = 19;
rc_MesopelagicSquid             = 20;
rc_SmlEpipelagicFish            = 21;
rc_PlanktovoreMicronektonFish	= 22;
rc_PredatorMicronektonFish      = 23;
rc_BathypelagicFish             = 24;
rc_SnakeMackerel                = 25;
rc_OtherLrgPelagicFish          = 26;
rc_mahimahi                     = 27;
rc_lancetfish                   = 28;
rc_opah                         = 29;
rc_BigeyeTuna                   = 30;
rc_JuvBigeyeTuna                = 31;
rc_YellowfinTuna                = 32;
rc_JuvYellowfinTuna             = 33;
rc_AlbacoreTuna                 = 34;
rc_JuvAlbacoreTuna              = 35;
rc_SkipjackTuna                 = 36;
rc_JuvSkipjackTuna              = 37;
rc_StripedMarlin                = 38;
rc_BlueMarlin                   = 39;
rc_OtherBillfish                = 40;
rc_JuvBillfish                  = 41;
rc_BroadbillSwordfish         	= 42;
rc_BlueSharks                   = 43;
rc_OtherSharks                  = 44;
rc_SeaTurtles                   = 45;
rc_seabirds                     = 46;
rc_BaleenWhales                 = 47;
rc_odontocetes                  = 48;
rc_PelagicBacteria              = 49;
rc_BenthicBacteria              = 50;
rc_detritusA                    = 51;
rc_detritusB                    = 52; % <<--terminal pelagic detritus
rc_detritusC                    = 53; % <<--terminal benthic detritus
rc_fleetA                       = 54; % purse seine <<--note that fleets come AFTER terminal detritus in ECOTRAN-2
rc_fleetB                       = 55; % <<--note that fleets come AFTER terminal detritus in ECOTRAN-2
rc_fleetC                       = 56; % <<--note that fleets come AFTER terminal detritus in ECOTRAN-2
rc_fleetD                       = 57; % HI deep set longline <<--note that fleets come AFTER terminal detritus in ECOTRAN-2
rc_fleetE                       = 58; % <<--note that fleets come AFTER terminal detritus in ECOTRAN-2

rc_mesopelagics             = [rc_MesoCrustaceans, rc_pteropods, rc_MesopelagicSquid, rc_PlanktovoreMicronektonFish, rc_PredatorMicronektonFish, rc_SnakeMackerel, rc_lancetfish, rc_opah, rc_BlueSharks, rc_OtherSharks];
rc_MiningPlume_effectsGrp	= [rc_MesoCrustaceans, rc_pteropods, rc_euphausiids, rc_MicronektonCrustaceans, rc_decapods, rc_MesopelagicSquid, rc_PlanktovoreMicronektonFish, rc_PredatorMicronektonFish];

rc_bathypelagics    = [rc_BathypelagicFish];
rc_benthics         = [rc_BenthicBacteria];

rc_EPIpelagic       = 1;
rc_MESOpelagic      = 2;
rc_BATHYpelagic     = 3;
rc_BENTHIC          = 4;

% *************************************************************************





%% *************************************************************************
% STEP 5: define TransferEfficiency terms----------------------------------
%         SET MANUALLY: TransferEfficiency = 1 for all groups because we are
%                       working with the consumption budget matrix that tracks the fate of ALL consumption (not just predation)
%                       system losses (groups where TransferEfficiency < 1 include benthic detritus, fishery, (others?)
%         NOTE: terminal benthic detritus column in EnergyBudget sums to 1 and has NO OtherMortality, 
%               TE is the only means of removing terminal benthic detritus from system
TransferEfficiency                              = ones(num_MC, num_grps);	% re-initialize as rows of horizontal vector of ones
% *************************************************************************





%% *************************************************************************
% STEP 6: prepare physical parameters--------------------------------------
% step 6a: input time-frame info ------------------------------------------
%          NOTE: prepared BUI_ERD .csv file covers: 01-Jan-1969 through 30-Sep-2015 and has 6-hr resolution
%          NOTE: ERDcuti file covers:               28-Feb-1990 through 28-Feb-2022  and has 1-day resolution
%          NOTE: NH Line temperature covers:        05-May-1997 through 08-Jul-2021 and has 1-day resolution

switch switch_PhysicalModel
    case '2D_physics'
        
        disp('--->>> 2D physics')

        datestart                       = datenum('01-Jan-2001'); % SSS --> enter starting date
        dateend                         = datenum('31-Dec-2025'); % SSS --> enter ending date (default for dynamic runs tests ('31-Dec-2020'))
        dt                              = 3/24; % t-step; (days); (dt = 1 = 1 d; dt = 4/24 = 4 hours)
                                          % NOTE: take care to select good dt values for diel vertical migration 
                                          %       (other values do not scale well between 1 & -1 in sin diel cycle (probably due to rounding error of pi() function)
        datestart_OrdinalDate           = f_OrdinalDate(datestart);
        % dateend_OrdinalDate             = f_OrdinalDate(dateend); % QQQ dateend_OrdinalDate is not used
        min_t                           = datestart_OrdinalDate;
        max_t                           = (dateend - datestart) + datestart_OrdinalDate;
        t_grid                          = linspace(min_t, (max_t+1-dt), ((dateend - datestart + 1)/dt))'; % QQQ NEW VERSION!!!; (vertical vector); t_grid runs from datestart to dateend and is inclusive of dateend; intervals = dt
        num_t                           = length(t_grid); % length of t_grid; (scaler)
        PHYSICSinput.datestart          = datestart;
        PHYSICSinput.dateend            = dateend;
        PHYSICSinput.dt                 = dt;
        PHYSICSinput.t_grid             = t_grid;

        PHYSICSinput.smoothing_window	= 5; % window for smoothing before and after time-point (e.g., 2 = a window of 5 days centered on time-point t)
        % -------------------------------------------------------------------------


        % step 6b: prepare advection & mixing time-series for each model box ------
        ECOTRANphysics                  = f_ECOTRANphysics_CNP_06172022(PHYSICSinput);
        
        % p_plotECOTRANphysics_06242019(ECOTRANphysics)
        % -------------------------------------------------------------------------


        % SSS apply any desired change to the temperature time-series (DO NOT change temperature_reference)
        ECOTRANphysics.temperature_timeseries(:, 1)	= ECOTRANphysics.temperature_timeseries(:, 1)	 + 0; % (2D matrix: num_t X num_boxes)
        % -------------------------------------------------------------------------

        % SSS define Mining Plume conditions
        ECOTRANphysics.SSCreference_timeseries	= ones(num_t, ECOTRANphysics.num_boxes) * 0.01; % reference (background) Suspended Sediment Concentration; (mg/L); (2D matrix: num_t X num_boxes)
        ECOTRANphysics.SSC_timeseries           = ones(num_t, ECOTRANphysics.num_boxes) * 0.01; % local plume Suspended Sediment Concentration; (mg/L); (2D matrix: num_t X num_boxes)
        % -------------------------------------------------------------------------
       

        % step 6c: migration of each group (shared box face areas) ----------------
        %          NOTE: code does not handle migration flux uncertainty (nor physical flux uncertainty)
        ECOTRANmigration            = f_ECOTRANmigration_CNP_07212021(ECOTRANphysics); % area of overlap of neighboring model sub-domains
        % ODEinput.biomass_migrator = ECOTRANmigration.biomass_migrator;  % SSS special definition of boundary biomasses for migrators; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes); NOTE: de-comment migrator biomass lines within ODE code
        % -------------------------------------------------------------------------

        
        ADVECTION                       = ECOTRANphysics.ADVECTION;         % volume advected per day from source box to destination box; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
        HORIZONTALMIXING                = ECOTRANphysics.HORIZONTALMIXING;	% volume mixed horizontally per day from source box to destination box; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
        VERTICALMIXING                  = ECOTRANphysics.VERTICALMIXING;	% volume mixed vertically per day from source box to destination box; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
        SINKING                         = ECOTRANphysics.SINKING;           % box floor area between sinking source box and destination box; (m2); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
        MIGRATION                       = ECOTRANmigration.MIGRATION;       % migration as boundary area between boxes; (m2); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
        
        CompactFlux_ADVECTION           = f_CompactFluxTimeSeries_11182019(ADVECTION);        % compact ADVECTION 3D matrix
        CompactFlux_HORIZONTALMIXING	= f_CompactFluxTimeSeries_11182019(HORIZONTALMIXING); % compact HORIZONTALMIXING 3D matrix
        CompactFlux_VERTICALMIXING      = f_CompactFluxTimeSeries_11182019(VERTICALMIXING);   % compact VERTICALMIXING 3D matrix
        CompactFlux_SINKING             = f_CompactFluxTimeSeries_11182019(SINKING);          % compact SINKING 3D matrix
        CompactFlux_MIGRATION           = f_CompactFluxTimeSeries_11182019(MIGRATION);        % compact MIGRATION 3D matrix

	% end case 2D_physics -------------------------------------------------

    
    case '3D_ROMS'
        
        disp('--->>> 3D ROMS physics')

% % %         datestart                       = datenum('01-Jan-2008'); % SSS --> enter starting date
% % %         dateend                         = datenum('31-Dec-2009'); % SSS --> enter ending date (default for dynamic runs tests ('31-Dec-2020'))
% % %         dt                              = 24/24; % t-step; (days); (dt = 24/24 = 1 d; dt = 3/24 = 3 hours)
% % %                                           % NOTE: take care to select good dt values for diel vertical migration 
% % %                                           %       (other values do not scale well between 1 & -1 in sin diel cycle (probably due to rounding error of pi() function)
% % %         
% % %         datestart_OrdinalDate           = f_OrdinalDate(datestart);
% % %         % dateend_OrdinalDate             = f_OrdinalDate(dateend); % QQQ dateend_OrdinalDate is not used
% % %         min_t                           = datestart_OrdinalDate;
% % %         max_t                           = (dateend - datestart) + datestart_OrdinalDate;
% % %         t_grid                          = linspace(min_t, (max_t+1-dt), ((dateend - datestart + 1)/dt))'; % QQQ NEW VERSION!!!; (vertical vector); t_grid runs from datestart to dateend and is inclusive of dateend; intervals = dt
% % %         num_t                           = length(t_grid); % length of t_grid; (scaler)
% % %         PHYSICSinput.datestart          = datestart;
% % %         PHYSICSinput.dateend            = dateend;
% % %         PHYSICSinput.dt                 = dt;
% % %         PHYSICSinput.t_grid             = t_grid;
% % % 
% % %         PHYSICSinput.smoothing_window	= 5; % window for smoothing before and after time-point (e.g., 2 = a window of 5 days centered on time-point t)
% % %         % -------------------------------------------------------------------------
% % % 
% % % 
% % %         % step 6b: prepare advection & mixing time-series for each model box ------
% % %         % 3D ROMS driver
% % %         ECOTRANphysics                  = f_ECOTRANphysics_NCC_ROMS_05282022(PHYSICSinput); % use for NCC ROMS
% % % 
% % %         % p_plotECOTRANphysics_06242019(ECOTRANphysics)
% % %         % -------------------------------------------------------------------------
% % % 
% % % 
% % %         % step 6c: migration of each group (shared box face areas) ----------------
% % %         %          NOTE: code does not handle migration flux uncertainty (nor physical flux uncertainty)
% % %         ECOTRANmigration                = f_ECOTRANmigration_NCC_02222022(ECOTRANphysics);  % SSS; use for NCC ROMS
% % % 
% % %         % ODEinput.biomass_migrator       = ECOTRANmigration.biomass_migrator;  % SSS special definition of boundary biomasses for migrators; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes); NOTE: de-comment migrator biomass lines within ODE code
% % %         % -------------------------------------------------------------------------
% % % 
% % %         
% % %         CompactFlux_ADVECTION           = ECOTRANphysics.CompactFlux_ADVECTION; % (structure)
% % %         CompactFlux_HORIZONTALMIXING	= ECOTRANphysics.CompactFlux_HORIZONTALMIXING; % (structure)
% % %         CompactFlux_VERTICALMIXING      = ECOTRANphysics.CompactFlux_VERTICALMIXING; % (structure)
% % %         CompactFlux_SINKING             = ECOTRANphysics.CompactFlux_SINKING; % (structure); compact SINKING as box floor areas and connectivity information; apply functional group sinking speeds in ECOTRANdynamic_ code
% % %         CompactFlux_MIGRATION           = ECOTRANmigration.CompactFlux_MIGRATION; % (structure)
% % % 
% % %         num_boxes                       = ECOTRANphysics.num_boxes;
% % %         spatial_BiomassScalers         	= ones(1, num_boxes);	% NCC scalers for estimating initial (or mean) primary producer biomasses across model domain; NOTE: these values are assumed; NOTE: x2 in Box I used to compensate for 30m depth relative to 15 m depths in Boxes II & IV; FFF apply NPZD scalers here
        
	% end case 3D_ROMS ----------------------------------------------------

end % (switch switch_PhysicalModel)
% -------------------------------------------------------------------------


% step 6d: unpack & process physics variables -----------------------------
num_boxes                    	= ECOTRANphysics.num_boxes;
BoxVolume                    	= ECOTRANphysics.BoxVolume;     	% (m3); (2D matrix: num_t X num_boxes)
% BoxLength                    	= ECOTRANphysics.BoxLength;        	% (m); (2D matrix: num_t X num_boxes)
BoxHeight                     	= ECOTRANphysics.BoxHeight;       	% (m); (2D matrix: num_t X num_boxes)
% BoxWidth                     	= ECOTRANphysics.BoxWidth;       	% (m); (2D matrix: num_t X num_boxes)
% BoxFloorArea                    = ECOTRANphysics.BoxFloorArea;      % (m2); (2D matrix: num_t X num_boxes)

Io                            	= ECOTRANphysics.Io;             	% time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (vertical vector: num_t X 1)
current_light                   = ECOTRANphysics.current_light;     % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day; (vertical vector: num_t X 1)
sunrise                         = ECOTRANphysics.sunrise;           % time of sunrise (time in h from midnight; 12 = noon, set as a default)
sunset                          = ECOTRANphysics.sunset;            % time of sunset  (time in h from midnight; 12 = noon, set as a default)
Kw                            	= ECOTRANphysics.Kw;              	% Light attenuation_seawater; (scalar)
Kp                             	= ECOTRANphysics.Kp;              	% Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N); (scalar)
MLD                           	= ECOTRANphysics.MLD;            	% mixed-layer depth; (m); (vertical vector: num_t X 1)
EuphoticDepth                  	= repmat(MLD, [1 num_boxes]);     	% FFF (eventually move to f_physics code); depth of euphotic zone, used when converting the vertically-integrated EwE primary producer biomass to biomass/volume; depth; (m); (2D matrix: num_t X num_boxes)

NO3timeseries_conc              = ECOTRANphysics.NO3timeseries_conc; % NO3 + NO2 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
NH4timeseries_conc              = ECOTRANphysics.NH4timeseries_conc; % NH4 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
NO3initial_rate                 = ECOTRANphysics.NO3initial_rate;	 % initial NO3 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))
NH4initial_rate                 = ECOTRANphysics.NH4initial_rate;	 % initial NH4 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))

WWT_to_C                      	= ECOTRANphysics.WWT_to_C;          % (scalar)
atomic_mass_C                  	= ECOTRANphysics.atomic_mass_C;     % (scalar)
C_to_N_phytoplankton          	= ECOTRANphysics.C_to_N_phytoplankton; % (scalar)

spatial_BiomassScalers         	= [1 1 1 1];	% QQQ GoMex scalers for estimating initial (or mean) primary producer biomasses across model domain; NOTE: these values are assumed; NOTE: x2 in Box I used to compensate for 30m depth relative to 15 m depths in Boxes II & IV; FFF apply NPZD scalers here

grp_row                      	= 1:num_grps;
% -------------------------------------------------------------------------


%% step 6e: compact fluxes -------------------------------------------------
%           Remove information defining non-existing box links
%       CompactFlux
%           compact_flux            (2D matrix: num_t X num_fluxes)
%                                       NOTE: fluxes include all linked boxes +1 for external links
%           looky_flux              (2D matrix: num_fluxes X 3)
%                                       clm 1: (destiny box) = list of boxes importing water volume (+ boundary)
%                                       clm 2: (source box) = list of boxes exporting water volume (+ boundary)
%                                       clm 3: flux address in non-compacted flux 2D matrix: destiny X source
%                                       NOTE: fluxes include all linked boxes +1 for external links
%                                       NOTE: values constant independent of t
%           looky_boundary_import	(2D matrix: num_fluxes_BoundaryImport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume
%                                       clm 2: (source box) = identity of boxes exporting water volume (always the boundary flux number)
%                                       clm 3: (import flux address) = addresses of import flux clm in compact_flux)
%           looky_boundary_export   (2D matrix: num_fluxes_BoundaryExport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume (always the boundary flux number)
%                                       clm 2: (source box) = identity of boxes exporting water volume
%                                       clm 3: (export flux address) = addresses of export fluxes clm in compact_flux
%           unique_source           (vertical vector: list of source boxes (+ boundary))
%           unique_destiny          (vertical vector: list of destiny boxes (+ boundary))
%           num_fluxes              number of realized fluxes between boxes (or boundary) over full time-series
%           fname_CompactFlux       name of this FuncName_CompactFlux function


% ADVECTION -----
ODEinput.num_fluxes_advection                   = CompactFlux_ADVECTION.num_fluxes;
ODEinput.ADVECTION_compact                      = CompactFlux_ADVECTION.compact_flux;                         % (m3/d); (2D matrix: num_t X num_fluxes_advection)
ODEinput.looky_AdvectionFlux                    = CompactFlux_ADVECTION.looky_flux;                           % (2D matrix: num_fluxes_advection X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
ODEinput.looky_AdvectionBoundary_import         = CompactFlux_ADVECTION.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
ODEinput.looky_AdvectionBoundary_export         = CompactFlux_ADVECTION.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);

switch switch_ODEsolver
    case 'MatlabSolver'
        ODEinput.repmat_looky_ADVECTION_destiny         = repmat(CompactFlux_ADVECTION.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
        ODEinput.repmat_looky_ADVECTION_source          = repmat(CompactFlux_ADVECTION.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
        ODEinput.repmat_GrpRow_ADVECTION                = repmat(grp_row', [1, CompactFlux_ADVECTION.num_fluxes]);	% addressses; (2D matrix: num_grps X num_fluxes); NOTE transpose
    case 'CppSolver'
        repmat_GrpRow_ADVECTION                         = repmat(grp_row', [1, CompactFlux_ADVECTION.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes); NOTE transpose
        repmat_looky_ADVECTION_source                   = repmat(CompactFlux_ADVECTION.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
        repmat_looky_ADVECTION_destiny                  = repmat(CompactFlux_ADVECTION.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
        ODEinput.repmat_looky_ADVECTION_source          = [repmat_GrpRow_ADVECTION(:) repmat_looky_ADVECTION_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_advection) X 2)
        ODEinput.repmat_looky_ADVECTION_destiny         = [repmat_GrpRow_ADVECTION(:) repmat_looky_ADVECTION_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_advection) X 2)
end % (switch switch_ODEsolver)

% HORIZONTALMIXING -----
ODEinput.num_fluxes_HorizontalMixing            = CompactFlux_HORIZONTALMIXING.num_fluxes;
ODEinput.HORIZONTALMIXING_compact               = CompactFlux_HORIZONTALMIXING.compact_flux;                         % (m3/d); (2D matrix: num_t X num_fluxes_HorizontalMixing)
ODEinput.looky_HorizontalMixingFlux             = CompactFlux_HORIZONTALMIXING.looky_flux;                           % (2D matrix: num_fluxes_HorizontalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
ODEinput.looky_HorizontalMixingBoundary_import	= CompactFlux_HORIZONTALMIXING.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
ODEinput.looky_HorizontalMixingBoundary_export	= CompactFlux_HORIZONTALMIXING.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);

switch switch_ODEsolver
    case 'MatlabSolver'
        ODEinput.repmat_looky_HORIZONTALMIXING_destiny  = repmat(CompactFlux_HORIZONTALMIXING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
        ODEinput.repmat_looky_HORIZONTALMIXING_source	= repmat(CompactFlux_HORIZONTALMIXING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
        ODEinput.repmat_GrpRow_HORIZONTALMIXING         = repmat(grp_row', [1, CompactFlux_HORIZONTALMIXING.num_fluxes]);	% addressses; (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
    case 'CppSolver'
        repmat_GrpRow_HORIZONTALMIXING                	= repmat(grp_row', [1, CompactFlux_HORIZONTALMIXING.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
        repmat_looky_HORIZONTALMIXING_source          	= repmat(CompactFlux_HORIZONTALMIXING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
        repmat_looky_HORIZONTALMIXING_destiny         	= repmat(CompactFlux_HORIZONTALMIXING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
        ODEinput.repmat_looky_HORIZONTALMIXING_source	= [repmat_GrpRow_HORIZONTALMIXING(:) repmat_looky_HORIZONTALMIXING_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_HorizontalMixing) X 2)
        ODEinput.repmat_looky_HORIZONTALMIXING_destiny	= [repmat_GrpRow_HORIZONTALMIXING(:) repmat_looky_HORIZONTALMIXING_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_HorizontalMixing) X 2)
end % (switch switch_ODEsolver)

% VERTICALMIXING -----
ODEinput.num_fluxes_VerticalMixing              = CompactFlux_VERTICALMIXING.num_fluxes;
ODEinput.VERTICALMIXING_compact                 = CompactFlux_VERTICALMIXING.compact_flux;                         % (m3/d); (2D matrix: num_t X num_fluxes_VerticalMixing)
ODEinput.looky_VerticalMixingFlux               = CompactFlux_VERTICALMIXING.looky_flux;                           % (2D matrix: num_fluxes_VerticalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
ODEinput.looky_VerticalMixingBoundary_import	= CompactFlux_VERTICALMIXING.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
ODEinput.looky_VerticalMixingBoundary_export	= CompactFlux_VERTICALMIXING.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);

switch switch_ODEsolver
    case 'MatlabSolver' 
        ODEinput.repmat_looky_VERTICALMIXING_destiny	= repmat(CompactFlux_VERTICALMIXING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
        ODEinput.repmat_looky_VERTICALMIXING_source     = repmat(CompactFlux_VERTICALMIXING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
        ODEinput.repmat_GrpRow_VERTICALMIXING           = repmat(grp_row', [1, CompactFlux_VERTICALMIXING.num_fluxes]);	% addressses; (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
    case 'CppSolver'
        repmat_GrpRow_VERTICALMIXING                 	= repmat(grp_row', [1, CompactFlux_VERTICALMIXING.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
        repmat_looky_VERTICALMIXING_source          	= repmat(CompactFlux_VERTICALMIXING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
        repmat_looky_VERTICALMIXING_destiny         	= repmat(CompactFlux_VERTICALMIXING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
        ODEinput.repmat_looky_VERTICALMIXING_source  	= [repmat_GrpRow_VERTICALMIXING(:) repmat_looky_VERTICALMIXING_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_VerticalMixing) X 2)
        ODEinput.repmat_looky_VERTICALMIXING_destiny	= [repmat_GrpRow_VERTICALMIXING(:) repmat_looky_VERTICALMIXING_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_VerticalMixing) X 2)
end % (switch switch_ODEsolver)

% SINKING -----
ODEinput.num_fluxes_sinking                     = CompactFlux_SINKING.num_fluxes;
num_fluxes_sinking                              = CompactFlux_SINKING.num_fluxes;
SinkingArea_compact                             = CompactFlux_SINKING.compact_flux;                         % box floor area between sinking SOURCE box and DESTINY box; (m2); (2D matrix: num_t X num_fluxes_sinking)
SinkingArea_compact                             = reshape(SinkingArea_compact, [num_t, 1, num_fluxes_sinking]); % (m2); (3D matrix: num_t X 1 X num_fluxes_sinking)
SinkingArea_compact                             = repmat(SinkingArea_compact, [1, num_grps, 1]);	% (m2); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% ODEinput.SINKING_compact <<--sinking as volumetric flux (m3/d) for each functional group is calculated below
ODEinput.looky_SinkingFlux                      = CompactFlux_SINKING.looky_flux;                           % (2D matrix: num_fluxes_sinking X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
ODEinput.looky_SinkingBoundary_import           = CompactFlux_SINKING.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
ODEinput.looky_SinkingBoundary_export           = CompactFlux_SINKING.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);

switch switch_ODEsolver
    case 'MatlabSolver'
        ODEinput.repmat_looky_SINKING_destiny           = repmat(CompactFlux_SINKING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
        ODEinput.repmat_looky_SINKING_source            = repmat(CompactFlux_SINKING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
        ODEinput.repmat_GrpRow_SINKING                  = repmat(grp_row', [1, num_fluxes_sinking]);        % addressses; (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
    case 'CppSolver'
        repmat_GrpRow_SINKING                           = repmat(grp_row', [1, CompactFlux_SINKING.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
        repmat_looky_SINKING_source                     = repmat(CompactFlux_SINKING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
        repmat_looky_SINKING_destiny                    = repmat(CompactFlux_SINKING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
        ODEinput.repmat_looky_SINKING_source            = [repmat_GrpRow_SINKING(:) repmat_looky_SINKING_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_sinking) X 2)
        ODEinput.repmat_looky_SINKING_destiny           = [repmat_GrpRow_SINKING(:) repmat_looky_SINKING_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_sinking) X 2)
end % (switch switch_ODEsolver)

% MIGRATION -----
ODEinput.num_fluxes_migration                 	= CompactFlux_MIGRATION.num_fluxes;                           % migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
num_fluxes_migration                            = CompactFlux_MIGRATION.num_fluxes;                           % migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
MigrationArea_compact                           = CompactFlux_MIGRATION.compact_flux;                         % migration as boundary area between boxes; (m2); (3D matrix: num_t X num_fluxes_migration)
% ODEinput.MIGRATION_compact <<--migration as volumetric flux (m3/d) for each functional group is calculated below
looky_MigrationFlux                             = CompactFlux_MIGRATION.looky_flux;                           % (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
ODEinput.looky_MigrationFlux                    = CompactFlux_MIGRATION.looky_flux;                           % (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
ODEinput.looky_MigrationBoundary_import         = CompactFlux_MIGRATION.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
ODEinput.looky_MigrationBoundary_export         = CompactFlux_MIGRATION.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);

switch switch_ODEsolver
    case 'MatlabSolver'
        ODEinput.repmat_looky_MIGRATION_destiny         = repmat(CompactFlux_MIGRATION.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
        ODEinput.repmat_looky_MIGRATION_source          = repmat(CompactFlux_MIGRATION.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
        ODEinput.repmat_GrpRow_MIGRATION                = repmat(grp_row', [1, num_fluxes_migration]);      % addressses; (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
    case 'CppSolver'
        repmat_GrpRow_MIGRATION                         = repmat(grp_row', [1, CompactFlux_MIGRATION.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
        repmat_looky_MIGRATION_source                   = repmat(CompactFlux_MIGRATION.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
        repmat_looky_MIGRATION_destiny                  = repmat(CompactFlux_MIGRATION.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
        ODEinput.repmat_looky_MIGRATION_source          = [repmat_GrpRow_MIGRATION(:) repmat_looky_MIGRATION_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_migration) X 2)
        ODEinput.repmat_looky_MIGRATION_destiny         = [repmat_GrpRow_MIGRATION(:) repmat_looky_MIGRATION_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_migration) X 2)
end % (switch switch_ODEsolver)
% -------------------------------------------------------------------------


%% step 6f: sinking rate of each group (m/d) -------------------------------
%          NOTE: sinking is treated like a physical term (i.e., not incorporated into EnergyBudget)
%          NOTE: apply this factor whether or not using Michaelis-Menten for primary producers
SinkingSpeed                        	= zeros(num_t, num_grps, num_fluxes_sinking);	% initialze sinking speed time-series; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)

SinkingSpeed(:, rc_detritusA, :)      	= repmat((5),    [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
SinkingSpeed(:, rc_detritusB, :)       	= repmat((200),  [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
SinkingSpeed(:, rc_detritusC, :)       	= repmat((1500), [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)

% MichaelisMenten_w                               = [0.6 1.0];                                    % SSS; sinking speed; (m/d); [Sm Phytoplankton, Lg Phytoplankton]
% SinkingSpeed(:, looky_ANYPrimaryProducer, :)	= repmat(MichaelisMenten_w, [num_t, 1, num_fluxes_sinking]);	% sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% -------------------------------------------------------------------------


%% step 6g: define migration speeds & duration for each group (m/d) --------
%          NOTE: dusk phase are negative (in "standard" DVM)
DVMinput.MigrationArea_compact          = MigrationArea_compact; % migration as boundary area between boxes; (m2); (3D matrix: num_t X num_fluxes_migration)
DVMinput.MigrationSpeed                 = zeros(num_grps, (num_boxes+1), (num_boxes+1)); % initialize; (m/d); (3D matrix: num_grps, source (num_boxes+1) X destiny (num_boxes+1))

DVMspeed_meso2epi                       = zeros(num_grps, 2); % initialize; DVM speed MESO<<-->>EPI; (2D matrix: num_grps X 2 [dusk dawn]);
DVMspeed_bathy2meso                     = zeros(num_grps, 2); % initialize; DVM speed BATHY<<-->>MESO; (2D matrix: num_grps X 2 [dusk dawn]);

% DVMspeeds 
%       NOTE: determination of DVM speed calculations are done outside of
%       ECOTRAN, based on the mean day & night vertical distribution and the fraction of the biomass pool
%       that undergoes migration, using code ECOTRANdynamic_DVMcalibration_08202021.m & p_plotMigrationFraction_08232021.m
%       Excel file DVMcalibration_CNP_08112021.xlsx "summary" tab has the
%       final DVM speed parameters to provide here. For each functional group, 
%       the first value is the dusk migration speed and the second is the
%       dawn migration speed.
DVMspeed_meso2epi(rc_copepods, 1:2)                  	= [313 27]; % (m/d)
DVMspeed_meso2epi(rc_MesoCrustaceans, 1:2)            	= [946 441]; % (m/d)
DVMspeed_meso2epi(rc_pteropods, 1:2)                 	= [2753 469]; % (m/d)
DVMspeed_meso2epi(rc_chaetognaths, 1:2)               	= [235 25]; % (m/d)
DVMspeed_meso2epi(rc_OtherMesozooplankton, 1:2)      	= [235 44]; % (m/d)
DVMspeed_meso2epi(rc_euphausiids, 1:2)               	= [3091 242]; % (m/d)
DVMspeed_meso2epi(rc_MicronektonCrustaceans, 1:2)    	= [1709 178]; % (m/d)
DVMspeed_meso2epi(rc_decapods, 1:2)                     = [2312 212]; % (m/d)
DVMspeed_meso2epi(rc_MesopelagicSquid, 1:2)         	= [1708 1378]; % (m/d)
DVMspeed_meso2epi(rc_PlanktovoreMicronektonFish, 1:2)	= [666 1525]; % (m/d)
DVMspeed_meso2epi(rc_SnakeMackerel, 1:2)             	= [1708 251]; % (m/d)
DVMspeed_meso2epi(rc_opah, 1:2)                     	= [235 620]; % (m/d)
DVMspeed_meso2epi(rc_BigeyeTuna, 1:2)                 	= [7852 1460]; % (m/d)
DVMspeed_meso2epi(rc_JuvBigeyeTuna, 1:2)              	= [7855 890]; % (m/d)
DVMspeed_meso2epi(rc_AlbacoreTuna, 1:2)              	= [7855 786]; % (m/d)
DVMspeed_meso2epi(rc_JuvAlbacoreTuna, 1:2)          	= [7855 399]; % (m/d)

% MESO<<-->>EPI migration
DVMinput.MigrationSpeed(:, 1, 2)           = DVMspeed_meso2epi(:, 2); % dawn migration, (epi-->>meso)
DVMinput.MigrationSpeed(:, 2, 1)           = DVMspeed_meso2epi(:, 1) * (-1); % dusk migration, (meso-->>epi), noon to midnight; NOTE: should be negative

% BATHY<<-->>MESO migration
DVMinput.MigrationSpeed(:, 2, 3)           = DVMspeed_bathy2meso(:, 2); % dawn migration, (meso-->>bathy)
DVMinput.MigrationSpeed(:, 3, 2)           = DVMspeed_bathy2meso(:, 1) * (-1); % dusk migration, (bathy-->>meso), noon to midnight; NOTE: should be negative


DVM     	= f_DVMsinusoid_08122021(ECOTRANphysics, ODEinput, DVMinput, t_grid); % calculate Diel Vertical Migration terms; (structure)
% -------------------------------------------------------------------------


%% step 6i: define production loss fractions -------------------------------
%          NOTE: this used for initial conditions only (for dynamic models)
%          NOTE: the ODE accounts for physical removal (or addition) at each time-step; while the static model accounts for physical loss as a reduction in Transfer Efficiency
%                this is the way it was handled in John's original code outline and in the pubs
PhysicalLossFraction        = repmat(ProductionLossScaler', [num_t, 1]);	% (2D matrix: num_t X num_grps); NOTE transpose
PhysicalLossFraction        = PhysicalLossFraction * 0;                     % SSS set to 0 for time-dynamic runs
% -------------------------------------------------------------------------


%% step 6j: pack physics values for ODE into ODEinput ----------------------
ODEinput.num_boxes                    	= num_boxes;
% ODEinput.BoxLength                    	= BoxLength;                         % (m);  (2D matrix: time X num_boxes)
% ODEinput.BoxHeight                    	= BoxHeight;                         % (m);  (2D matrix: time X num_boxes)
% ODEinput.BoxWidth                      	= BoxWidth;                          % (m);  (2D matrix: time X num_boxes)
ODEinput.BoxVolume                  	= BoxVolume;                         % (m3); (2D matrix: time X num_boxes)
ODEinput.t_grid                      	= t_grid;
ODEinput.dt                             = dt;
ODEinput.num_t                          = num_t;
ODEinput.MLD                        	= MLD;                               % (vertical vector: length = time)
ODEinput.Io                           	= Io;
ODEinput.Kw                           	= Kw;
ODEinput.Kp                           	= Kp;
ODEinput.RetentionScaler              	= repmat(RetentionScaler, [1, (num_boxes)]);	% (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
ODEinput.SinkingSpeed               	= SinkingSpeed;                      % sinking speed; (m/d); (3D matrix: num_t X num_grps X num_boxes)
ODEinput.flux_domain_import_t        	= zeros(num_grps, (num_boxes+1));	 % initialized variable used in intraODE; (2D matrix: num_grps X num_boxes+1)
ODEinput.flux_domain_export_t       	= zeros(num_grps, (num_boxes+1));	 % initialized variable used in intraODE; (2D matrix: num_grps X num_boxes+1)
ODEinput.flux_domain_import_driver_t	= zeros(num_grps, num_boxes);        % initialized variable used in intraODE; (2D matrix: num_grps X num_boxes)
ODEinput.biomass_plus1                  = zeros(num_grps, 1, (num_boxes+1)); % initialized as all zeros; add 1 clm for boundary fluxes; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
ODEinput.PhysicalLossFraction        	= PhysicalLossFraction;              % (2D matrix: num_t X num_grps)
ODEinput.SINKING_compact             	= SinkingArea_compact .* SinkingSpeed;  % sinking fluxes; (m3/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)                    % box floor area between sinking source box and destination box; (m2); (2D matrix: num_t X num_fluxes_sinking)
ODEinput.MIGRATION_compact              = DVM.MIGRATION_compact;             	% migration fluxes; (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
% *************************************************************************





%% *************************************************************************
% % STEP QQQ: calculate temperature response scalers-------------------------
% %
% % Q10_scaler = ones(num_t, num_grps, num_boxes); % place-holder default = 1; (3D matrix: num_t X num_grps X num_boxes)
% % QQQ need to load the Q10 parameter each groups
% % 
% % Q10 shape (horizontal vector: 1 X num_grps)
% % T1 & T2 shape (2D matrices: time X sub-region_box)

% step QQQa: define Q10 parameters ------------------------------------------
Q10base = ones(1, num_grps); % (horizontal vector: 1 X num_grps)

Q10base(1, 8:45) = 2.3; % apply default value to all living, non-warm-blooded groups


%          LESS than trophic level 3 = [8 9:11 13 14 17];
%          MORE than trophic level 3     = [12 15 16 18 19 20 21 22 23 24 27 32 30 34 36];

% Migrators        
        Q10base(rc_copepods)                             = 2;    %9  TL=2.74
        Q10base(rc_MesoCrustaceans)                      = 2.2;  %10 TL=2.76
        Q10base(rc_pteropods)                            = 2;    %11 TL=2.68
        Q10base(rc_chaetognaths)                         = 3.5;  %12 TL=3.07
        Q10base(rc_OtherMesozooplankton)                 = 1.5;  %13  TL=2.76
        Q10base(rc_euphausiids)                          = 1.7;  %14 TL=2.91
        Q10base(rc_MicronektonCrustaceans)               = 1.55; %15 TL=3.30
        Q10base(rc_decapods)                             = 2;    %16 TL=3.20
        Q10base(rc_MesopelagicSquid)                     = 2.55; %20 TL=3.49
        Q10base(rc_PlanktovoreMicronektonFish)           = 2.1;  %22 TL=3.50                              
        Q10base(rc_BigeyeTuna)                           = 1.88; %30 TL=5.05            
        Q10base(rc_AlbacoreTuna)                         = 1.67; %34 TL=4.91

% Non-Migrators
        Q10base(rc_microzooplankton)                    = 1.5;  %8   TL=2.25  
        Q10base(rc_FilterJellies)                     	= 3.45; %17 TL=2.21
        Q10base(rc_PredatorJellies)                     = 1.67; %18 TL=3.72
        Q10base(rc_EpipelagicSquid)                     = 2.1;  %19  TL=4.29
        Q10base(rc_SmlEpipelagicFish)                   = 1.67; %21 TL=3.78
        Q10base(rc_PredatorMicronektonFish)             = 2.1;  %23  TL=4.34
        Q10base(rc_BathypelagicFish)                    = 1.9;  %24  TL=3.87
        Q10base(rc_mahimahi)                            = 2.86; %27 TL=4.91
        Q10base(rc_YellowfinTuna)                       = 3.16; %32 TL=4.83
        Q10base(rc_SkipjackTuna)                        = 2.44; %36 TL=4.56


% % % QQQ exlude TL groups
low_TL      = [1:18, 20 22 49:53];
high_TL     = [19 21 23:48 54:58];
Q10base(1, low_TL) = 1;


        
% step QQQb: define reference temperatures
temp_temperature_reference                      = ECOTRANphysics.temperature_reference'; % (vertical vector: num_boxes X 1); NOTE transpose

% temp_temperature_reference(rc_EPIpelagic) = 20; %change the epipelagic temp here

temp_temperature_reference                      = repmat(temp_temperature_reference, [1, num_grps]); % (2D matrix: num_boxes X num_grps)
build_temperature_reference                     = temp_temperature_reference(rc_EPIpelagic, :); % by default, the EPIpelagic is the reference temperature for all groups; (horizontal vector: 1 X num_grps)
build_temperature_reference(rc_mesopelagics)    = temp_temperature_reference(rc_MESOpelagic, rc_mesopelagics); % use the mesopelagic as the reference temperature for groups with >=50% in the mesopelagic; (horizontal vector: 1 X num_grps)
build_temperature_reference(rc_bathypelagics)   = temp_temperature_reference(rc_BATHYpelagic, rc_bathypelagics); % use the bathypelagic as the reference temperature for groups with >=50% in the bathypealgic; (horizontal vector: 1 X num_grps)
build_temperature_reference(rc_benthics)        = temp_temperature_reference(rc_BENTHIC, rc_benthics); % use the benthic as the reference temperature for groups with >=50% in the benthic; (horizontal vector: 1 X num_grps)
ECOTRANphysics.temperature_reference            = build_temperature_reference; % (horizontal vector: 1 X num_grps)

% % % ECOTRANphysics.temperature_reference(rc_MesopelagicSquid) = 24;
% % % ECOTRANphysics.temperature_timeseries(:, rc_EPIpelagic) = ECOTRANphysics.temperature_timeseries(:, rc_EPIpelagic) + 10; % change EPIpelagic temperature



if length(Q10base) ~= num_grps
    warning('The number of Q10 parameters does not equal the number of groups')
end

looky_undefined     = find(Q10base == 1);

switch switch_physiology
    case 'physiology_ON'
        disp('Calculating Q10 scaler')
        
        switch switch_physiology_scenario
            case 'base'
                Q10_parameter = Q10base;
            % end (case 'base') --------------

            case 'up10perc' 
                disp('Q10 scaled UP 110%')
                Q10_parameter = Q10base * 1.1;
                Q10_parameter(looky_undefined) = 1; % set undefined Q10 values back to 1
            % end (case 'up10perc') ----------

             case 'up25perc' 
                disp('Q10 scaled UP 125%')
                Q10_parameter = Q10base * 1.25;
                Q10_parameter(looky_undefined) = 1; % set undefined Q10 values back to 1
            % end (case 'up25perc') ----------

            case 'up40perc' 
                disp('Q10 scaled UP 140%')
                Q10_parameter = Q10base * 1.4;
                Q10_parameter(looky_undefined)
            % end (case 'up40perc') ----------

            case 'dwn10perc'
                disp('Q10 scaled DOWN to 90%')
                Q10_parameter = Q10base * 0.9;
                Q10_parameter(looky_undefined) = 1; % set undefined Q10 values back to 1
            % end (case 'dwn10perc') ---------

            case 'dwn25perc'
                disp('Q10 scaled DOWN to 75%')
                Q10_parameter = Q10base * 0.75;
                Q10_parameter(looky_undefined) = 1; % set undefined Q10 values back to 1
            % end (case 'dwn10perc') ---------

        end % (switch phys_scenario)        
        
        [Q10_scaler, fname_physiology_Q10] = f_physiology_Q10_06172022(ECOTRANphysics, Q10_parameter); % (3D matrix: num_t X num_grps X num_boxes)

	% end (case 'physiology_ON') ---------
        
    case 'physiology_OFF'
        disp('No Q10 scaling is used.')
        Q10_scaler = [];
	% end (case 'physiology_OFF') --------
        
end % (switch switch_physiology) -----------


% step qqq: Thornton & Lessem temperature adjustment model for consumption
%           consumption scaling factor (relative to q_max); (3D matrix: num_t X num_grps X num_boxes)
ThorntonLessem_parameters                       = ones(8, num_grps); %  initialize (2D matrix: 8 X num_grps)

% % define ThorntonLessem_parameters for individual groups 
% ThorntonLessem_parameters(1, rc_herring)        = 1;    % te1; Temperature for xk1 (in ºC) 
% ThorntonLessem_parameters(2, rc_herring)        = 13;   % te2; Temperature for xk2 (in ºC)
% ThorntonLessem_parameters(3, rc_herring)        = 15;   % te3; Temperature for xk3 (in ºC) 
% ThorntonLessem_parameters(4, rc_herring)        = 23;   % te4; Temperature for xk4 (in ºC)
% ThorntonLessem_parameters(5, rc_herring)        = 0.10; % xk1; Proportion of CMAX at te1      
% ThorntonLessem_parameters(6, rc_herring)        = 0.98; % xk2; Proportion of CMAX at te2      
% ThorntonLessem_parameters(7, rc_herring)        = 0.98; % xk3; Proportion of CMAX at te3      
% ThorntonLessem_parameters(8, rc_herring)        = 0.01; % xk4; Proportion of CMAX at te4      


switch switch_ThorntonLessem_physiology
    case 'ThorntonLessem_ON'
        disp('Calculating Thornton-Lessem consumption rate scaler')
        [q_TemperatureScaler, fname_ThorntonLessem]     = f_ThorntonLessem_05262022(ECOTRANphysics, ThorntonLessem_parameters); % value between 0 and 1; (3D matrix: num_t X num_grps X num_boxes)
	% end (case 'physiology_ON') ---------
        
    case 'ThorntonLessem_OFF'
        disp('Physiology Thornton-Lessem is OFF')
        fname_ThorntonLessem	= 'f_ThorntonLessem not used';
        q_TemperatureScaler     = ones(num_t, num_grps, num_boxes); 
        
	% end (case 'physiology_OFF') --------
end % (switch switch_physiology) -----------


% step qqq: Mining Plume adjustments to consumption (q) and senescence ----
PlumeEffect_parameters_feeding = zeros(1, num_grps); % initialize; (horizontal vector: 1 X num_grps)
PlumeEffect_parameters_feeding(1, [rc_MiningPlume_effectsGrp]) = -0.532; % plume effect on feeding; response slope

PlumeEffect_parameters_survival = zeros(1, num_grps); % initialize; (horizontal vector: 1 X num_grps)
PlumeEffect_parameters_survival(1, [rc_MiningPlume_effectsGrp]) = -0.128; % plume effect on survival; response slope

switch switch_MiningPlume_physiology
    
    case 'MiningPlume_ON'
        disp('Mining Plume effect is ON')
        [q_PlumeEffect_scaler, fname_physiology_MiningPlume]     	= f_physiology_MiningPlume_08212022(ECOTRANphysics, PlumeEffect_parameters_feeding); % (3D matrix: num_t X num_grps X num_boxes)
        [survival_PlumeEffect_scaler, fname_physiology_MiningPlume]	= f_physiology_MiningPlume_08212022(ECOTRANphysics, PlumeEffect_parameters_survival); % (3D matrix: num_t X num_grps X num_boxes)
        mortality_PlumeEffect_scaler                                = 1 + survival_PlumeEffect_scaler; % NOTE: survival_PlumeEffect_scaler should be negative
	% end (case 'MiningPlume_ON') ---------
        
    case 'MiningPlume_OFF'
        disp('Mining Plume effect is OFF')
        fname_physiology_MiningPlume	= 'f_physiology_MiningPlume not used';
        q_PlumeEffect_scaler            = ones(num_t, num_grps, num_boxes);
        mortality_PlumeEffect_scaler 	= ones(num_t, num_grps, num_boxes);
	% end (case 'MiningPlume_OFF') --------
    
end % (switch switch_MiningPlume_physiology) -----------

% *************************************************************************





%% *************************************************************************
% STEP 7: prepare external_driver time-series------------------------------
%         NOTE: This is the external input that enters the model domain via 
%               advection & mixing.
%         NOTE: DEFAULT: The reflective boundary assumption is that biomass of non-external_driver 
%               groups are the same on either side of model domain outer boundaries.

% step 7a: define external driver group(s) & driver time-series -----------
%          NOTE: driver time-series (external_driver) is a biomass density and not a rate (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
%          FFF: defining for all boxes now, but will try to trim down to just the boxes with fluxes

looky_driver                            = [looky_NO3; looky_plgcNH4]; % SSS; identify driver group(s)
num_drivers                             = length(looky_driver);
NO3timeseries_conc                      = reshape(NO3timeseries_conc, [num_t, 1, num_boxes]);	% (mmole N/m3); (3D matrix: num_t X 1 X num_boxes)
NH4timeseries_conc                      = reshape(NH4timeseries_conc, [num_t, 1, num_boxes]);	% (mmole N/m3); (3D matrix: num_t X 1 X num_boxes)

external_driver(:, looky_NO3, :)     	= NO3timeseries_conc; % SSS; external boundary driver (e.g. NO3) biomass for each box; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes)
external_driver(:, looky_plgcNH4, :)	= NH4timeseries_conc; % SSS; external boundary driver (e.g. NH4) biomass for each box; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes)

external_driver(:, :, (end+1))          = 0;              	    % external boundary driver (e.g. NO3) biomass for each box; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
% -------------------------------------------------------------------------


%% step 7b: define boundary biomass time-series ----------------------------
%             NOTE: boundary conditions are defined for each domain box
%                   whether or not that box is on the edge of the domain 
%                   (if not, the boundary values are not used)
%	OPTION 1: reflective boundary
%             NOTE: biomass is imported into each boundary box from external environment of same biomass densities
%             NOTE: boundary conditions defined by external_driver & by box biomass @ t within the ODE
%
% %   OPTION 2: defined boundary conditions
% %             NOTE: manually define a time-series of conditions
% %             NOTE: could define via function, via load external file, etc...
% %             NOTE: could update and use f_StaticProductionTimeseries_09042017(ODEinput, ECOTRANphysics, production_input, looky_OffshoreSurfaceBox)
%     biomass_boundary       	= zeros(num_t, num_grps, (num_boxes+1));	% initialize; (mmole N/m3); (3D matrix: num_t X num_grps X num_boxes+1)
% -------------------------------------------------------------------------


%% step 7c: pack external_driver & biomass_boundary into ODEinput ----------
ODEinput.looky_driver        	= looky_driver;	% driver group address(es)
ODEinput.num_drivers            = num_drivers;
ODEinput.external_driver     	= external_driver;      % (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
% ODEinput.biomass_boundary       = biomass_boundary;	  % NOTE: use for defined boundary conditions (OPTION 2); (mmole N/m3); (3D matrix: num_t X num_grps X num_boxes+1)
% *************************************************************************





%% *************************************************************************
% STEP 8: define functional predator-prey relations------------------------
%         select a "linear" or "ECOSIM" response --------------------------
%         NOTE: FunctionalResponseParams is a function of the CONSUMER (not the producer) and needs to be aligned with the consumer ROWS in ECOTRAN (not producer clms)
%         NOTE: to change one half-sat constant for individual groups, change FunctionalResponseParams([looky_grp])
%                FunctionalResponseParams = 0 is "constant donor-driven"
%                FunctionalResponseParams = 1 is "non-linear" ECOSIM default
FunctionalResponse_CV       = ones(num_grps, 4) * 0;        % SSS --> set uncertainty level for functional response terms; (2D matrix: num_grps->CONSUMERS X 4)
[FunctionalResponseParams, fname_FunctionalResponse_MonteCarlo]	= f_FunctionalResponse_MonteCarlo_05132021(ECOTRAN, ODEinput, FunctionalResponse_CV); % producer "vulnerabilities", m_p; (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)

switch switch_FunctionalResponse

    case 'Linear' % constant (predation independent of predator biomass)
        FunctionalResponseParams	= FunctionalResponseParams * 0;     % SSS force to linear for testing & default; (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)
        disp('NOTE: DONOR-DRIVEN FUNCTIONAL RESPONSE')
	% end (case 'Linear') ----------
    
    case 'NonLinear_default' % non-linear default
        disp('NOTE: DEFAULT NON-LINEAR FUNCTIONAL RESPONSE')
	% end (case 'NonLinear_default') ----------
    
    case 'NonLinear_alt' % constant (predation independent of predator biomass)
        FunctionalResponseParams	= FunctionalResponseParams * 2;     % SSS alternate (steeper) Functional Response; (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)
        disp('NOTE: ALTERNATE NON-LINEAR FUNCTIONAL RESPONSE')
	% end (case 'NonLinear_alt') ----------
    
end % switch (switch_FunctionalResponse)-----------------------------------
% *************************************************************************





%% *************************************************************************
% STEP 9: define Michaelis-Menten functional predator-prey relations-------
%         NOTE: for primary producers

% step 9a: initialize Michaelis-Menten parameters -------------------------
MichaelisMenten_Vmax     = zeros(num_ANYPrimaryProd, 1, num_MC); % Vmax  = maximum nutrient uptake rate;          (1/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_KNO3     = zeros(num_ANYPrimaryProd, 1, num_MC); % K     = NO3 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_KNH4     = zeros(num_ANYPrimaryProd, 1, num_MC); % K     = NH4 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_alpha    = zeros(num_ANYPrimaryProd, 1, num_MC); % alpha = initial slope of light response curve; (m2/W/d);     (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_psi      = zeros(num_ANYPrimaryProd, 1, num_MC); % psi   = NO3 uptake inhibition by NH4;          (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_w        = zeros(num_ANYPrimaryProd, 1, num_MC); % w     = sinking rate;                          (m/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_eta      = zeros(num_ANYPrimaryProd, 1, num_MC); % eta   = non-grazing mortality;                 (1/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% -------------------------------------------------------------------------


% %% step 9b: assign values --------------------------------------------------
% %	       NOTE: allow no nutrient uptake in sub-surface boxes; MichaelisMenten_Vmax(:, :, looky_SubSurfaceBoxes) = 0;
% %	       NOTE: for this example, row 1 = small phyto, row 2 = large phyto, row 3 = macroalgae
% MichaelisMenten_Vmax(1, 1, :)	= 0.42;     % Vmax  = maximum nutrient uptake rate;          (1/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_KNO3(1, 1, :)	= 0.14;     % K     = NO3 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_KNH4(1, 1, :)	= 0.14;     % K     = NH4 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_alpha(1, 1, :)	= 0.025;    % alpha = initial slope of light response curve; (m2/W/d);     (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_psi(1, 1, :)	= 0;        % psi   = NO3 uptake inhibition by NH4;          (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_w(1, 1, :)      = 0.6;      % NOTE: ALREADY DEFINED & APPLIED IN STEP 6e as SinkingSpeed; w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_eta(1, 1, :)	= 0.11;     % QQQ made up value for testing; eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)

% MichaelisMenten_Vmax(2, 1, :)	= 1.5;      % Vmax  = maximum nutrient uptake rate;          (1/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_KNO3(2, 1, :)	= 1.2;      % K     = NO3 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_KNH4(2, 1, :)	= 1.2;      % K     = NH4 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_alpha(2, 1, :)	= 0.0375;   % alpha = initial slope of light response curve; (m2/W/d);     (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_psi(2, 1, :)	= 1.46;     % psi   = NO3 uptake inhibition by NH4;          (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_w(2, 1, :)      = 1;        % NOTE: ALREADY DEFINED & APPLIED IN STEP 6e as SinkingSpeed; w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% MichaelisMenten_eta(2, 1, :)	= 0.12;     % QQQ made up value for testing; eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% *************************************************************************





%% *************************************************************************
% STEP 10: adjustments to ConsumptionBudget--------------------------------
%           ConsumptionBudget_MC:
%                               1) feces
%                               2) metabolism
%                               3) eggs (reproduction)
%                               4) predation
%                               5) senescence
%                               6) ba (biomass accumulation)
%                               7) em (emigration); NOTE: negative for immigration

%% step 10a: Adjust terminal benthic detritus in ConsumptionBudget
%           Removal (sequestration) of terminal benthic detritus is accounted for via emigration (em)
%           Add senescence term to em & set senescence term to 0
%           NOTE: ConsumptionBudget_MC dimensions (3D matrix: 7 X num_grps X num_MC)
ConsumptionBudget_MC(7, looky_terminalBNTHdetritus, :)	= ConsumptionBudget_MC(7, looky_terminalBNTHdetritus, :) + ConsumptionBudget_MC(5, looky_terminalBNTHdetritus, :);
ConsumptionBudget_MC(5, looky_terminalBNTHdetritus, :)	= 0;
% -------------------------------------------------------------------------


% % step 10b: apply and/or test alternate ConsumptionBudget terms------------
% ConsumptionBudget_MC(4, InertTracer2_RC, :)	= 0;	% CB_predation QQQ
% ConsumptionBudget_MC(5, [InertTracer2_RC InertTracer3_RC], :)	= 0;    % CB_senescence QQQ
% ConsumptionBudget_MC(1, InertTracer2_RC, :)	= 0.1;    % CB_feces
% ConsumptionBudget_MC(2, InertTracer2_RC, :)	= 0.1;    % CB_metabolism
% ConsumptionBudget_MC(3, InertTracer2_RC, :)	= 0;    % CB_eggs
% ConsumptionBudget_MC(4, InertTracer2_RC, :)	= 0.1;	% CB_predation
% ConsumptionBudget_MC(5, InertTracer2_RC, :)	= 0.1;    % CB_senescence
% ConsumptionBudget_MC(6, InertTracer2_RC, :)	= 0;	% CB_ba
% ConsumptionBudget_MC(7, InertTracer2_RC, :)	= -0.4;    % CB_em
% % -------------------------------------------------------------------------


% % step 10c: set senescence term(s) for primary producers to NPZD model values
% %           NOTE: activate this step if using Michaelis-Menten definition of primary producer senescence
% %           NOTE: this coding adpated especially for models which also have a macroalgae primary producer
% %           NOTE: in NPZD: non-grazing mortality rate (mmoles N/m3/d) = biomass (mmoles N/m3) * eta (1/d)
% %                 in ECOTRAN: senescence rate (mmoles N/m3/d)         = ConsumptionBudget_senescence (unitless) * ProductionRates (mmoles N/m3/d)
% %                 Substitution into ConsumptionBudget_senescence is therefore = eta (1/d) * 1/pb (d)
% %           NOTE: there are a lot of line of code here, most of this is just to
% %                 renormalize each ConsumptionBudget after replacing senescence with eta
% pb_PrimaryProducer                          = pb(looky_ANYPrimaryProducer)';	% (1/y); (horizontal vector: 1 X num_ANYPrimaryProd); NOTE transpose
% pb_PrimaryProducer                          = pb_PrimaryProducer / 365;         % (1/d); (horizontal vector: 1 X num_ANYPrimaryProd)
% pb_PrimaryProducer                          = repmat(pb_PrimaryProducer, [1, 1, num_MC]);                       % (1/d); (3D matrix: 1 X num_ANYPrimaryProd X num_MC); FFF replace this line if physicological parameters have MonteCarlo options
% eta_PrimaryProducer                         = reshape(MichaelisMenten_eta, [1, num_ANYPrimaryProd, num_MC]);	% reshape MichaelisMenten_eta; eta = non-grazing mortality; (1/d); (3D matrix: 1 X num_ANYPrimaryProd X num_MC)
% 
% % substitute eta for senescence and renormailze ConsumptionBudget_MC for each group to 1
% senescence_PrimaryProducer                  = eta_PrimaryProducer .* (1./pb_PrimaryProducer);	% (fraction); (3D matrix: 1 X num_ANYrimaryProducer X num_MC)
% NONsenescence_PrimaryProducer               = 1 - senescence_PrimaryProducer;                   % (fraction); (3D matrix: 1 X num_ANYrimaryProducer X num_MC)
% NONsenescence_PrimaryProducer               = repmat(NONsenescence_PrimaryProducer, [7, 1, 1]); % (3D matrix: 7 X num_ANYrimaryProducer X num_MC)
% 
% ConsumptionBudget_PrimaryProducer           = ConsumptionBudget_MC(:, looky_ANYPrimaryProducer, :);	% (fraction); (3D matrix: 7 X num_ANYrimaryProducer X num_MC)
% ConsumptionBudget_PrimaryProducer(5, :, :)	= 0;                                                    % temporarily set senescence to 0
% sum_ConsumptionBudget_PrimaryProducer       = sum(ConsumptionBudget_PrimaryProducer, 1);            % (fraction); (3D matrix: 1 X num_ANYrimaryProducer X num_MC)
% sum_ConsumptionBudget_PrimaryProducer       = repmat(sum_ConsumptionBudget_PrimaryProducer, [7, 1, 1]); % (3D matrix: 7 X num_ANYrimaryProducer X num_MC)
% ConsumptionBudget_PrimaryProducer           = ConsumptionBudget_PrimaryProducer ./ sum_ConsumptionBudget_PrimaryProducer; % (3D matrix: 7 X num_ANYrimaryProducer X num_MC)
% 
% ConsumptionBudget_PrimaryProducer(isnan(ConsumptionBudget_PrimaryProducer)) = 0;	% filter div/0 NaNs
% ConsumptionBudget_PrimaryProducer           = ConsumptionBudget_PrimaryProducer .* NONsenescence_PrimaryProducer; % (3D matrix: 7 X num_ANYrimaryProducer X num_MC)
% ConsumptionBudget_PrimaryProducer(5, :, :)	= senescence_PrimaryProducer;           % (fraction); (3D matrix: 7 X num_ANYrimaryProducer X num_MC)
% 
% sum_ConsumptionBudget_PrimaryProducer       = sum(ConsumptionBudget_PrimaryProducer, 1);	% (fraction); (3D matrix: 1 X num_ANYrimaryProducer X num_MC)
% [looky_grp, looky_MC]                       = find(squeeze(sum_ConsumptionBudget_PrimaryProducer) == 0); % find primary producer groups and Monte Carlo cases where there are no defined ConsumptionBudget fates
% ConsumptionBudget_PrimaryProducer(4, looky_grp, looky_MC) = 0.8;	% in cases where there are no defined ConsumptionBudget fates assume 80% goes to grazing and 20% goes to senescence
% ConsumptionBudget_PrimaryProducer(5, looky_grp, looky_MC) = 0.2;    % in cases where there are no defined ConsumptionBudget fates assume 80% goes to grazing and 20% goes to senescence
% ConsumptionBudget_MC(:, looky_ANYPrimaryProducer, :) = ConsumptionBudget_PrimaryProducer;	% substitute new ConsumptionBudget_PrimaryProducer back into ConsumptionBudget_MC; (fraction); (3D matrix: 7 X num_grps X num_MC)
% 
% % give warning that predation & senescence terms in ConsumptionBudget_MC have been redefined
% display('-->WARNING in ECOTRAN_DynamicScenario: primary producer predation & senescence terms in ConsumptionBudget are redefined with Michaelis-Menten eta values.')
% 
% % clear temporary terms that are not used again
% clear eta_PrimaryProducer pb_PrimaryProducer
% clear ConsumptionBudget_PrimaryProducer sum_ConsumptionBudget_PrimaryProducer
% clear looky_grp looky_MC
% % *************************************************************************





% *************************************************************************
%% STEP 11: define individual sub-regional food webs------------------------
%         NOTE: this can be replaced by the deliberate definition of
%               individual sub-regional food webs

% step 11a: define specific box types -------------------------------------
%          ConsumptionBudget_BoxType
%                               1) feces
%                               2) metabolism
%                               3) eggs (reproduction)
%                               4) predation
%                               5) senescence
%                               6) ba (biomass accumulation)
%                               7) em (emigration); NOTE: negative for immigration

switch switch_PhysicalModel
    case '2D_physics'
        % use this for vertically-resolved 4-box oceanic models
        EnergyBudget_BoxType         	= ones(num_grps, num_grps, num_boxes);	% (3D matrix: num_grps X num_grps X num_boxes)
        ConsumptionBudget_BoxType       = ones(7, num_grps, num_boxes);         % (3D matrix: 7 X num_grps X num_boxes)
    % end case (2D_physics) ---------
        
    case '3D_ROMS'
        % use this for vertically-resolved models and 3D models
        EnergyBudget_BoxType        = ones(num_grps, num_grps, num_boxes);	% (3D matrix: num_grps X num_grps X num_boxes)
        ConsumptionBudget_BoxType	= ones(7, num_grps, num_boxes);         % (3D matrix: 7 X num_grps X num_boxes)
 	% end case (3D_ROMS) ---------

end % switch (switch_PhysicalModel) ---------
% -------------------------------------------------------------------------


%% step 11b: pack BoxType definitions into ODEinput ------------------------
ODEinput.EnergyBudget_BoxType         	= EnergyBudget_BoxType;         % distinguish between surface & sub-surface EnergyBudget; (3D matrix: num_grps X num_grps X num_boxes)
ODEinput.ConsumptionBudget_BoxType      = ConsumptionBudget_BoxType;	% distinguish between surface & sub-surface ConsumptionBudget; (3D matrix: 7 X num_grps X num_boxes)
% *************************************************************************





% *************************************************************************
%% STEP 12: run ODE for each MonteCarlo model------------------------------

% step 12a: loop through each MonteCarlo model or treatment ---------------

% % QQQ turn off MonteCarlo_loop to run only the "type" model or for debugging 
MonteCarlo_loop     = 1; % MonteCarlo_loop is off for debugging
saveFile            = [SaveFile_directory SaveFile_label '_' date '.mat'];

% for MonteCarlo_loop = 1:num_MC
%     display(['MonteCarlo run ' num2str(MonteCarlo_loop) ' of ' num2str(num_MC)])
%     saveFile                        = [SaveFile_directory SaveFile_label '_' num2str(MonteCarlo_loop) '.mat'];
    
    % pick the current MonteCarlo food web --------------------------------
    current_biomass                	= biomass(:, 1);	% (t WWT/km2); NOTE: these are INITIAL biomass conditions; (vertical vector: num_grps X 1); NOTE: nutrients are zeros
    
    current_pb                      = pb(:, 1)';      	% (1/y); (horizontal vector: 1 X num_grps); NOTE transpose
    current_pb                      = current_pb / 365;	% (1/d); specific growth rate per day; (horizontal vector: 1 X num_grps)
    current_pb([looky_nutrients; looky_eggs; looky_ANYdetritus])	= 1;	% pb for nutrients & detritus are always 1 regardless of time-frame; (horizontal vector: 1 X num_grps)
    
    current_qb                      = qb(:, 1)';      	% (1/y); (horizontal vector: 1 X num_grps); NOTE transpose
    current_qb                      = current_qb / 365;	% (1/d); specific consumption rate per day; (horizontal vector: 1 X num_grps)
    current_qb([looky_nutrients; looky_eggs; looky_ANYdetritus; looky_fleets])	= 1;	% qb for nutrients, detritus, & fleets are always 1 regardless of time-frame; (horizontal vector: 1 X num_grps)
    
    current_EnergyBudget            = EnergyBudget_MC(:, :, MonteCarlo_loop);     	% (2D matrix: num_grps X num_grps)
    current_ConsumptionBudget       = ConsumptionBudget_MC(:, :, MonteCarlo_loop);	% (2D matrix: 7 X num_grps)
    
    current_fate_feces              = fate_feces(:, :, 1);      % (2D matrix: num_ANYdetritus X num_grps); FFF: fates are constant across Monte Carlo models as of now QQQ NO! I have new fates in ECOTRAN_MC, use them now???
	current_fate_metabolism       	= fate_metabolism(:, :, 1);	% (2D matrix: num_nutrients X num_grps); FFF: fates are constant across Monte Carlo models as of now    
    current_fate_eggs              	= fate_eggs(:, :, 1);    	% (2D matrix: num_eggs X num_grps); FFF: fates are constant across Monte Carlo models as of now
    % NOTE fate_predation is calculated below
    current_fate_senescence         = fate_senescence(:, :, 1);	% (2D matrix: num_ANYdetritus X num_grps); FFF: fates are constant across Monte Carlo models as of now
    
    current_TransferEfficiency     	= TransferEfficiency(MonteCarlo_loop, 1:num_grps);	% (horizontal vector: 1 X num_grps)
    
%     current_FunctionalResponseParams      = FunctionalResponseParams(:, :, MonteCarlo_loop); % QQQ DEACTIVATE IF NOT USING MONTE CARO FOR FUNCTIONAL RESPONSE; (2D matrix: CONSUMERS (num_grps) X prey group (num_grps)) replicated across clms (= producers)
	current_FunctionalResponseParams      = FunctionalResponseParams(:, :, 1); % QQQ USE IF NOT USING MONTE CARLO FOR FUNCTIONAL RESPONSE; (2D matrix: CONSUMERS (num_grps) X prey group (num_grps)) replicated across clms (= producers)

    current_MichaelisMenten_Vmax	= MichaelisMenten_Vmax(:,  1, MonteCarlo_loop);	% Vmax = maximum nutrient uptake rate; (1/d); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_KNO3	= MichaelisMenten_KNO3(:,  1, MonteCarlo_loop);	% KNO3 = NO3 half-saturation constant; (mmol N/m3); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_KNH4	= MichaelisMenten_KNH4(:,  1, MonteCarlo_loop);	% KNH4 = NH4 half-saturation constant; (mmol N/m3);  (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_alpha	= MichaelisMenten_alpha(:, 1, MonteCarlo_loop);	% alpha = initial slope of light response curve; (m2/W/d); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_psi     = MichaelisMenten_psi(:,   1, MonteCarlo_loop);	% psi = NO3 uptake inhibition by NH4; (m3/mmole N); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_w       = MichaelisMenten_w(:,     1, MonteCarlo_loop);	% w	= sinking rate; (m/d); (vertical vector: num_ANYPrimaryProd X 1)
	current_MichaelisMenten_eta     = MichaelisMenten_eta(:,   1, MonteCarlo_loop);	% eta = non-grazing mortality; (1/d); (vertical vector: num_ANYPrimaryProd X 1)
    % ---------------------------------------------------------------------

    
    % step 12b: build-up spatial sub-models -------------------------------
    %           NOTE: at this step, variable layers define spatial boxes and no longer define Monte Carlo models
    switch switch_SubModel

        case 'identical_SubModel'  % OPTION 1: use the same food web across the shelf ----------------
    
            disp('NOTE: IDENTICAL regional biological models')
            
            current_biomass                 = repmat(current_biomass, [1, num_boxes]);                          % (t WWT/km2); (2D matrix: num_grps X num_boxes)
            current_biomass                 = current_biomass .* repmat(spatial_BiomassScalers, [num_grps 1]);	% (t WWT/km2); (2D matrix: num_grps X num_boxes)

            current_pb                      = repmat(current_pb,                    [1, 1, num_boxes]);	% specific growth rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)
            current_qb                      = repmat(current_qb,                    [1, 1, num_boxes]);	% specific consumption rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)

            current_EnergyBudget            = repmat(current_EnergyBudget,          [1, 1, num_boxes]); % use for a single spatial definition of the EnergyBudget
            current_EnergyBudget            = current_EnergyBudget .* EnergyBudget_BoxType;             % adjust sub-surface EnergyBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: num_grps X num_grps X num_boxes)

            current_ConsumptionBudget       = repmat(current_ConsumptionBudget,     [1, 1, num_boxes]); % use for a single spatial definition of the ConsumptionBudget
            current_ConsumptionBudget       = current_ConsumptionBudget .* ConsumptionBudget_BoxType;	% adjust sub-surface ConsumptionBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: 7 X num_grps X num_boxes)

            current_fate_feces              = repmat(current_fate_feces,            [1, 1, num_boxes]); % (3D matrix: num_ANYdetritus X num_grps X num_boxes)
            current_fate_metabolism         = repmat(current_fate_metabolism,       [1, 1, num_boxes]); % (3D matrix: num_nutrients X num_grps X num_boxes)
            current_fate_eggs               = repmat(current_fate_eggs,             [1, 1, num_boxes]); % (3D matrix: num_eggs X num_grps X num_boxes)
            % NOTE fate_predation is calculated below
            current_fate_senescence         = repmat(current_fate_senescence,       [1, 1, num_boxes]); % (3D matrix: num_ANYdetritus X num_grps X num_boxes)

            current_TransferEfficiency      = repmat(current_TransferEfficiency,    [1, 1, num_boxes]); % (3D matrix: 1 X num_grps X num_boxes)        
            current_FunctionalResponseParams      = repmat(current_FunctionalResponseParams,    [1, 1, num_boxes]);	% (3D matrix: CONSUMERS (num_grps) X prey group (num_grps) X num_boxes) replicated across clms (= producers)

            current_MichaelisMenten_Vmax	= repmat(current_MichaelisMenten_Vmax,  [1, 1, num_boxes]);	% Vmax = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_KNO3	= repmat(current_MichaelisMenten_KNO3,  [1, 1, num_boxes]);	% KNO3 = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_KNH4	= repmat(current_MichaelisMenten_KNH4,  [1, 1, num_boxes]);	% KNH4 = NH4 half-saturation constant; (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_alpha	= repmat(current_MichaelisMenten_alpha, [1, 1, num_boxes]);	% alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_psi     = repmat(current_MichaelisMenten_psi,   [1, 1, num_boxes]);	% psi = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_w       = repmat(current_MichaelisMenten_w,     [1, 1, num_boxes]); % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_eta     = repmat(current_MichaelisMenten_eta,   [1, 1, num_boxes]);	% eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
        % end case (identical_SubModel) ------------
            
        case 'independent_SubModel'  % OPTION 2: use independently defined webs for each shelf zone ----
                                     %           use for multiple definitions of the EnergyBudget
                                     % FFF future code could easily allow for Monte Carlo models for each sub-region

            disp('NOTE: LOADING independently defined regional biological models')
                                      
            % load invidividual regional models
            % SSS use for CNP----
            epipelagic_CNP_06122022_B
            mesopelagic_CNP_06122022_B
            bathypelagic_CNP_06122022_B
            benthic_CNP_06122022_B
            % SSS ----
             
            current_EnergyBudget                = repmat(current_EnergyBudget, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_grps X num_grps X num_boxes)
            current_EnergyBudget(:, :, 1)       = EPIPELAGIC.EnergyBudget;
            current_EnergyBudget(:, :, 2)       = MESOPELAGIC.EnergyBudget;
            current_EnergyBudget(:, :, 3)       = BATHYPELAGIC.EnergyBudget;
            current_EnergyBudget(:, :, 4)       = BENTHIC.EnergyBudget;
            current_EnergyBudget                = current_EnergyBudget .* EnergyBudget_BoxType; % adjust sub-surface EnergyBudget matrices; allow only certain processes in sub-surface boxes; (fraction); (3D matrix: num_grps X num_grps X num_boxes)

            current_ConsumptionBudget           = repmat(current_ConsumptionBudget, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: 7 X num_grps X num_boxes)
            current_ConsumptionBudget(:, :, 1) 	= EPIPELAGIC.ConsumptionBudget;
            current_ConsumptionBudget(:, :, 2) 	= MESOPELAGIC.ConsumptionBudget;
            current_ConsumptionBudget(:, :, 3) 	= BATHYPELAGIC.ConsumptionBudget;
            current_ConsumptionBudget(:, :, 4)	= BENTHIC.ConsumptionBudget;
            current_ConsumptionBudget         	= current_ConsumptionBudget .* ConsumptionBudget_BoxType; % adjust sub-surface ConsumptionBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: 7 X num_grps X num_boxes)

            current_biomass                     = repmat(current_biomass, [1, num_boxes]);	% replicate box 1 model to initialize; (mmoles N/m3); NOTE: these are INITIAL biomass conditions; (2D matrix: num_grps X num_boxes); NOTE: nutrients are zeros; NOTE: only the definition for primary producer biomass is used, it is used to drive initial model production values
            current_biomass(:, 1)               = EPIPELAGIC.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;
            current_biomass(:, 2)               = MESOPELAGIC.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;
            current_biomass(:, 3)               = BATHYPELAGIC.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;
            current_biomass(:, 4)               = BENTHIC.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;

            current_pb                          = repmat(current_pb, [1, 1, num_boxes]); % replicate box 1 model to initialize; specific growth rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)
            current_pb(:, :, 1)                 = (EPIPELAGIC.pb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose
            current_pb(:, :, 2)                 = (MESOPELAGIC.pb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose
            current_pb(:, :, 3)                 = (BATHYPELAGIC.pb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose
            current_pb(:, :, 4)                 = (BENTHIC.pb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose

            current_qb                          = repmat(current_qb, [1, 1, num_boxes]); % replicate box 1 model to initialize; specific consumption rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)
            current_qb(:, :, 1)                 = (EPIPELAGIC.qb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose
            current_qb(:, :, 2)                 = (MESOPELAGIC.qb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose
            current_qb(:, :, 3)                 = (BATHYPELAGIC.qb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose
            current_qb(:, :, 4)                 = (BENTHIC.qb)'; % NOTE: conversion to 1/d was already done manually; NOTE transpose

            current_fate_feces                  = repmat(current_fate_feces, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_ANYdetritus X num_grps X num_boxes)
            current_fate_feces(:, :, 1)         = EPIPELAGIC.fate_feces;
            current_fate_feces(:, :, 2)         = MESOPELAGIC.fate_feces;
            current_fate_feces(:, :, 3)         = BATHYPELAGIC.fate_feces;
            current_fate_feces(:, :, 4)         = BENTHIC.fate_feces;

            current_fate_metabolism             = repmat(current_fate_metabolism, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_nutrients X num_grps X num_boxes)
            current_fate_metabolism(:, :, 1)	= EPIPELAGIC.fate_metabolism;
            current_fate_metabolism(:, :, 2)	= MESOPELAGIC.fate_metabolism;
            current_fate_metabolism(:, :, 3)	= BATHYPELAGIC.fate_metabolism;
            current_fate_metabolism(:, :, 4)	= BENTHIC.fate_metabolism;

            current_fate_eggs                   = repmat(current_fate_eggs, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_eggs X num_grps X num_boxes)
            current_fate_eggs(:, :, 1)          = EPIPELAGIC.fate_eggs;
            current_fate_eggs(:, :, 2)          = MESOPELAGIC.fate_eggs;
            current_fate_eggs(:, :, 3)          = BATHYPELAGIC.fate_eggs;
            current_fate_eggs(:, :, 4)          = BENTHIC.fate_eggs;

            % NOTE: fate_predation is calculated below

            current_fate_senescence             = repmat(current_fate_senescence, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_ANYdetritus X num_grps X num_boxes)
            current_fate_senescence(:, :, 1)	= EPIPELAGIC.fate_senescence;
            current_fate_senescence(:, :, 2)	= MESOPELAGIC.fate_senescence;
            current_fate_senescence(:, :, 3)	= BATHYPELAGIC.fate_senescence;
            current_fate_senescence(:, :, 4)	= BENTHIC.fate_senescence;

            % these terms assumed to be the same in all sub-regions
            current_TransferEfficiency      = repmat(current_TransferEfficiency, [1, 1, num_boxes]);  	% (3D matrix: 1 X num_grps X num_boxes)        
            current_FunctionalResponseParams      = repmat(current_FunctionalResponseParams, [1, 1, num_boxes]);	% (3D matrix: CONSUMERS X prey group X num_boxes) replicated across clms (= producers)

            current_MichaelisMenten_Vmax	= repmat(current_MichaelisMenten_Vmax, [1, 1, num_boxes]);	% Vmax = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_KNO3	= repmat(current_MichaelisMenten_KNO3, [1, 1, num_boxes]);	% KNO3 = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_KNH4	= repmat(current_MichaelisMenten_KNH4, [1, 1, num_boxes]);	% KNH4 = NH4 half-saturation constant; (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_alpha	= repmat(current_MichaelisMenten_alpha, [1, 1, num_boxes]);	% alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_psi     = repmat(current_MichaelisMenten_psi, [1, 1, num_boxes]);	% psi = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_w       = repmat(current_MichaelisMenten_w, [1, 1, num_boxes]);     % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_eta     = repmat(current_MichaelisMenten_eta, [1, 1, num_boxes]);	% eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            
            
% % %             % QQQ QQQ QQQ QQQ QQQ QQQ QQQ
% % %             % QQQ modifications for nutrient input calibrations
% % %             % 1) Set ConsumptionBudget predation term to 0.5 to 1.0 in increments of 0.1
% % %             CB_predation                = 0.99; % SSS predation size in ConsumptionBudget (0.8 is value of type models); NOTE: using calibrated GoMexOcn value = 0.99
% % % 
% % %             current_ConsumptionBudget(4, [rc_detritusA rc_detritusB], :)	= CB_predation;
% % %             current_ConsumptionBudget(5, [rc_detritusA rc_detritusB], :)	= 1 - CB_predation;
% % % 
% % %             % 2) Set EnergyBudget flow from detritus to pelagicBacteria to 0 to 1.0 in increments of 0.1
% % %             detritus_inflow     = 0; % SSS set EnergyBudget flow from detritus to pelagicBacteria
% % %             current_EnergyBudget([rc_PelagicBacteria rc_detritusA rc_detritusB rc_detritusC], [rc_detritusA rc_detritusB], :) = 0; % set flow from detritus to pelagicBacteria or to other detritus to 0
% % %             EB_detritusClms     = current_EnergyBudget(:, [rc_detritusA rc_detritusB], :); % copy detritus clms from E.B. in order to modify them
% % %             sum_EB              = sum(EB_detritusClms, 1); % un-modified sum of E.B. excluding flows to detritus and bacteria (= sum of flow to other consumers)
% % %             EB_detritusClms     = EB_detritusClms ./ repmat(sum_EB, [num_grps, 1, 1]); % normalize EB_detritusClms to 1
% % %             EB_detritusClms(isnan(EB_detritusClms)) = 0; % correct for any div/0 error
% % % 
% % %             EB_predation        = max(0, (CB_predation - detritus_inflow)); % allowed total predation by non-bacteria
% % %             EB_detritusClms     = EB_detritusClms .* EB_predation; % scale EB_detritusClms to allowed total predation by non-bacteria
% % % 
% % %             detritus_inflow     = min(CB_predation, detritus_inflow); % cap detritus_inflow at value defined for CB_predation
% % %             EB_detritusClms(rc_PelagicBacteria, :, :) = detritus_inflow;
% % %             EB_senescence       = fate_senescence(:, [rc_detritusA rc_detritusB], :) * (1 - CB_predation);
% % %             EB_detritusClms(looky_ANYdetritus, :, :) = repmat(EB_senescence, [1, 1, num_boxes]);
% % %             current_EnergyBudget(:, [rc_detritusA rc_detritusB], :) = EB_detritusClms;
% % % 
% % %             % 3) Redefined fate_predation is done below
% % % 
% % %             % 4) clear temporary variables
% % %             clear EB_detritusClms EB_senescence EB_predation sum_EB detritus_inflow
% % % 
% % %             % QQQ QQQ QQQ QQQ QQQ QQQ QQQ
            
        % end case ('independent_SubModel') ------------
            
    end % switch switch_SubModel
	% ---------------------------------------------------------------------

                
	% step 12c: calculate fate_predation ----------------------------------
    %           NOTE: multiple MC versions must replace the single "type" fate_predation from the main ECOTRAN code; 
    sum_predation                                           = sum(current_EnergyBudget(looky_livingANDfleets, :, :)); % (3D matrix: 1 X num_grps X num_boxes)
    current_fate_predation                                  = current_EnergyBudget(looky_livingANDfleets, :, :) ./ repmat(sum_predation, [num_livingANDfleets, 1, 1]); % (3D matrix: num_livingANDfleets X num_grps X num_boxes)
    current_fate_predation(isnan(current_fate_predation))	= 0; % correct div/0 errors
    % ---------------------------------------------------------------------
    
    
    % step 12d: build time-series of varying physiologies -----------------
    %      (rows = time, clms = ECOTRAN groups, layers = spatial boxes)
    %      NOTE: FFF at this point in the code, we can read in time-series changes for each of these terms
    %            this would allow for seasonal reproduction differences (ConsumptionBudget_eggs), & 
    %            seasonal migration changes (ConsumptionBudget_em)
    %      NOTE: the need to break out separate rows of ConsumptionBudget
    %            is to prevent having to deal with 4D matrices:
    %                                                   row 1: feces
    %                                                   row 2: metabolism
    %                                                   row 3: eggs
    %                                                   row 4: predation
    %                                                   row 5: senescence
    %                                                   row 6: ba
    %                                                   row 7: em

	current_pb                      = repmat(current_pb, [num_t, 1, 1]);	% specific growth rate per day; (3D matrix: num_t X num_grps X num_boxes)
	current_qb                      = repmat(current_qb, [num_t, 1, 1]);	% specific consumption rate per day; (3D matrix: num_t X num_grps X num_boxes)

    %     % OPTION: calculate primary producer production rate (pb)
    %     % FFF need to update this option for ECOTRAN 2; move outside MC_loop?
    %     [current_pb_PrimaryProducers, biomass_PrimaryProducers]	= f_SeasonalPB_PrimaryProducers(ECOTRANphysics, ODEinput, EwEResult); % (1/d); (3D matrix: num_t X num_ANYPrimaryProd X num_boxes)
    %     current_pb(:, looky_ANYPrimaryProducer, :)              = current_pb_PrimaryProducers; % plug in new primary producer timeseries; QQQ -->??? use Box 1 as the full domain definition
    %     current_qb(:, looky_ANYPrimaryProducer, :)              = current_pb_PrimaryProducers; % plug in new primary producer timeseries; QQQ -->??? use Box 1 as the full domain definition
    
    ConsumptionBudget_feces      	= repmat(current_ConsumptionBudget(1, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_metabolism	= repmat(current_ConsumptionBudget(2, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_eggs        	= repmat(current_ConsumptionBudget(3, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_predation     = repmat(current_ConsumptionBudget(4, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_senescence	= repmat(current_ConsumptionBudget(5, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_ba            = repmat(current_ConsumptionBudget(6, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_em            = repmat(current_ConsumptionBudget(7, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    
    
    % step 12e: rescale ConsumptionBudget_metabolism according to Q10 temperature responsee
    switch switch_physiology
        case 'physiology_ON'
            
            disp('NOTE: ConsumptionBudget_metabolism response to temperature change is ON')

            ConsumptionBudget_metabolism_new	= ConsumptionBudget_metabolism .* Q10_scaler; % (3D matrix: num_t X num_grps X num_boxes)
            ConsumptionBudget_metabolism        = ConsumptionBudget_metabolism_new; % (3D matrix: num_t X num_grps X num_boxes)

        % end case 'physiology_ON' ----------------------------------------
            
        case 'physiology_OFF'
            
            disp('NOTE: ConsumptionBudget_metabolism response to temperature change is OFF')
            
        % end case 'physiology_OFF' ---------------------------------------
        
    end % (switch switch_physiology) --------------------------------------
    
    
    % step 12f: rescale ConsumptionBudget_senescence according to Q10 temperature response AND Mining Plume response
    switch switch_MiningPlume_physiology
        case 'MiningPlume_ON'
            
            disp('NOTE: ConsumptionBudget_senescence response to Mining Plume is ON')

            ConsumptionBudget_senescence_new	= ConsumptionBudget_senescence .* mortality_PlumeEffect_scaler; % (3D matrix: num_t X num_grps X num_boxes)
            ConsumptionBudget_senescence        = ConsumptionBudget_senescence_new; % (3D matrix: num_t X num_grps X num_boxes)

        % end case 'MiningPlume_ON' ---------------------------------------
            
        case 'MiningPlume_OFF'
            
            disp('NOTE: ConsumptionBudget_senescence response to Mining Plume is OFF')
            
        % end case 'MiningPlume_OFF' --------------------------------------
        
    end % (switch switch_MiningPlume_physiology) --------------------------
    
    % ---------------------------------------------------------------------
    
    
    % % step 12f: OPTIONAL: make scenario changes to ConsumptionBudget ----
    % %           NOTE: just a test scenario
    % ConsumptionBudget_ba(10000:20000, rc_fleet_midWater_trawls, 1)              = 0.5;    % QQQ NCC: high mortality for 1 month
    % ConsumptionBudget_senescence(10000:20000, rc_fleet_midWater_trawls, 1)      = 0.1307; % QQQ NCC: high mortality for 1 month
    % ConsumptionBudget_predation(950:1000, rc_anchovy, 1)                        = 0.1136; % QQQ NCC: high mortality for 1 month
    % 
    % ConsumptionBudget_ba(10000:20000, rc_fleet_midWater_trawls, 4)              = 0.9;    % QQQ NCC: high mortality for 1 month
    % ConsumptionBudget_ba(950:1000, str2num(SPECIES), 1)                         = 0.5;    % QQQ NCC: high mortality for 1 month
    % % -------------------------------------------------------------------


    %% step 12g: pack variables needed for ODE ----------------------------
    ODEinput.pb                             = current_pb;                       % (1/d); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.qb                             = current_qb;                       % (1/d); (3D matrix: num_t X num_grps X num_boxes)
    
    ODEinput.EnergyBudget                   = current_EnergyBudget;             % (proportions); (3D matrix: num_grps X num_grps X num_boxes)
    
    ODEinput.ConsumptionBudget_feces        = ConsumptionBudget_feces;          % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_metabolism	= ConsumptionBudget_metabolism;     % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_eggs         = ConsumptionBudget_eggs;           % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_predation	= ConsumptionBudget_predation;      % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_senescence	= ConsumptionBudget_senescence;     % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_ba           = ConsumptionBudget_ba;             % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_em           = ConsumptionBudget_em;             % (proportions); (3D matrix: num_t X num_grps X num_boxes)

    ODEinput.fate_feces                     = current_fate_feces;               % (proportions); (3D matrix: num_ANYdetritus X num_grps X num_boxes)
    ODEinput.fate_metabolism                = current_fate_metabolism;          % (proportions); (3D matrix: num_nutrients X num_grps X num_boxes)
    ODEinput.fate_eggs                  	= current_fate_eggs;                % (proportions); (3D matrix: num_eggs X num_grps X num_boxes)
    ODEinput.fate_predation                 = current_fate_predation;           % (proportions); (3D matrix: num_livingANDfleets X num_grps X num_boxes)
    ODEinput.fate_senescence                = current_fate_senescence;          % (proportions); (3D matrix: num_ANYdetritus X num_grps X num_boxes)
    
    ODEinput.TransferEfficiency           	= current_TransferEfficiency;       % (3D matrix: 1 X num_grps X num_boxes)
    ODEinput.FunctionalResponseParams     	= current_FunctionalResponseParams;       % (vulnerability of producer, m_p); (3D matrix: CONSUMERS (num_grps) X prey group (num_grps) X num_boxes) replicated across clms (= producers)

    ODEinput.MichaelisMenten_Vmax           = current_MichaelisMenten_Vmax;     % Vmax = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_KNO3           = current_MichaelisMenten_KNO3;     % KNO3 = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_KNH4           = current_MichaelisMenten_KNH4;     % KNH4 = NH4 half-saturation constant; (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_alpha          = current_MichaelisMenten_alpha;	% alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_psi            = current_MichaelisMenten_psi;      % psi = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_w              = current_MichaelisMenten_w;        % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_eta            = current_MichaelisMenten_eta;      % eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    
    ODEinput.q_TemperatureScaler            = q_TemperatureScaler .* q_PlumeEffect_scaler; % Thornton-Lessem temperature and MiningPlume effects are MULTIPLIED for adjustment to consumption rate (ASSUMES THESE 2 RESPONSES ARE ADDITIVE); value between 0 and 1; (3D matrix: num_t X num_grps X num_boxes)
    % *********************************************************************
    
    
    
    
    
	% *********************************************************************
	%% STEP 13: calculate INITIAL production rate conditions --------------
    %           NOTE: (INITIAL production = consumption inflow)
	%           NOTE: several options are provided
    %           QQQ after debugging, use a consistent summer season base time (consolodate t_initial & t_mean calls) 
    
    t_initial                               = 1;
    
	% step 13a: run 1 of several options for defining initial condition rates
    switch switch_INITIALproduction

        case 'INITIALproduction_MichaelisMenton' % METHOD 1: use for driving with primary production defined by Michaelis-Menton uptake
                                                  % QQQ this section needs proofing in ECOTRAN 2
        
            biomass_PrimaryProducer_t                           = current_biomass(looky_ANYPrimaryProducer, :);         % biomass of primary producers in current MonteCarlo model; (t WWT/km2); (2D matrix: num_ANYPrimaryProd X num_boxes)
            biomass_macroalgae_t                                = current_biomass(looky_macroalgae, :);                 % biomass of Macroalgae in current MonteCarlo model; (t WWT/km2); (2D matrix: num_macroalgae X num_boxes)
            pb_PrimaryProducer_t                                = current_pb(t_initial, looky_ANYPrimaryProducer, :);	% pb of primary producers in current MonteCarlo model; (1/d); (3D matrix: 1 X num_ANYPrimaryProd X num_boxes)
            pb_macroalgae_t                                     = current_pb(t_initial, looky_macroalgae, :);           % pb of Macroalgae in current MonteCarlo model; (1/d); (3D matrix: 1 X num_macroalgae X num_boxes)

            % special case to accomodate any macroalgae
            production_PrimaryProducer_WWT                      = reshape(biomass_PrimaryProducer_t, [1, num_ANYPrimaryProd, num_boxes]) .* pb_PrimaryProducer_t;       % (t WWT/km2/d); (3D matrix: 1 X num_ANYPrimaryProd X num_boxes)
            production_macroalgae_WWT                           = reshape(biomass_macroalgae_t, [1, num_macroalgae, num_boxes])          .* pb_macroalgae_t;            % (t WWT/km2/y); (3D matrix: 1 X num_macroalgae X num_boxes)
            ProductionFraction_macroalgae                       = production_macroalgae_WWT   ./ repmat(sum(production_PrimaryProducer_WWT), [1, num_macroalgae, 1]);	% macroalgae fraction(s) of total primary production based on EwE p=b*pb; (3D matrix: 1 X num_macroalgae X num_boxes)
            ProductionFraction_macroalgae(isnan(ProductionFraction_macroalgae)) = 0;                                                         % fix div/0 NaNs
            ProductionFraction_macroalgae                       = reshape(ProductionFraction_macroalgae, [num_macroalgae, 1, num_boxes]);    % (3D matrix: num_macroalgae X 1 X num_boxes)

            % primary producer biomass conversion from (t WWT/km2) to (mmole N/m3)
            EuphoticDepth_t                                     = interp1(t_grid, EuphoticDepth, t_initial);            % euphotic zone depth @ t; (m); (horizontal vector: 1 X num_boxes)
            EuphoticDepth_t                                     = repmat(EuphoticDepth_t, [num_ANYPrimaryProd, 1]);     % replicate box heights for each primary producer; (m); (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/(1000*1000));          % area to volumetric conversion; (t WWT/m2); (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t .* (1./EuphoticDepth_t);    % area to volumetric conversion; (t WWT/m3); (2D matrix: num_ANYPrimaryProd X num_boxes); NOTE use of EuphoticDepth instead of MLD because we are starting out with phyto biomasses defined as vertically integrated over whole boxes
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * 1000000;                  % (g WWT/m3);   (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * 1000;                     % (mg WWT/m3);  (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/WWT_to_C);             % (mg C/m3);    (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/atomic_mass_C);        % (mmole C/m3); (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/C_to_N_phytoplankton); % (mmole N/m3); (2D matrix: num_ANYPrimaryProd X num_boxes);

            biomass_PrimaryProducer_t                           = reshape(biomass_PrimaryProducer_t, [num_ANYPrimaryProd, 1, num_boxes]); % (mmole N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            NO3_t                                               = repmat(NO3initial_rate, [num_ANYPrimaryProd, 1]);       	% (mmole NO3/m3); (2D matrix: num_ANYPrimaryProd X num_boxes)
            % NH4_t                                               = repmat(NH4initial_rate, [num_ANYPrimaryProd, 1]);        	% (mmole NH4/m3); (2D matrix: num_ANYPrimaryProd X num_boxes)
            NO3_t                                               = reshape(NO3_t, [num_ANYPrimaryProd, 1, num_boxes]);	% (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes)
            % NH4_t                                               = reshape(NH4_t, [num_ANYPrimaryProd, 1, num_boxes]);	% (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes)
            Io_t                                                = interp1(t_grid, Io, t_mean);                       	% surface PAR light intensity @ t; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (scaler); NOTE: using t_mean rather than t_initial
            MLD_t                                               = interp1(t_grid, MLD, t_mean);                      	% MLD @ t; (m); (scaler)
            NH4fraction_bnth                                    = zeros(num_ANYPrimaryProd, 1, num_boxes);            	% value 0 - 1; (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); (for Michaelis-Menten calcs); for initial conditions, any NH4 is asigned as pelagic

            [Q_cp_NO3, Q_cp_NH4, Q_cp_plgcNH4, Q_cp_bnthNH4, PB_MichaelisMenten] = f_MichaelisMenten_05152016(ODEinput, biomass_PrimaryProducer_t, NO3_t, NH4_t, NH4fraction_bnth, Io_t, MLD_t); % phytoplankton uptake rate of NO3 & NH4 (mmole N/m3/d) & pb @ t as calculated from MichaelisMenton uptake (1/y); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes) QQQ proof code for ECOTRAN II
            production_initial_PrimaryProducer                  = Q_cp_NO3 + Q_cp_NH4;                                                      % phytoplankton production rate; (mmole N/m3/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            production_initial_PrimaryProducer               	= squeeze(production_initial_PrimaryProducer);                                        % (mmole N/m3/d); (2D matrix: num_ANYPrimaryProd X num_boxes)
            production_initial_driver                        	= zeros(num_grps, num_boxes);                                          % initialize; (2D matrix: num_grps X num_boxes)
            production_initial_driver(looky_ANYPrimaryProducer, :)	= production_initial_PrimaryProducer;                                                 % (2D matrix: num_grps X num_boxes)
            if num_macroalgae > 0
                production_initial_driver(looky_macroalgae, :)      = repmat(sum(production_initial_PrimaryProducer), [num_macroalgae, 1]) .* ProductionFraction_macroalgae; % (2D matrix: num_grps X num_boxes)
            end
            production_initial_driver                           = reshape(production_initial_driver, [num_grps, 1, num_boxes]);          % (3D matrix: num_grps X 1 X num_boxes)
            production_initial_driver                           = reshape(production_initial_driver, [1, num_grps, num_boxes]);          % (3D matrix: 1 X num_grps X num_boxes) 
            
            % finalize initial condition calculations
            [production_initial, fname_InitialProductionRates]	= f_InitialProductionRates_11202019(ODEinput, production_initial_driver, t_initial);  % initial or mean production rates (actually consumption inflow); (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE: code does NOT make transfer of bnthNH4 from surface to sub-surface boxes
            
            % paste in initial NO3 & NH4 input rates
            production_initial(looky_NO3, :)      	= NO3initial_rate;   % append initial NO3 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
            % production_initial(looky_plgcNH4, :)    = NH4initial_rate;   % append initial NH4 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
        % end (case 'INITIALproduction_MichaelisMenton') ------------------
        
        case 'INITIALproduction_pb'	% METHOD 2: use for driving initial model conditions with primary production defined by p = [(p/b) * b]
                                        %           NOTE: actually q = [(q/b) * b] because p & q are same thing for primary producers
        
            biomass_PrimaryProducer_t                           = current_biomass(looky_ANYPrimaryProducer, :);      	% biomass of primary producers in current MonteCarlo model; (t WWT/km2); (2D matrix: num_ANYPrimaryProd X num_boxes)
            pb_PrimaryProducer_t                             	= current_pb(t_initial, looky_ANYPrimaryProducer, :);	% pb of primary producers in current MonteCarlo model;      (1/d);       (3D matrix: 1 X num_ANYPrimaryProd X num_boxes); NOTE: pb = qb for primary producers; FFF allow for pb changes throughout the year

            % primary producer biomass conversion from (t WWT/km2) to (mmole N/m3)
            EuphoticDepth_t                                     = interp1(t_grid, EuphoticDepth, t_initial);            % euphotic zone depth @ t; (m); (horizontal vector: 1 X num_boxes)
            EuphoticDepth_t                                     = repmat(EuphoticDepth_t, [num_ANYPrimaryProd, 1]);     % replicate box heights for each primary producer; (m); (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/(1000*1000));          % area to volumetric conversion; (t WWT/m2); (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t .* (1./EuphoticDepth_t);    % area to volumetric conversion; (t WWT/m3); (2D matrix: num_ANYPrimaryProd X num_boxes); NOTE use of EuphoticDepth instead of MLD because we are starting out with phyto biomasses defined as vertically integrated over whole boxes
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * 1000000;                  % (g WWT/m3);   (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * 1000;                     % (mg WWT/m3);  (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/WWT_to_C);             % (mg C/m3);    (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/atomic_mass_C);        % (mmole C/m3); (2D matrix: num_ANYPrimaryProd X num_boxes);
            biomass_PrimaryProducer_t                           = biomass_PrimaryProducer_t * (1/C_to_N_phytoplankton); % (mmole N/m3); (2D matrix: num_ANYPrimaryProd X num_boxes);

            % initial production (actually consumption inflow) driver vector for each Box; (mmole N/m3/d); (3D matrix: 1 X num_groups X num_boxes)
            production_initial_PrimaryProducer                        	= reshape(biomass_PrimaryProducer_t, [1, num_ANYPrimaryProd, num_boxes]) .* pb_PrimaryProducer_t;    % primary production rate; (mmole N/m3/d); (3D matrix: 1 X num_ANYPrimaryProd X num_boxes)
            production_initial_driver                                   = zeros(1, num_grps, num_boxes);        % initialize DriverProductionVector; (3D matrix: 1 X num_grps X num_boxes)
            production_initial_driver(1, looky_ANYPrimaryProducer, :)	= production_initial_PrimaryProducer;	% plug in primary production rates; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
            % end Method 2 ----

            % finalize initial condition calculations
            [production_initial, fname_InitialProductionRates] = f_InitialProductionRates_02012022(ODEinput, production_initial_driver, t_initial);  % initial or mean production rates (actually consumption inflow); (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE: code does NOT make transfer of bnthNH4 from surface to sub-surface boxes
            
            % paste in initial NO3 & NH4 input rates
            production_initial(looky_NO3, :)      	= NO3initial_rate;   % append initial NO3 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
            % production_initial(looky_plgcNH4, :)    = NH4initial_rate;   % append initial NH4 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
        % end (case 'INITIALproduction_pb') -------------------------------
            
        case 'INITIALproduction_SubModel'	% METHOD 3: use for driving initial model conditions with values loaded along with regional sub-model definitions 
                                            % For example, use for oceanic GoMexOcn, CNP models----
                                            % CBB: This is where we add in initial conditions from a long run from a different file (save end 5-10 yrs re_Y average consumption rate).
            
            disp('USING INITIAL CONDITIONS FROM EXTERNAL FILE: ')
                                            
            fname_InitialProductionRates	= 'pre-defined in independent sub-model';
            production_initial              = [EPIPELAGIC.production_initial MESOPELAGIC.production_initial BATHYPELAGIC.production_initial BENTHIC.production_initial]; % initial q; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes)
           
            % % paste in initial NO3 & NH4 input rates
            % production_initial(looky_NO3, :)      	= NO3initial_rate;   % append initial NO3 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
            % production_initial(looky_plgcNH4, :)    = NH4initial_rate;   % append initial NH4 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
        
        % end case (INITIALproduction_SubModel) ---------------------
    
            
        case 'INITIALproduction_nutrients'	% METHOD 4: use for driving initial model conditions with mean annual nutrient input rates
                                            %           NOTE: NH4 uptake is deactivated in f_InitialProductionRates_02012022, so production_initial will be lower for METHOD 4 than for METHOD 2 (the latter implicitly includes recylced primary production)
            
            disp('INITIAL CONDITIONS: calculated from nutrient input rate')
                                            
            production_initial_driver                           = zeros(1, num_grps, num_boxes);        % initialize DriverProductionVector; (3D matrix: 1 X num_grps X num_boxes)
            production_initial_driver(1, looky_NO3, :)          = NO3initial_rate;	% plug in NO3 input rate; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
            production_initial_driver(1, rc_plgcNH4, :)         = NH4initial_rate;	% plug in NH4 input rate; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)

            % finalize initial condition calculations
            [production_initial, fname_InitialProductionRates]	= f_InitialProductionRates_02012022(ODEinput, production_initial_driver, t_initial);  % QQQ NH4 uptake turned OFF; initial or mean production rates (actually consumption inflow); (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE: code does NOT make transfer of bnthNH4 from surface to sub-surface boxes

            % paste in initial NO3 & NH4 input rates
            production_initial(looky_NO3, :)      	= NO3initial_rate;   % append initial NO3 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
            production_initial(looky_plgcNH4, :)    = NH4initial_rate;   % append initial NH4 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
        % end (case 'INITIALproduction_nutrients') --------------------
            
    end % (switch switch_INITIALproduction) -------------------------------
    
    production_initial                  	= reshape(production_initial, [num_grps, 1, num_boxes]);	% reshape boxes as layers; (3D matrix: num_grps X 1 X num_boxes); (mmole N/m3/d)
    productionC_initial_repmat            	= repmat(production_initial,   [1, num_grps, 1]);           % replicate groups across columns; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: this is needed for the ODE
    % ---------------------------------------------------------------------
    
    
    % step 13b: pack more initial conditions for ODE ----------------------
    ODEinput.production_initial            	= production_initial;               % production rates to use as initial conditions; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes); NOTE: used in ODE only for special cases of step-thru debugging or for quadratic functional responses
    ODEinput.productionC_initial_repmat   	= productionC_initial_repmat;       % initial conditions reshaped; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: replicated vertical vectors across columns
    %% ********************************************************************



        
    
    %% ********************************************************************
    % STEP 14: solve the dynamic model-------------------------------------
    switch switch_ODEsolver
  
        case 'CppSolver' % C++ solver
            % NOTE: >> mex mex_ECOTRANode_11272020.cpp -I/usr/local/include/ % compile mex function
            
            disp('NOTE: using C++ ODE solver')

            % step 14a: prepare ECOTRAN variables for using C++ ODE solver mex function
            %           Pack parameters & drivers along proper dimensions.
            disp('prepare variables for C++ ODE solver...')
            
%             [AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct] = f_PrepMexODE_10272020(ODEinput); % withOUT Thornton-Lessem
%             [AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct] = f_PrepMexODE_05312022(ODEinput); % WITH Thornton-Lessem
            [AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct] = f_PrepMexODE_08212022(ODEinput); % pack variables to send to ODE solver; WITH Thornton-Lessem

    
%             fname_ECOTRANode                 = 'mex_ECOTRANode_03092021_D';	% withOUT Thornton-Lessem
            fname_ECOTRANode                 = 'mex_ECOTRANode_05312022'; % WITH Thornton-Lessem
            fname_PhysicalFlux_intraODE      = 'NA';	% SSS            
            % -------------------------------------------------------------
    
            % step 14b: run the model in C++ ------------------------------
            disp(['Running C++ solver: ' fname_ECOTRANode])
            tic            
            % [output_Cpp]            = mex_ECOTRANode_03092021_D(AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct); % withOUT Thornton-Lessem
            [output_Cpp]            = mex_ECOTRANode_05312022(AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct); % WITH Thornton-Lessem
            time_ODE                = toc
            store_T                 = output_Cpp.mat_obs_times; % (d); (vertical vector: num_t X 1)
            store_ProductionRates	= output_Cpp.mat_obs_states; % (mmole N/m3/d); (2D matrix: num_t X (num_grps*num_boxes))
            % -------------------------------------------------------------
    
            % step 13c: unstack result (store_ProductionRates) to retrieve spatial boxes
            re_Y    = reshape(store_ProductionRates, [num_t, num_grps, num_boxes]); % (mmole N/m3/d); (3D matrix: time X groups X num_boxes)
        % end (case 'CppSolver') ------------------------------------------

    
        case 'MatlabSolver' % Matlab solver
            % NOTE: calculation uses ode23t (trial-and-error suggests a bit better performance than ODE45)
            % NOTE: size of ProductionRates is 2D matrix: time X (num_grps*num_boxes)
            
            disp('NOTE: using MATLAB ODE solver')

            % step 13a: solve the ODE ---------------------------------------------
            fname_ECOTRANode                = 'f_ECOTRANode_05262022';              % name of the MATLAB ODE solver
            fname_PhysicalFlux_intraODE     = 'f_PhysicalFlux_intraODE_09092019';   % name of the MATLAB ODE solver sub-function

            disp(['Running MATLAB solver: ' fname_ECOTRANode])

            tic
            [T, ProductionRates] = ode23t(@(t, ProductionRates_t) f_ECOTRANode_05262022(t, ProductionRates_t, ODEinput), t_grid, production_initial(:)); % run MATLAB ODE solver; use this for getting soln at each time-point; NOTE: uses Thornton-Lessem scaler
            time_ODE = toc

            store_T                                             = T;
            store_ProductionRates(:, 1:(num_grps*num_boxes))	= ProductionRates; % (mmole N/m3/d)
            % ---------------------------------------------------------------------

            % step 14b: unstack result to recover variable structure of spatial boxes
            re_Y    = reshape(store_ProductionRates, [num_t, num_grps, num_boxes]); % (mmole N/m3/d); (3D matrix: time X groups X num_boxes)
        % end (case 'MatlabSolver') ---------------------------------------
            
    end % (switch switch_ODEsolver) 
    % *********************************************************************
    
    
    
    
    
    % *********************************************************************
    %% STEP 15: save run results--------------------------------------------
    
    %% step 15a: RUNlog of called functions --------------------------------
    RUNlog.fname_ECOTRANdynamic             = fname_ECOTRANdynamic;
    RUNlog.BiologicalModel                  = BiologicalModel_name;
    RUNlog.PhysicalModel                    = ECOTRANphysics.fname_PhysicalModel;
    RUNlog.MigrationModel                   = ECOTRANmigration.fname_ECOTRANmigration;
    RUNlog.DVMmodel                         = DVM.fname_DVMmodel;
    RUNlog.fname_ReadEwE                    = dat.fname_ReadEwE;
    RUNlog.fname_AggregateBiologicalModel	= EwEResult.fname_AggregateBiologicalModel;
    RUNlog.fname_ECOTRANheart               = ECOTRAN.fname_ECOTRANheart;
    RUNlog.fname_ECOfunction                = ECOTRAN.fname_ECOfunction;                % name of this version of f_ECOfunction
%     RUNlog.fname_RedistributeCannibalism	= ECOTRAN.fname_RedistributeCannibalism;	% FFF this function does not yet pass its name along; name of this version of f_RedistributeCannibalism
    RUNlog.fname_calcEE                     = ECOTRAN.fname_calcEE;                     % file name of this f_calcEEsub-function
%     RUNlog.fname_CalcPredationBudget        = ECOTRAN.fname_CalcPredationBudget;        % FFF this function does not yet pass its name along; file name of the sub-function f_CalcPredationBudget
%     RUNlog.fname_E2Epedigree                = ECOTRAN_PEDIGREE.fname_E2Epedigree;       % name of this f_E2Epedigree function
%     RUNlog.fname_E2E_MonteCarlo             = ECOTRAN_MC.fname_E2E_MonteCarlo;        	% name of this f_E2E_MonteCarlo function
    RUNlog.fname_LightIntensity             = ECOTRANphysics.fname_LightIntensity;      % name of physical model sub-function
    RUNlog.fname_EvaluateFluxBalance        = ECOTRANphysics.fname_EvaluateFluxBalance;	% name of this f_EvaluateFluxBalance sub-function
    RUNlog.fname_UnCompactFluxTimeSeries	= ECOTRANphysics.fname_UnCompactFluxTimeSeries; % name of this f_UnCompactFluxTimeSeries sub-function
    RUNlog.fname_CalcNetFlux                = ECOTRANphysics.fname_CalcNetFlux;          % name of this f_CalcNetFlux function
    RUNlog.fname_CompactFlux                = CompactFlux_ADVECTION.fname_CompactFlux;            % file name of f_CompactFluxTimeSeries function
%     RUNlog.fname_StaticProductionTimeseries	= fname_StaticProductionTimeseries;         % (not yet converted to ECOTRAN II) name of f_StaticProductionTimeseries function
    RUNlog.fname_FunctionalResponse_MonteCarlo = fname_FunctionalResponse_MonteCarlo;	% name of this f_FunctionalResponse_MonteCarlo function
% 	  RUNlog.fname_SeasonalPB_PrimaryProducers	= fname_SeasonalPB_PrimaryProducers;	% (deactivated); name of f_SeasonalPB_PrimaryProducers function
%     RUNlog.fname_MichaelisMenten            = fname_MichaelisMenten;                    % (optional); name of this f_MichaelisMenten function
    RUNlog.fname_InitialProductionRates     = fname_InitialProductionRates;             % QQQ name of f_InitialProductionRates function
    RUNlog.fname_WebProductivity            = 'f_WebProductivity_03272019';             % SSS name of f_WebProductivity function
    RUNlog.fname_ThorntonLessem             = fname_ThorntonLessem;                     % name of f_ThorntonLessem function
    RUNlog.fname_ECOTRANode                 = fname_ECOTRANode;
    RUNlog.fname_PhysicalFlux_intraODE      = fname_PhysicalFlux_intraODE;
	RUNlog.fname_VarianceDivision           = 'f_VarianceDivision_12132018';            % SSS name of f_VarianceDivision function
	RUNlog.fname_VarianceMultiplication     = 'f_VarianceMultiplication_12132018';      % SSS name of f_VarianceMultiplication function
    RUNlog.time_ODE                         = time_ODE;                                 % time to run ODE (mins?)
    % ---------------------------------------------------------------------
    
    
    % step 15b: save run results ------------------------------------------
    save(saveFile, 'RUNlog', 're_Y', 'store_T', 'PHYSICSinput', 'ECOTRANphysics', 'ODEinput'); % save model run
    % *************************************************************************

% end % (MonteCarlo_loop) % QQQ MC loop is off for debugging


% end m-file***************************************************************