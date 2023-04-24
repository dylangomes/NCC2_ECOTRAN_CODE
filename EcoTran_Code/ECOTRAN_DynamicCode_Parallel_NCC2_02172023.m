function ECOTRANdynamic_NCC2_09262022_batch(setWD,Model_name,START,END,Region,upwelling_driver,CUTI_LAT,CUTI_YEARS,run_Treatments,switch_FoodWebScenario,switch_SubModel,switch_INITIALproduction,switch_MonteCarlo,num_MC,FileOffset,ShowOutput)
% a function version of ECOTRANdynamic for batch runs; run a dynamic model over time
%
% calls:
  %       f_readEwEcsv_10pp_07072021                	read in ECOPATH (EwE) model from VisualBasic .csv file and store as variable 'dat'; (use for VisualBasic food web files with 10 primary producers)
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
%       f_WebProductivity_03272019
%       f_ScenarioGenerator_10212020                for food web restructuring scenario
%
%       alternate physical models:
  %           f_OrdinalDate                         	calculate ordinal dates; day of year with January 1 of any year = 1
  %           f_ECOTRANphysics_CNP_06172022           NCC 5-box cross-shelf ecosystem: define the physical geometry of the model and prepare physical advection, mixing, and sinking rate time-series
  %               f_read_NWPO3data_02012022           read wind data from station NWPO3 on Newport South Jetty (no smoothing)	(m/s)                           (vertical vector)
  %               f_UpwellingIndex                    calculat Bakun Upwelling Index from wind data                           (m3/s per 100m of coastline)	(vertical vector)
  %               f_smooth                            mean smoothing over provided smoothing_window (in days before and after time-point)
  %               f_read_ERDdata_05122021             calculate daily median upwelling intensities from ERD products (http://http://www.pfeg.noaa.gov/products/PFELData/upwell/6_hourly/upwell45N125W)
  %               [deactivated] f_ekmandepth          calculate ekman depth
  %               f_prep_ERD_CUTI_08052022           process ERD CUTI timeseries for upwelling flux cropped to datestart:dateend time period	(m2/s)	(vertical vector: num_t X 1)
  %                   f_smooth                            mean smoothing over provided smoothing_window (in days before and after time-point)
  %               round2                              round to a specified number of decimal places
  %               f_CompactFluxTimeSeries_11182019	compact flux time series when arranged as (3D matrix: time X source box X destiny box)
  %               f_UnCompactFluxTimeSeries_12112019  UnCompact a flux time series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole
  %                   f_calcNetFlux_12112019              calculate net flux into and net flux out of each model box and across outer domain boundaries
  %               f_EvaluateFluxBalance_11262021      examine for flux time-series for imbalances IN & OUT of invidual boxes and IN & OUT of the overall domain
  %               calcur_res.mat                      (dataset with monthly mean salinity, temperature, and NO3+NO2 data)
  %               f_LightIntensity_12112020           instantaneous (W/m2), daily mean averaged across 24 h (W m^-2 h^-1), & daily integrated (W m^-2 d^-1) solar raditation at ocean surface; (vertical vector: num_t X 1)
  %
  %       f_ECOTRANmigration_NCC_11252020             area of overlap of neighboring model sub-domains
  %
  %       f_CompactFluxTimeSeries_11182019            compact physical flux 3D matrix
  %
  %       f_DVMsinusoid_02192021                      Calculate DVM flux rates between all model domain boxes for each functional group at each time point; Uses a daily sinusoidal migration pattern
  %
  %       f_FunctionalResponse_MonteCarlo_05132021    prepare array of vulnerability terms and allows for random generation of functional response terms within a predefined uncertainty level
  %
  %       f_SeasonalPB_PrimaryProducers               (deactivated)
  %       f_MichaelisMenten_05152016                  (optional) phytoplankton uptake rate of NO3 & NH4 (mmole N/m3/d) & p/b @ t as calculated from Michaelis-Menten uptake kinetics
  %
  %       f_InitialProductionRates_02012022          	calculate initial or mean production conditions
  %           f_WebProductivity_03272019                calculate production rates of all groups under a given driver (e.g., NO3 or primary production); also accounts for defined rates of group production export when running static scenarios
  %
  %   2 ODE solver options:
    %       f_PrepMexODE_2D_08212022                  	prepare ECOTRAN variables for using C++ ODE solver mex function. Pack parameters & drivers along proper dimensions
  %           f_unspoolMATRIX_04282020                  linearize multidimenional matrices up to 4-D for use in C++
    %       mex_ECOTRANode_2D_08222022                 	solve the ecosystem ODE in C++
    %
  %       f_ECOTRANode_2D_08202022                 	Ordinary Differential Equation (use this for getting soln at each time-point; default is for reflective boundary, but has built-in options for defined boundary conditions)
  %           f_PhysicalFlux_intraODE_09092019          calculate biomass fluxes into and out of each box and across domain boundaries (mmoles N/m3/d)	(3D matrix: 1 X num_grps X num_boxes DESTINATION)
  %           f_MichaelisMenten_05152016                (optional) Michaelis-Menton uptake kinetics for phytoplankton
  %
  % returns:
    %       store_T                     time-series of day umbers
  %       store_ProductionRates       time-series of production rates for each functional group; (t/km2/d) 
  %
  % NOTE: code cannot currently accomodate seasonal changes in flows to non-terminal detritus pools (e.g., fisheries and detritus columns of EnergyBudget)
  %
  % to compile (MUST have boost installed):
    % mex mex_QQQ.cpp -I/usr/local/boost_1_75_0/include/
    %
  CHANGEGROUP = run_Treatments(1);
  SCALE = run_Treatments(2);
  
  % detritus recycling parameters
  current_benthicDetritusSequestration	  = 0.1;
  current_pelagicDetritusMetabolism       = 0.1;
  current_benthicDetritusMetabolism       = 0.1;
  
  test_budget = (current_benthicDetritusMetabolism + current_benthicDetritusSequestration);
  if test_budget > 0.99
  % if this is true, do sequestration first, then adjust metabolism to leave a minimum 1% available to predation
  current_benthicDetritusMetabolism = 0.99 - current_benthicDetritusSequestration;
  
  if ShowOutput
  disp('NOTICE: benthic detritus sequestration & metabolism use exceeeds 99%. Metabolism will be reduced.')   
  end
  end
  
  current_benthicDetritusPredation        = 1 - (current_benthicDetritusMetabolism + current_benthicDetritusSequestration);
  
  
  % *************************************************************************
    % STEP 1: load & aggregate EwE results-------------------------------------
    % step 1a: set operating conditions ---------------------------------------
    switch_ODEsolver            = 'CppSolver';                          % OPTION 1: C++
    % switch_ODEsolver            = 'MatlabSolver';                       % OPTION 2: Matlab
  
  % switch_FunctionalResponse	= 'Linear';             % linear functional response; NOTE: STRICTLY DONER-DRIVEN DYNAMICS (rate of consumption by each consumer is a direct proportion of the production by each of its prey groups)
  switch_FunctionalResponse	= 'NonLinear_default';	% NonLinear_default functional response
  % switch_FunctionalResponse	= 'NonLinear_alt';      % NonLinear_default functional response
  
  switch_PhysicalModel        = '2D_physics';
  % switch_PhysicalModel        = '3D_ROMS';
  
  
  % -------------------------------------------------------------------------
    
    
% step 1b: define food web model to use -----------------------------------
ReadFile_directory      = setWD;  % QQQ activate for Jim's directory

BiologicalModel_name	= Model_name; % Dylan's post-heatwave model
  
readFile            	= [ReadFile_directory BiologicalModel_name];
  
% step 1c: load ECOPATH (EwE) model from Aydin VisualBasic file (.csv format)
dat                  	= f_readEwEcsv_10pp_07072021(readFile);	% use for models with up to 10 primary producers
  
% step 1d: aggregate model results & prep EwEResult for analysis ----------
[EwEResult, PEDIGREE] 	= f_AggregateBiologicalModel_02052021(dat);
  
  % step 1e: define filename and directory for saving results ---------------
    SaveFile_directory      = setWD; 
  SaveFile_label          = strcat(regexprep(Model_name,".csv",""),"_", num2str(CHANGEGROUP,'%03d'),"_", num2str(SCALE)); 
  % *************************************************************************
    
    
    
    
    % *************************************************************************
    % STEP 2: ECOTRAN conversion-----------------------------------------------
    MonteCarloStore         = [];
  [ECOTRAN]            	= ECOTRANheart_09032021(EwEResult, MonteCarloStore);
  % *************************************************************************
    
    
    
    
    
    % *************************************************************************
    % STEP 3: Generate E2E Monte Carlo models based on ECOTRAN EnergyBudget----
    %         Start with the one original ECOTRAN base model and generate a set of 
  %           Monte Carlo models from the ECOTRAN EnergyBudget & ConsumptionBudget
  %           matrices using predefined CV values
  
  switch switch_MonteCarlo
  
  case 'MonteCarlo_build' % generate (and save) a stack of Monte Carlo models
  
  
  if ShowOutput
  disp(['MonteCarlo: building stack of ' num2str(num_MC) ' food webs'])
  end
  
  ECOTRAN.num_MC      = num_MC;
  
  PEDIGREE.ee_eggs_CV                               = 0.01; % SSS egg pedigree W.R.T. production budget for all groups; (CV); (scaler)
  PEDIGREE.BacterialMTBLSM_CV                       = 0.01; % SSS pedigree for implicit bacterial metabolism of terminal detritus (CV)
  PEDIGREE.Oxidation_NH4_CV                         = 0.01;	% fraction of NH4 produced oxidized directly back to NO3 abiologically; QQQ scaler?? (vertical vector: num_NH4 X 1)??
  PEDIGREE.NutrientUptake_CV                        = 0.01; % SSS pedigree for nutrient uptake by primary producers (CV)
  ECOTRAN_PEDIGREE                                  = f_E2Epedigree_08042020(ECOTRAN, PEDIGREE); % NEW!!!
    
    % SSS use for standardized pedigree
  %     overwrite the pedigree values from the ECOPATH (EwE) model from VisualBasic file (.csv format)
  [rows, clms]                                = size(ECOTRAN_PEDIGREE.EnergyBudget_CV);
  ECOTRAN_PEDIGREE.EnergyBudget_CV            = 0.001 * ones(rows, clms); % changes how important predators are relative to eachother
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV       = zeros(7, clms);
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV(1, :) = 0.05; % feces
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV(2, :) = 0.05; % metabolism
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV(3, :) = 0.05; % eggs
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV(4, :) = 0.05; % predation
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV(5, :) = 0.05; % senescence
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV(6, :) = 0.05; % ba
  ECOTRAN_PEDIGREE.ConsumptionBudget_CV(7, :) = 0.05; % em
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
  %         if ShowOutput
  %         disp(['MonteCarlo: SAVING stack: ' filename_MC])
  %         end
  %         save(filename_MC, 'ECOTRAN_MC')
  
  % end (case build_MonteCarlo) -----------------------------------------
    
    
    case 'MonteCarlo_load' % load a set of Monte Carlo models
  filename_MC = [SaveFile_directory 'MonteCarlo_NCC_stack_27-Jun-2022.mat']; % SSS be sure to give correct saved file name here
  if ShowOutput
  disp(['MonteCarlo: LOADING stack: ' filename_MC])
  end
  load(filename_MC, 'ECOTRAN_MC')
  
  num_MC                              = ECOTRAN_MC.num_MC;
  EnergyBudget_MC                     = ECOTRAN_MC.EnergyBudget_MC;
  ConsumptionBudget_MC                = ECOTRAN_MC.ConsumptionBudget_MC;
  DiscardFraction_MC                  = ECOTRAN_MC.DiscardFraction_MC;
  % end (case load_MonteCarlo) ------------------------------------------
    
    
    case 'MonteCarlo_TypeModel'
  
  if ShowOutput
  disp('MonteCarlo: using the defining TypeModel')
  end
  
  num_MC                          = 1;       % only the "type" model is used
  EnergyBudget_MC                 = ECOTRAN.EnergyBudget;
  ConsumptionBudget_MC            = ECOTRAN.ConsumptionBudget;
  DiscardFraction_MC              = ECOTRAN.DiscardFraction;
  
  ECOTRAN_MC.num_MC               = num_MC;
  ECOTRAN_MC.EnergyBudget_MC      = ECOTRAN.EnergyBudget; % (3D matrix: num_grps (consumers) X num_grps (producers) X num_MC (1))
  ECOTRAN_MC.ConsumptionBudget_MC = ECOTRAN.ConsumptionBudget; % (3D matrix: 7 X num_grps (producers) X num_MC (1))
  ECOTRAN_MC.DiscardFraction_MC	  = ECOTRAN.DiscardFraction;
  
  ECOTRAN_MC.fate_metabolism      = ECOTRAN.fate_metabolism;	% (3D matrix: num_nutrients X num_grps (producers) X num_MC (1))
  ECOTRAN_MC.fate_eggs            = ECOTRAN.fate_eggs;        % (3D matrix: num_eggs X num_grps (producers) X num_MC (1))
  ECOTRAN_MC.fate_feces           = ECOTRAN.fate_feces;       % (3D matrix: num_ANYdetritus X num_grps (producers) X num_MC (1))
  ECOTRAN_MC.fate_senescence      = ECOTRAN.fate_senescence;  % (3D matrix: num_ANYdetritus X num_grps (producers) X num_MC (1))
  ECOTRAN_MC.fate_predation       = ECOTRAN.fate_predation;   % (3D matrix: num_livingANDfleets X num_grps (producers) X num_MC (1))
  % end (case MonteCarlo_TypeModel) -------------------------------------
    
    end % (switch_MonteCarlo)
  % *************************************************************************
    
    
    
    
    
    % *************************************************************************
    % STEP 4: prep and pack ECOTRAN model parameters---------------------------
    
    % step 4a: read in ECOTRAN structure variables ----------------------------
    %          (so that no changes are made to original values)
  GroupType                         = ECOTRAN.GroupType;
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
  num_grps                          = ECOTRAN.num_grps;             % number of model groups
  % num_MC                              = ECOTRAN.num_MC;               % number of Monte Carlo models
  % TransferEfficiency                  = ECOTRAN.TransferEfficiency;	  % gets redefined manually below
  % -------------------------------------------------------------------------
    
    
    %% step 4b: find detritus, nutrients, ba & em -----------------------------
    %           row addresses in EnergyBudget_MC
  %           NOTE: ba = biomass accumulation term, em = emigration term
  
  looky_NO3                       	= find(GroupType        == ECOTRAN.GroupTypeDef_NO3);
  looky_plgcNH4                     = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4);
  looky_bnthNH4                     = find(GroupType        == ECOTRAN.GroupTypeDef_bnthNH4);
  looky_NH4                         = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4 | GroupType == ECOTRAN.GroupTypeDef_bnthNH4);
  looky_nutrients                   = find(floor(GroupType)	== ECOTRAN.GroupTypeDef_ANYNitroNutr);	% row addresses of nutrients
  looky_ANYPrimaryProducer        	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYPrimaryProd);
  looky_macroalgae                  = find(GroupType        == ECOTRAN.GroupTypeDef_Macrophytes);
  looky_phytoplankton               = looky_ANYPrimaryProducer(~ismember(looky_ANYPrimaryProducer, looky_macroalgae));
  looky_ANYconsumer                 = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYConsumer);
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
  num_NO3                           = length(looky_NO3);
  num_NH4                          	= length(looky_NH4);
  num_plgcNH4                       = length(looky_plgcNH4); % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
  num_bnthNH4                       = length(looky_bnthNH4); % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
  num_ANYPrimaryProd                = length(looky_ANYPrimaryProducer);
  num_phytoplankton                	= length(looky_phytoplankton);
  num_macroalgae                   	= length(looky_macroalgae);
  num_ANYconsumers                	= length(looky_ANYconsumer);
  num_fleets                        = length(looky_fleets);
  num_predators                     = num_ANYconsumers + num_fleets;
  num_livingANDfleets               = length(looky_livingANDfleets);
  num_eggs                         	= length(looky_eggs);
  num_ANYdetritus                   = length(looky_ANYdetritus);
  % -------------------------------------------------------------------------
    
    
    % step 4c: pack variables for ODE solver ----------------------------------
  ODEinput.looky_nutrients              = looky_nutrients;
  ODEinput.looky_NO3                    = looky_NO3;
  ODEinput.looky_NH4                    = looky_NH4;
  ODEinput.looky_plgcNH4                = looky_plgcNH4;
  ODEinput.looky_bnthNH4                = looky_bnthNH4;
  ODEinput.looky_ANYPrimaryProducer     = looky_ANYPrimaryProducer;
  ODEinput.looky_phytoplankton          = looky_phytoplankton;
  ODEinput.looky_macroalgae             = looky_macroalgae;
  ODEinput.looky_fleets                 = looky_fleets;
  ODEinput.looky_ANYconsumer            = looky_ANYconsumer;
  ODEinput.looky_livingANDfleets        = looky_livingANDfleets;
  ODEinput.looky_eggs                   = looky_eggs;
  ODEinput.looky_terminalPLGCdetritus	= looky_terminalPLGCdetritus;
  ODEinput.looky_terminalBNTHdetritus	= looky_terminalBNTHdetritus;
  ODEinput.looky_ANYdetritus            = looky_ANYdetritus;
  ODEinput.looky_nonNO3             	= looky_nonNO3;
  
  ODEinput.num_grps                     = num_grps;
  ODEinput.num_nutrients                = num_nutrients;
  ODEinput.num_NO3                      = num_NO3;
  ODEinput.num_NH4                      = num_NH4;
  ODEinput.num_plgcNH4              	= num_plgcNH4; % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
  ODEinput.num_bnthNH4              	= num_bnthNH4; % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
  ODEinput.num_ANYPrimaryProd           = num_ANYPrimaryProd;
  ODEinput.num_phytoplankton            = num_phytoplankton;
  ODEinput.num_macroalgae               = num_macroalgae;
  ODEinput.num_ANYconsumers             = num_ANYconsumers;
  ODEinput.num_predators                = num_predators;
  ODEinput.num_livingANDfleets          = num_livingANDfleets;
  ODEinput.num_eggs                 	= num_eggs;
  ODEinput.num_ANYdetritus              = num_ANYdetritus;
  % -------------------------------------------------------------------------
    
    
    %% *************************************************************************
% STEP 5: define TransferEfficiency terms----------------------------------
%         SET MANUALLY: TransferEfficiency = 1 for all groups because we are
%                       working with the consumption budget matrix that tracks the fate of ALL consumption (not just predation)
%                       system losses (groups where TransferEfficiency < 1 include benthic detritus, fishery, (others?)
%         NOTE: terminal benthic detritus column in EnergyBudget sums to 1 and has NO OtherMortality, 
%               TE is the only means of removing terminal benthic detritus from system
TransferEfficiency                              = ones(num_MC, num_grps);	% re-initialize as rows of horizontal vector of ones
                                         
% % QQQ CHANGE THIS DURING BATCH PROCESSING FOR f-ratio
% TransferEfficiency(:, find(contains(ECOTRAN.label,'benthic detritus')))         = 0.1; % set TE for terminal benthic detritus to < 1
% if ShowOutput
%disp('QQQ NOTICE: LOWERED TE OF BENTHIC DETRITUS')
%end
                                         
% *************************************************************************


                                           
                                           
% *************************************************************************
% SCENARIO STEP: Force changes to the food web structure via Static Scenario method
switch switch_FoodWebScenario
                                         
case 'FoodWebScenario_ON' % force changes to food web before dynamic run
                                         
% step 6a: define TransferEfficiency terms ------------------------
%         SET MANUALLY: TransferEfficiency = 1 for all groups because we are
%                       working with the consumption budget matrix that tracks the fate of ALL consumption (not just predation)
%                       system losses (groups where TransferEfficiency < 1 include benthic detritus, fishery, (others?)
%          NOTE: terminal benthic detritus column in EnergyBudget sums to 1 and has NO OtherMortality, 
%                TE is the only means of removing terminal benthic detritus from system
TransferEfficiency_scenario                 = ones(num_MC, num_grps); % re-initialize as rows of horizontal vector of ones
TransferEfficiency_scenario(:, [looky_terminalBNTHdetritus]) = 0.1; % NOTE: I use 0.1 TE for terminal benthic detritus as a standard default (JRuz 9/27/2018)
% -----------------------------------------------------------------
                                                                                  
                                                                                  
% step 6b: define production loss fractions------------------------
%         (FFF in future, might use vector read in from .csv file; but for now, this is all handled here)
LX_base                           	= ones(1, num_grps);
ProductionLossScaler_scenario       = zeros(1, num_grps); % initialize
PhysicalLossFraction_scenario       = ProductionLossScaler_scenario .* LX_base;
% -----------------------------------------------------------------

                                                                                  
% step 6c: initialize InputProductionVector------------------------
LrgPhy=find(contains(ECOTRAN.label,'large phytoplankton'));
SmPhy=find(contains(ECOTRAN.label,'small phytoplankton'));

InputProductionVector_scenario                    = zeros(1, num_grps);   % initialize InputProductionVector_scenario
InputProductionVector_scenario(LrgPhy)            = pb(LrgPhy)	.* biomass(LrgPhy);	% (t/km2/y); (vertical vector: num_grps X 1)
InputProductionVector_scenario(SmPhy)             = pb(SmPhy)	.* biomass(SmPhy);	% (t/km2/y); (vertical vector: num_grps X 1)
% -----------------------------------------------------------------
                                                                                  
                                                                                  
% step 6d: calculate base productivity-----------------------------
EnergyBudget_MC_temp                = EnergyBudget_MC; % this is used only to calculate static production
TestInput_scenario                  = InputProductionVector_scenario;
TestInput_scenario(looky_nutrients) = [];

if max(TestInput_scenario > 0) 
EnergyBudget_MC_temp(:, looky_nutrients, :)	= 0; % shut off recycling; % deactivate nutrient recycling; set nutrient columns in EnergyBudget to 0 so that recycled nutrients don't flow into any other group
end
        
Production_scenario                                 = zeros(num_grps, num_MC); % initialize; (2D matrix: num_grps X num_MC)
        for MonteCarlo_loop = 1:num_MC    
            current_EnergyBudget                    = EnergyBudget_MC_temp(:, :, MonteCarlo_loop);
            current_TransferEfficiency              = TransferEfficiency_scenario(MonteCarlo_loop, :);
            Production_scenario(:, MonteCarlo_loop)	= f_WebProductivity_03272019(current_TransferEfficiency, current_EnergyBudget, InputProductionVector_scenario, PhysicalLossFraction_scenario); % production (actually CONSUMPTION) rates as "initial" or "mean" conditions; (mmole N/m3/d); (2D matrix: num_grps X num_MC)
        end % (MonteCarlo_loop)
        % -----------------------------------------------------------------
        
        
        % step 6e: force change(s) to the food web ------------------------

        ScenarioConditions.modify_consumer         	= [CHANGEGROUP];           % row number(s) of consumer group(s) to force-modify
        ScenarioConditions.ScaleFactor            	= SCALE;                     	% value > 1 means increase flow to modify_consumer (and reduced flow to offset_consumer); < 1 means decrease flow to modify_consumer; 
        %                                                                               NOTE: keep ScaleFactor >= 0 (negative value makes no sense)
        ScenarioConditions.offset_consumer       	= [1:num_grps];                 % row number(s) of consumer group(s) to modify as offset to force-modified grp(s)
        %                                                                               NOTE: change in modify_consumer grp(s); can be [] if offset is to be distributed among ALL consumer groups
        ScenarioConditions.offset_consumer([looky_nutrients; looky_ANYdetritus]) = [];	% (assume no change to detritus or to nutrients; --YOU ARE FREE TO CHANGE THIS ASSUMPTION--!!!)
        ScenarioConditions.trgt_producer            = 1:num_grps;                   % column number(s) of producer group(s)
        % ScenarioConditions.EnergyBudget_base   	  = EnergyBudget_MC;            % MonteCarlo set of EnergyBudget_base models; (3D matrix: num_grps X num_grps X num_MC); 
        % %                                                                             NOTE: no changes made directly to EnergyBudget_base, everything is done via the fates

        ScenarioConditions.ConsumptionBudget_base	= ConsumptionBudget_MC;         % MonteCarlo set of ConsumptionBudget_base models; (3D matrix: 7 X num_grps X num_MC)

        ScenarioConditions.fate_metabolism_base     = ECOTRAN_MC.fate_metabolism;	% fate of metabolism in base model; (3D matrix: num_nutrients X num_grps X num_MC); NOTE: using fate_metabolism with nutrient recycling deactivated
        ScenarioConditions.fate_eggs_base           = ECOTRAN_MC.fate_eggs;         % fate of eggs (reproduction) in base model; (3D matrix: num_eggs X num_grps X num_MC)
        ScenarioConditions.fate_feces_base      	= ECOTRAN_MC.fate_feces;        % fate of feces detritus in base model; (3D matrix: num_ANYdetritus X num_grps X num_MC)
        ScenarioConditions.fate_predation_base   	= ECOTRAN_MC.fate_predation; 	% fate of production among all predators in base model; (3D matrix: num_livingANDfleets X num_grps X num_MC); NOTE: using fate_predation with nutrient recycling deactivated
        ScenarioConditions.fate_senescence_base   	= ECOTRAN_MC.fate_senescence;	% fate of senescence detritus in base model; (3D matrix: num_ANYdetritus X num_grps X num_MC)

        % ScenarioConditions.DiscardFraction_base     = DiscardFraction_MC;           % Monte Carlo set of fractions of catch discarded for each group and fleet; (3D matrix: num_grps X num_fleets X num_MC)
        ScenarioConditions.DiscardFraction_base     = repmat(ECOTRAN.DiscardFraction, [1, 1, num_MC]);	% QQQ "type" model until I fix algorithm for adjusting individual >1 fractions; Monte Carlo set of fractions of catch discarded for each group and fleet; (3D matrix: num_grps X num_fleets X num_MC)

        ScenarioConditions.TransferEfficiency    	= TransferEfficiency;           % (2D matrix: num_MC X num_grps)
        ScenarioConditions.InputProductionVector   	= InputProductionVector_scenario;	% NOTE: use SCENARIO input production driver rates TO REDUCE primary production; (horizontal vector: 1 X num_grps)
        ScenarioConditions.PhysicalLossFraction    	= PhysicalLossFraction_scenario;	% (horizontal vector: 1 X num_grps)

        StaticScenario_results_1                   	= f_ScenarioGenerator_10212020(ScenarioConditions, ECOTRAN);
        % -----------------------------------------------------------------


        % step 6f: replace EnergyBudget_MC & ConsumptionBudget_MC with the modified food web
        EnergyBudget_MC             = StaticScenario_results_1.EnergyBudget_scenario;
        ConsumptionBudget_MC        = StaticScenario_results_1.ConsumptionBudget_scenario;
        % -----------------------------------------------------------------
        
	% end (case 'FoodWebScenario_ON') -------------------------------------

end % (switch_FoodWebScenario)
% *************************************************************************





%% *************************************************************************
% STEP 6: prepare physical parameters--------------------------------------
% step 6a: input time-frame info ------------------------------------------
%          NOTE: prepared ERD_BUI .csv file covers: 01-Jan-1969 through 30-Sep-2015 and has 6-hr resolution
%          NOTE: ERD_CUTI file covers:              28-Feb-1990 through 28-Feb-2022  and has 1-day resolution
%          NOTE: NH Line temperature covers:        05-May-1997 through 08-Jul-2021 and has 1-day resolution

if ShowOutput
disp('--->>> 2D upwelling physics')
end

% ERD_CUTI dates (QQQ - QQQ)
datestart                       = datenum(START); % SSS --> enter starting date
dateend                         = datenum(END); % SSS --> enter ending date

dt                              = 24/24; % t-step; (days); (dt = 1 = 1 d; dt = 4/24 = 4 hours)
                                  % NOTE: take care to select good dt values for diel vertical migration 
                                  %       (other values do not scale well between 1 & -1 in sin diel cycle (probably due to rounding error of pi() function)
datestart_OrdinalDate           = f_OrdinalDate(datestart);
min_t                           = datestart_OrdinalDate;
max_t                           = (dateend - datestart) + datestart_OrdinalDate;
t_grid                          = linspace(min_t, (max_t+1-dt), ((dateend - datestart + 1)/dt))'; % QQQ NEW VERSION!!!; (vertical vector); t_grid runs from datestart to dateend and is inclusive of dateend; intervals = dt
num_t                           = length(t_grid); % length of t_grid; (scaler)
PHYSICSinput.datestart          = datestart;
PHYSICSinput.dateend            = dateend;
PHYSICSinput.dt                 = dt;
PHYSICSinput.t_grid             = t_grid;

% auto select latitude
if strcmp(CUTI_LAT,"Auto")
if strcmp(Region, 'WA')
CUTI_LAT = 47;
elseif strcmp(Region, 'CR')
CUTI_LAT = 46;
elseif strcmp(Region, 'NOR')
CUTI_LAT = 45;
elseif strcmp(Region, 'SOR')
CUTI_LAT = 43;
elseif strcmp(Region, 'NCA')
CUTI_LAT = 41;
else
  if ShowOutput
disp('ERROR: Choose a valid subregion: (WA,CR,NOR,SOR,NCA)')
end
end
end
CUTI_LAT
PHYSICSinput.target_latitudes	= CUTI_LAT;
PHYSICSinput.smoothing_window	= 5; % window for smoothing before and after time-point (e.g., 2 = a window of 5 days centered on time-point t)

spatial_BiomassScalers         	= [1 1 0 1 0];	% NCC scalers for estimating initial (or mean) primary producer biomasses across model domain; NOTE: these values are assumed; NOTE: x2 in Box I used to compensate for 30m depth relative to 15 m depths in Boxes II & IV; FFF apply NPZD scalers here
% -------------------------------------------------------------------------


  % step 6b: prepare advection & mixing time-series for each model box ------
  % 2D upwelling driver

ECOTRANphysics                  = f_ECOTRANphysics_NCC2_upwelling_09042022(PHYSICSinput, upwelling_driver,ShowOutput,CUTI_YEARS); % SSS specify physical flux time-series to use: 'Brink_BUI', 'NWPO3_BUI', 'ERD_BUI', 'ERD_CUTI', or 'Fake_Upwelling')


 % step 6c: migration of each group (shared box face areas) ----------------
   %          NOTE: code does not handle migration flux uncertainty (nor physical flux uncertainty)
 % ECOTRANmigration            = f_ECOTRANmigration_NCC_11252020(ECOTRANphysics);
     ECOTRANmigration            = f_ECOTRANmigration_NCC_2D_09042022(ECOTRANphysics);

 % ODEinput.biomass_migrator = ECOTRANmigration.biomass_migrator;  % SSS special definition of boundary biomasses for migrators; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes); NOTE: de-comment migrator biomass lines within ODE code
 % -------------------------------------------------------------------------


   %% step 6d: compact fluxes -------------------------------------------------
   %           remove information defining non-existing box links
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
  
  ADVECTION                       = ECOTRANphysics.ADVECTION;         % volume advected per day from source box to destination box; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
  HORIZONTALMIXING                = ECOTRANphysics.HORIZONTALMIXING;	% volume mixed horizontally per day from source box to destination box; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
  VERTICALMIXING                  = ECOTRANphysics.VERTICALMIXING;	% volume mixed vertically per day from source box to destination box; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
  SINKING                         = ECOTRANphysics.SINKING;           % box floor area between sinking source box and destination box; (m2); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
  MIGRATION                       = ECOTRANmigration.MIGRATION;       % migration as boundary area between boxes; (m2); (3D matrix: num_t X SOURCE (num_boxes+1 boundary) X DESTINY (num_boxes+1 boundary))
  
  CompactFlux_ADVECTION           = f_CompactFluxTimeSeries_11182019(ADVECTION,ShowOutput);        % compact ADVECTION 3D matrix
  CompactFlux_HORIZONTALMIXING	= f_CompactFluxTimeSeries_11182019(HORIZONTALMIXING,ShowOutput); % compact HORIZONTALMIXING 3D matrix
  CompactFlux_VERTICALMIXING      = f_CompactFluxTimeSeries_11182019(VERTICALMIXING,ShowOutput);   % compact VERTICALMIXING 3D matrix
  CompactFlux_SINKING             = f_CompactFluxTimeSeries_11182019(SINKING,ShowOutput);          % compact SINKING 3D matrix
  CompactFlux_MIGRATION           = f_CompactFluxTimeSeries_11182019(MIGRATION,ShowOutput);        % compact MIGRATION 3D matrix
  
  
  % QQQ new 8/13/2022 (put this directly in the physics code)
  % SinkLink_surface                 = ECOTRANphysics.SinkLink_surface;                 % define sinking connections from surface to sub-surface boxes; (unitless); (2D matrix: destination boxes X source boxes)
  % SinkLink_surface                 = reshape(SinkLink_surface, [num_boxes, 1, num_boxes]);    % (3D matrix: destination boxes X 1 X source boxes)
  % SinkLink_benthos                 = ECOTRANphysics.SinkLink_benthos;                         % define sinking connections from sub-surface boxes to benthic detritus; (unitless); (2D matrix: destination boxes X source boxes)
  % SinkLink_benthos                 = reshape(SinkLink_benthos, [num_boxes, 1, num_boxes]);    % (3D matrix: destination boxes X 1 X source boxes)
  
  ODEinput.SinkLink_surface = ECOTRANphysics.SinkLink_surface; % QQQ only used for 5-box, 2D cross-shelf settings; (3D matrix: destination boxes X 1 X source boxes)
  ODEinput.SinkLink_benthos = ECOTRANphysics.SinkLink_benthos; % (3D matrix: destination boxes X 1 X source boxes)
  % -------------------------------------------------------------------------
    
    
    %% step 6e: unpack & process physics variables -----------------------------
  num_boxes                    	= ECOTRANphysics.num_boxes;
  BoxVolume                    	= ECOTRANphysics.BoxVolume;     	% (m3); (2D matrix: num_t X num_boxes)
  % BoxLength                    	= ECOTRANphysics.BoxLength;        	% (m); (2D matrix: num_t X num_boxes)
  % BoxHeight                     = ECOTRANphysics.BoxHeight;       	% (m); (2D matrix: num_t X num_boxes)
  % BoxWidth                     	= ECOTRANphysics.BoxWidth;       	% (m); (2D matrix: num_t X num_boxes)
  
  Io                            	= ECOTRANphysics.Io;             	% time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (vertical vector: num_t X 1)
  current_light                   = ECOTRANphysics.current_light;     % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day; (vertical vector: num_t X 1)
  sunrise                         = ECOTRANphysics.sunrise;           % time of sunrise (time in h from midnight; 12 = noon, set as a default)
  sunset                          = ECOTRANphysics.sunset;            % time of sunset  (time in h from midnight; 12 = noon, set as a default)
  Kw                            	= ECOTRANphysics.Kw;              	% Light attenuation_seawater; (scalar)
  Kp                             	= ECOTRANphysics.Kp;              	% Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N); (scalar)
  MLD                           	= ECOTRANphysics.MLD;            	% mixed-layer depth; (m); (vertical vector: num_t X 1)
  EuphoticDepth                  	= repmat(MLD, [1 num_boxes]);     	% FFF (eventually move to f_physics code); depth of euphotic zone, used when converting the vertically-integrated EwE primary producer biomass to biomass/volume; depth; (m); (2D matrix: num_t X num_boxes)
  
  NO3timeseries_conc              = ECOTRANphysics.NO3timeseries_conc; % NO3 + NO2 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
  % NH4timeseries_conc              = ECOTRANphysics.NH4timeseries_conc; % NH4 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
  NO3initial_rate                 = ECOTRANphysics.NO3initial_rate;	 % initial NO3 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))
  % NH4initial_rate                 = ECOTRANphysics.NH4initial_rate;	 % initial NH4 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))
  
  WWT_to_C                      	= ECOTRANphysics.WWT_to_C;          % (scalar)
  atomic_mass_C                  	= ECOTRANphysics.atomic_mass_C;     % (scalar)
  C_to_N_phytoplankton          	= ECOTRANphysics.C_to_N_phytoplankton; % (scalar)
  
  grp_row                      	= 1:num_grps;
  % % -----------------------------------------------------------------------
    
    
    %% step 6f: pack physics variables for ODE solver -------------------------
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
    
    
    %% step 6g: define sinking rate of individual groups (m/d) ----------------
    %          NOTE: sinking is treated like a physical term (i.e., not incorporated into EnergyBudget)
  %          NOTE: apply this factor whether or not using Michaelis-Menten for primary producers
  SinkingSpeed                                = zeros(num_t, num_grps, num_fluxes_sinking);	% initialze sinking speed time-series; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
  
  % % QQQ deactivate sinking for debugging
  if ShowOutput
  disp('QQQ NOTICE: SINKING IS ACTIVE ONLY FOR PLGC_DETRITUS')
  end
  SinkingSpeed(:, find(contains(ECOTRAN.label,'pelagic detritus')), :)      	= repmat((10.5), [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
  % SinkingSpeed(:, rc_bnth_detritus, :)       	= repmat((25),   [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
  % SinkingSpeed(:, rc_fishery_offal, :)       	= repmat((40),   [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
  
  % MichaelisMenten_w                               = [0.6 1.0];                                    % SSS; sinking speed; (m/d); [Sm Phytoplankton, Lg Phytoplankton]
  % SinkingSpeed(:, looky_ANYPrimaryProducer, :)	= repmat(MichaelisMenten_w, [num_t, 1, num_fluxes_sinking]);	% sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
  % -------------------------------------------------------------------------
    
    
    %% step 6h: define migration speeds & duration for each group (m/d) -------
    %          NOTE: dusk phase are negative (in "standard" DVM)
  DVMinput.MigrationArea_compact          = MigrationArea_compact; % migration as boundary area between boxes; (m2); (3D matrix: num_t X num_fluxes_migration)
  DVMinput.MigrationSpeed                 = zeros(num_grps, (num_boxes+1), (num_boxes+1)); % initialize; (m/d); (3D matrix: num_grps, source (num_boxes+1) X destiny (num_boxes+1))
  
  DVMspeed_meso2epi                       = zeros(num_grps, 2); % initialize; DVM speed MESO<<-->>EPI; (2D matrix: num_grps X 2 [dusk dawn]);
  DVMspeed_bathy2meso                     = zeros(num_grps, 2); % initialize; DVM speed BATHY<<-->>MESO; (2D matrix: num_grps X 2 [dusk dawn]);
  
  % % SSS set GoMexOcn DVMspeeds 
  % %       NOTE: DVM speed calculations done outside of ECOTRAN (see file "DVMcalibration_GoMexOcn_08232021.xlsx")
  % %       NOTE: first term = dusk migration speed; second term = dawn migration speed
  % DVMspeed_meso2epi(largeCopepods_RC, 1:2)            = [340    50]; % (m/d)
  % DVMspeed_meso2epi(smallCopepods_RC, 1:2)            = [732    53]; % (m/d)
  
  % DVMspeed_bathy2meso(WM_mesopelagics_RC, 1:2)        = [1119	 566]; % (m/d)
  
  
  % % MESO<<-->>EPI migration
  % DVMinput.MigrationSpeed(:, 1, 2)           = DVMspeed_meso2epi(:, 2); % dawn migration, (epi-->>meso)
  % DVMinput.MigrationSpeed(:, 2, 1)           = DVMspeed_meso2epi(:, 1) * (-1); % dusk migration, (meso-->>epi), noon to midnight; NOTE: should be negative
  % 
  % % BATHY<<-->>MESO migration
  % DVMinput.MigrationSpeed(:, 2, 3)           = DVMspeed_bathy2meso(:, 2); % dawn migration, (meso-->>bathy)
  % DVMinput.MigrationSpeed(:, 3, 2)           = DVMspeed_bathy2meso(:, 1) * (-1); % dusk migration, (bathy-->>meso), noon to midnight; NOTE: should be negative
  
  
  DVM                         = f_DVMsinusoid_08122021(ECOTRANphysics, ODEinput, DVMinput, t_grid); % calculate Diel Vertical Migration terms; (structure)
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
  % STEP QQQ: physiological modification terms-------------------------------
  
  % define default place-holder values
Q10_scaler              = []; % modifies metabolism time-series in ConsumptionBudget (modifications all done within this script)
q_TemperatureScaler     = ones(num_t, num_grps, num_boxes); % modify consumption rates (q) within the ODE solver

%% *************************************************************************
  
  
  
  
  
  %% *************************************************************************
  % STEP 7: prepare external_driver time-series------------------------------
  %         NOTE: This is the external input that enters the model domain via 
%               advection & mixing.
%         NOTE: DEFAULT: The reflective boundary assumption is that biomass of non-external_driver 
%               groups are the same on either side of model domain outer boundaries.

% step 7a: define external driver group(s) & driver time-series -----------
  %          NOTE: driver time-series (external_driver) is a biomass density and not a rate (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
%          FFF: defining for all boxes now, but will try to trim down to just the boxes with fluxes

looky_driver                            = looky_NO3; % identify driver group(s)
num_drivers                             = length(looky_driver);
NO3timeseries_conc                      = reshape(NO3timeseries_conc, [num_t, 1, num_boxes]);	% (mmole N/m3); (3D matrix: num_t X 1 X num_boxes)
% NH4timeseries_conc                      = reshape(NH4timeseries_conc, [num_t, 1, num_boxes]);	% (mmole N/m3); (3D matrix: num_t X 1 X num_boxes)

external_driver(:, looky_NO3, :)     	= NO3timeseries_conc; % SSS; external boundary driver (e.g. NO3) biomass for each box; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes)
% external_driver(:, looky_plgcNH4, :)	= NH4timeseries_conc; % SSS; external boundary driver (e.g. NH4) biomass for each box; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes)

external_driver(:, :, (end+1))          = 0; % add layer for external boundary driver (e.g. NO3) biomass for each box; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
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
if ShowOutput
disp('NOTE: DONOR-DRIVEN FUNCTIONAL RESPONSE')
end
% end (case 'Linear') ----------
  
  case 'NonLinear_default' % non-linear default
if ShowOutput
disp('NOTE: DEFAULT NON-LINEAR FUNCTIONAL RESPONSE')
end
% end (case 'NonLinear_default') ----------
  
  case 'NonLinear_alt' % constant (predation independent of predator biomass)
FunctionalResponseParams	= FunctionalResponseParams * 2;     % SSS force to linear for testing & default; (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)
if ShowOutput
disp('NOTE: ALTERNATE NON-LINEAR FUNCTIONAL RESPONSE')
end
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
  
  
  
  
  
  
  % QQQ *********************************************************************
  % QQQ apply batch tuning adjustments to the ConsumptionBudget (only at top of MC stack)
% EnergyBudget terms will be adjusted to match these changes in the ODE solver
ConsumptionBudget_MC(2, looky_terminalBNTHdetritus, :) = current_benthicDetritusMetabolism;
ConsumptionBudget_MC(4, looky_terminalBNTHdetritus, :) = current_benthicDetritusPredation;
ConsumptionBudget_MC(5, looky_terminalBNTHdetritus, :) = current_benthicDetritusSequestration;
     
% adjust pelagic Detritus ConsumptionBudget rate by modifying ConsumptionBudget immediately before STEP 10 as follows:
% set pelagic detritus predation  = (1 - current_pelagicDetritusMetabolism) * (ConsumptionBudget_MC(4, looky_terminalPLGCdetritus, 1) / (ConsumptionBudget_MC(4, looky_terminalPLGCdetritus, 1) + ConsumptionBudget_MC(5, looky_terminalPLGCdetritus, 1))
% set pelagic detritus senescence = (1 - current_pelagicDetritusMetabolism) * (ConsumptionBudget_MC(5, looky_terminalPLGCdetritus, 1) / (ConsumptionBudget_MC(4, looky_terminalPLGCdetritus, 1) + ConsumptionBudget_MC(5, looky_terminalPLGCdetritus, 1))
ConsumptionBudget_MC(2, looky_terminalPLGCdetritus, :)	= current_pelagicDetritusMetabolism;
original_predation                                      = ConsumptionBudget_MC(4, looky_terminalPLGCdetritus, 1);
original_senescence                                     = ConsumptionBudget_MC(5, looky_terminalPLGCdetritus, 1);
ConsumptionBudget_MC(4, looky_terminalPLGCdetritus, :)  = (1 - current_pelagicDetritusMetabolism) * (original_predation / (original_predation + original_senescence));
ConsumptionBudget_MC(5, looky_terminalPLGCdetritus, :)  = (1 - current_pelagicDetritusMetabolism) * (original_senescence / (original_predation + original_senescence));

% QQQ set nitrification of plgcNH4 to 20% following Yool et al 2007
ConsumptionBudget_MC(2, looky_plgcNH4, :)   = 0.2; % nitrification
ConsumptionBudget_MC(4, looky_plgcNH4, :)   = 0.8; % 1- nitrification
if ShowOutput
disp('QQQ NOTICE: pelagic NH4 conversion to NO3 set to 20% & predation set to 80% in ConsumptionBudget')
end
% QQQ *********************************************************************
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  

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
% if ShowOutput
%display('-->WARNING in ECOTRAN_DynamicScenario: primary producer predation & senescence terms in ConsumptionBudget are redefined with Michaelis-Menten eta values.')
%end
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

% use this for 5-box 2D physics
EnergyBudget_BoxType        = ones(num_grps, num_grps, num_boxes);	% (3D matrix: num_grps X num_grps X num_boxes)
ConsumptionBudget_BoxType	= ones(7, num_grps, num_boxes);         % (3D matrix: 7 X num_grps X num_boxes)
looky_SubSurfaceBoxes       = [3 5];

EnergyBudget_BoxType(:,                 :,                          looky_SubSurfaceBoxes)	= 0; % by default, allow no trophic interactions in sub-surface boxes; (3D matrix: num_grps X num_grps X num_boxes)
EnergyBudget_BoxType(looky_NO3,         looky_NH4,                  looky_SubSurfaceBoxes)	= 1; % allow oxidation of NH4 to NO3 (nitrification); (3D matrix: num_grps X num_grps X num_boxes)
EnergyBudget_BoxType(looky_NH4,         looky_ANYdetritus,          looky_SubSurfaceBoxes)	= 1; % allow implicit bacterial metabolism of detritus; (3D matrix: num_grps, num_grps, num_boxes)
EnergyBudget_BoxType(looky_ANYdetritus, looky_ANYPrimaryProducer,   looky_SubSurfaceBoxes)	= 1; % allow senescence flow of primary producers to detritus; (3D matrix: num_grps, num_grps, num_boxes)
EnergyBudget_BoxType(looky_ANYdetritus, looky_ANYdetritus,          looky_SubSurfaceBoxes)	= 1; % allow senescence flow between detritus pools (including flow of terminal pelagic detritus to terminal benthic detritus); (3D matrix: num_grps, num_grps, num_boxes)

ConsumptionBudget_BoxType(:,            :,                          looky_SubSurfaceBoxes)	= 0; % by default, allow no trophic interactions in sub-surface boxes; (3D matrix: 7 X num_grps X num_boxes)
ConsumptionBudget_BoxType(2,            looky_NH4,                  looky_SubSurfaceBoxes)	= 1; % allow oxidation of NH4 to NO3 (nitrification); QQQ check logic of this to match EnergyBudget_BoxType
ConsumptionBudget_BoxType(2,            looky_ANYdetritus,          looky_SubSurfaceBoxes)	= 1; % allow implicit bacterial metabolism of detritus
ConsumptionBudget_BoxType(5,            looky_ANYPrimaryProducer,   looky_SubSurfaceBoxes)	= 1; % allow senescence flow of primary producers to detritus
ConsumptionBudget_BoxType(5,            looky_ANYdetritus,          looky_SubSurfaceBoxes)	= 1; % allow senescence flow between detritus pools (including flow of terminal pelagic detritus to terminal benthic detritus)
ConsumptionBudget_BoxType(7,            looky_terminalBNTHdetritus,	looky_SubSurfaceBoxes)	= 1; % allow for sequestration of benthic detritus via emigration (em) in sub-surface boxes

% QQQ for 2D cross-shelf physics ONLY
EnergyBudget_BoxType(:,                 looky_bnthNH4,                                :)	= 0; % QQQ set all trophic flows of bnthNH4 to 0 in all boxes; bnthNH4 will get pooled into plgcNH4 and that is when it becomes available to the food web
if ShowOutput
disp('QQQ NOTICE: no trophic outflow of bnthNH4 in any boxes, EnergyBudget_BoxType rows in bnthNH4 column are set to 0')
end
ConsumptionBudget_BoxType(:,            looky_bnthNH4,                                :)	= 0; % QQQ set all trophic flows of bnthNH4 to 0 in all boxes; bnthNH4 will get pooled into plgcNH4 and that is when it becomes available to the food web
if ShowOutput
disp('QQQ NOTICE: no trophic outflow of bnthNH4 in any boxes, ConsumptionBudget_BoxType rows in bnthNH4 column are set to 0')
end
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
% saveFile            = strcat(SaveFile_directory,'Output/',SaveFile_label,'_',date,'.mat');

for MonteCarlo_loop = 1:num_MC
if ShowOutput
display(['MonteCarlo run ' num2str(MonteCarlo_loop) ' of ' num2str(num_MC)])
end
saveFile                        = strcat(SaveFile_directory,'Output/temp/',SaveFile_label,'_', num2str(MonteCarlo_loop+FileOffset,'%04d'),'_',date,'.mat'); 

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
  
  if ShowOutput
disp('NOTE: IDENTICAL regional biological models')
end

current_biomass                 = repmat(current_biomass, [1, num_boxes]);                          % (t WWT/km2); (2D matrix: num_grps X num_boxes)
current_biomass                 = current_biomass .* repmat(spatial_BiomassScalers, [num_grps 1]);	% (t WWT/km2); (2D matrix: num_grps X num_boxes)

current_pb                      = repmat(current_pb,                    [1, 1, num_boxes]);	% specific growth rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)
current_qb                      = repmat(current_qb,                    [1, 1, num_boxes]);	% specific consumption rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)

current_EnergyBudget            = repmat(current_EnergyBudget,          [1, 1, num_boxes]); % use for a single spatial definition of the EnergyBudget; (3D matrix: num_grps X num_grps X num_boxes)
current_EnergyBudget            = current_EnergyBudget .* EnergyBudget_BoxType;             % adjust sub-surface EnergyBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: num_grps X num_grps X num_boxes)

current_ConsumptionBudget       = repmat(current_ConsumptionBudget,     [1, 1, num_boxes]); % use for a single spatial definition of the ConsumptionBudget; (3D matrix: 7 X num_grps X num_boxes)
current_ConsumptionBudget       = current_ConsumptionBudget .* ConsumptionBudget_BoxType;	% adjust sub-surface ConsumptionBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: 7 X num_grps X num_boxes)

current_fate_feces              = repmat(current_fate_feces,            [1, 1, num_boxes]); % (3D matrix: num_ANYdetritus X num_grps X num_boxes)
current_fate_metabolism         = repmat(current_fate_metabolism,       [1, 1, num_boxes]); % (3D matrix: num_nutrients X num_grps X num_boxes)
current_fate_eggs               = repmat(current_fate_eggs,             [1, 1, num_boxes]); % (3D matrix: num_eggs X num_grps X num_boxes)
% NOTE fate_predation is calculated below
current_fate_senescence         = repmat(current_fate_senescence,       [1, 1, num_boxes]); % (3D matrix: num_ANYdetritus X num_grps X num_boxes)

current_TransferEfficiency      = repmat(current_TransferEfficiency,    [1, 1, num_boxes]); % (3D matrix: 1 X num_grps X num_boxes)
if ShowOutput
disp('QQQ NOTICE: only have bnth detritus TE < 1 in sub-surface boxes?')
end
current_FunctionalResponseParams      = repmat(current_FunctionalResponseParams,    [1, 1, num_boxes]);	% (3D matrix: CONSUMERS (num_grps) X prey group (num_grps) X num_boxes) replicated across clms (= producers)

current_MichaelisMenten_Vmax	= repmat(current_MichaelisMenten_Vmax,  [1, 1, num_boxes]);	% Vmax = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_KNO3	= repmat(current_MichaelisMenten_KNO3,  [1, 1, num_boxes]);	% KNO3 = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_KNH4	= repmat(current_MichaelisMenten_KNH4,  [1, 1, num_boxes]);	% KNH4 = NH4 half-saturation constant; (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_alpha	= repmat(current_MichaelisMenten_alpha, [1, 1, num_boxes]);	% alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_psi     = repmat(current_MichaelisMenten_psi,   [1, 1, num_boxes]);	% psi = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_w       = repmat(current_MichaelisMenten_w,     [1, 1, num_boxes]); % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_eta     = repmat(current_MichaelisMenten_eta,   [1, 1, num_boxes]);	% eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
% end case 'identical_SubModel' -----------------------------------
  
  case 'independent_SubModel'  % OPTION 2: use independently defined webs for each shelf zone ----
  %           use for multiple definitions of the EnergyBudget
% FFF future code could easily allow for Monte Carlo models for each sub-region

if ShowOutput
disp('NOTE: LOADING independently defined regional food web models')
end

% load invidividual regional models
if strcmp(Region, 'WA')
SR = [1 2 3];             
elseif strcmp(Region, 'CR')
SR = [4 5 6];
elseif strcmp(Region, 'NOR')
SR = [7 8 9];
elseif strcmp(Region, 'SOR')
SR = [10 11 12];
elseif strcmp(Region, 'NCA')
SR = [13 14 15];
else 
  if ShowOutput
disp('ERROR: Choose a valid subregion: (WA,CR,NOR,SOR,NCA)')
end
end   

if ShowOutput
disp(strcat("Region selected: ",Region ," = zones ", num2str(SR)))
end

SubRegion1 = NCC2_SubregionReadIn_12072022(ECOTRAN,SR(1));
SubRegion2 = NCC2_SubregionReadIn_12072022(ECOTRAN,SR(2));
SubRegion3 = SubRegion2;
SubRegion4 = NCC2_SubregionReadIn_12072022(ECOTRAN,SR(3));
SubRegion5 = SubRegion4;


current_EnergyBudget                = repmat(current_EnergyBudget, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_grps X num_grps X num_boxes)
current_EnergyBudget(:, :, 1)       = SubRegion1.EnergyBudget;
current_EnergyBudget(:, :, 2)       = SubRegion2.EnergyBudget;
current_EnergyBudget(:, :, 3)       = SubRegion3.EnergyBudget;
current_EnergyBudget(:, :, 4)       = SubRegion4.EnergyBudget;
current_EnergyBudget(:, :, 5)       = SubRegion5.EnergyBudget;
current_EnergyBudget                = current_EnergyBudget .* EnergyBudget_BoxType; % adjust sub-surface EnergyBudget matrices; allow only certain processes in sub-surface boxes; (fraction); (3D matrix: num_grps X num_grps X num_boxes)

current_ConsumptionBudget           = repmat(current_ConsumptionBudget, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: 7 X num_grps X num_boxes)
current_ConsumptionBudget(:, :, 1)       = SubRegion1.ConsumptionBudget;
current_ConsumptionBudget(:, :, 2)       = SubRegion2.ConsumptionBudget;
current_ConsumptionBudget(:, :, 3)       = SubRegion3.ConsumptionBudget;
current_ConsumptionBudget(:, :, 4)       = SubRegion4.ConsumptionBudget;
current_ConsumptionBudget(:, :, 5)       = SubRegion5.ConsumptionBudget;

% see step 11a
current_ConsumptionBudget(7, looky_terminalBNTHdetritus, :)	= current_ConsumptionBudget(7, looky_terminalBNTHdetritus, :) + current_ConsumptionBudget(5, looky_terminalBNTHdetritus, :); % emigration = 7 and senescence = 5 -> both go to emigration
current_ConsumptionBudget(5, looky_terminalBNTHdetritus, :)	= 0; % set senescence to 0
            
current_ConsumptionBudget         	= current_ConsumptionBudget .* ConsumptionBudget_BoxType; % adjust sub-surface ConsumptionBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: 7 X num_grps X num_boxes)

current_biomass                     = repmat(current_biomass, [1, num_boxes]);	% replicate box 1 model to initialize; (mmoles N/m3); NOTE: these are INITIAL biomass conditions; (2D matrix: num_grps X num_boxes); NOTE: nutrients are zeros; NOTE: only the definition for primary producer biomass is used, it is used to drive initial model production values
current_biomass(:, 1)               = SubRegion1.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;
current_biomass(:, 2)               = SubRegion2.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;
current_biomass(:, 3)               = SubRegion3.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;
current_biomass(:, 4)               = SubRegion4.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;
current_biomass(:, 5)               = SubRegion5.biomass; % NOTE: conversion to (mmoles N/m3) was already done manually;

current_pb                          = zeros(1, num_grps, num_boxes); % re-initialize; specific growth rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)
current_pb(:, :, 1)                 = SubRegion1.pb;
current_pb(:, :, 2)               	= SubRegion2.pb;
current_pb(:, :, 3)                 = SubRegion3.pb;
current_pb(:, :, 4)                 = SubRegion4.pb;
current_pb(:, :, 5)                 = SubRegion5.pb;
            
current_qb                          = zeros(1, num_grps, num_boxes); % re-initialize; specific consumption rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)
current_qb(:, :, 1)                 = SubRegion1.qb;
current_qb(:, :, 2)               	= SubRegion2.qb;
current_qb(:, :, 3)                 = SubRegion3.qb;
current_qb(:, :, 4)                 = SubRegion4.qb;
current_qb(:, :, 5)                 = SubRegion5.qb;
            
current_fate_feces                  = repmat(current_fate_feces, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_ANYdetritus X num_grps X num_boxes)
current_fate_metabolism             = repmat(current_fate_metabolism, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_nutrients X num_grps X num_boxes)
current_fate_eggs                   = repmat(current_fate_eggs, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_eggs X num_grps X num_boxes)
current_fate_senescence             = repmat(current_fate_senescence, [1, 1, num_boxes]); % replicate box 1 model to initialize; (3D matrix: num_ANYdetritus X num_grps X num_boxes)

current_TransferEfficiency       = repmat(current_TransferEfficiency, [1, 1, num_boxes]);  	% (3D matrix: 1 X num_grps X num_boxes)        
current_FunctionalResponseParams = repmat(current_FunctionalResponseParams, [1, 1, num_boxes]);	% (3D matrix: CONSUMERS X prey group X num_boxes) replicated across clms (= producers)

current_MichaelisMenten_Vmax	= repmat(current_MichaelisMenten_Vmax, [1, 1, num_boxes]);	% Vmax = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_KNO3	= repmat(current_MichaelisMenten_KNO3, [1, 1, num_boxes]);	% KNO3 = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_KNH4	= repmat(current_MichaelisMenten_KNH4, [1, 1, num_boxes]);	% KNH4 = NH4 half-saturation constant; (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_alpha	= repmat(current_MichaelisMenten_alpha, [1, 1, num_boxes]);	% alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_psi     = repmat(current_MichaelisMenten_psi, [1, 1, num_boxes]);	% psi = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_w       = repmat(current_MichaelisMenten_w, [1, 1, num_boxes]);     % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
current_MichaelisMenten_eta     = repmat(current_MichaelisMenten_eta, [1, 1, num_boxes]);	% eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)



% end case{'independent_SubModel'} --------------------------------
  
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
% ---------------------------------------------------------------------
  
  
  % % step 12e: OPTIONAL: make scenario changes to ConsumptionBudget ----
  % %           NOTE: just a test scenario
% ConsumptionBudget_ba(10000:20000, rc_fleet_midWater_trawls, 1)              = 0.5;    % QQQ NCC: high mortality for 1 month
% ConsumptionBudget_senescence(10000:20000, rc_fleet_midWater_trawls, 1)      = 0.1307; % QQQ NCC: high mortality for 1 month
% ConsumptionBudget_predation(950:1000, rc_anchovy, 1)                        = 0.1136; % QQQ NCC: high mortality for 1 month
% 
% ConsumptionBudget_ba(10000:20000, rc_fleet_midWater_trawls, 4)              = 0.9;    % QQQ NCC: high mortality for 1 month
% ConsumptionBudget_ba(950:1000, str2num(SPECIES), 1)                         = 0.5;    % QQQ NCC: high mortality for 1 month
% % -------------------------------------------------------------------
  
  
  %% step 12f: pack variables needed for ODE ----------------------------

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

ODEinput.q_TemperatureScaler            = q_TemperatureScaler; % Thornton-Lessem temperature effects for adjustment to consumption rate; value between 0 and 1; (3D matrix: num_t X num_grps X num_boxes)
% *********************************************************************
  
  
  
  
  
  % *********************************************************************
  %% STEP 13: calculate INITIAL production rate conditions---------------
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
[production_initial, fname_InitialProductionRates]	= f_InitialProductionRates_02012022(ODEinput, production_initial_driver, t_initial,ShowOutput);  % initial or mean production rates (actually consumption inflow); (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE: code does NOT make transfer of bnthNH4 from surface to sub-surface boxes

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
[production_initial, fname_InitialProductionRates] = f_InitialProductionRates_02012022(ODEinput, production_initial_driver, t_initial,ShowOutput);  % initial or mean production rates (actually consumption inflow); (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE: code does NOT make transfer of bnthNH4 from surface to sub-surface boxes

% paste in initial NO3 & NH4 input rates
production_initial(looky_NO3, :)      	= NO3initial_rate;   % append initial NO3 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
  % production_initial(looky_plgcNH4, :)    = NH4initial_rate;   % append initial NH4 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
  % end (case 'INITIALproduction_pb') -------------------------------
  
  case 'INITIALproduction_SubModel'	% METHOD 3: use for driving initial model conditions with values loaded along with regional sub-model definitions 
% For example, use for oceanic GoMexOcn, CNP models----
  % CBB: This is where we add in initial conditions from a long run from a different file (save end 5-10 yrs re_Y average consumption rate).

if ShowOutput
disp('USING INITIAL CONDITIONS FROM EXTERNAL FILE: ')
end

fname_InitialProductionRates	= 'pre-defined in independent sub-model';
production_initial              = [SubRegion1.production_initial SubRegion2.production_initial SubRegion3.production_initial SubRegion4.production_initial SubRegion5.production_initial]; % initial q; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes)

% % paste in initial NO3 & NH4 input rates
% production_initial(looky_NO3, :)      	= NO3initial_rate;   % append initial NO3 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
  % production_initial(looky_plgcNH4, :)    = NH4initial_rate;   % append initial NH4 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
  
  % end case (INITIALproduction_SubModel) ---------------------
  
  
  case 'INITIALproduction_nutrients'	% METHOD 4: use for driving initial model conditions with mean annual nutrient input rates
%           NOTE: NH4 uptake is deactivated in f_InitialProductionRates_02012022, so production_initial will be lower for METHOD 4 than for METHOD 2 (the latter implicitly includes recylced primary production)

if ShowOutput
disp('INITIAL CONDITIONS: calculated from nutrient input rate')
end

production_initial_driver                           = zeros(1, num_grps, num_boxes);        % initialize DriverProductionVector; (3D matrix: 1 X num_grps X num_boxes)
production_initial_driver(1, looky_NO3, :)          = NO3initial_rate;	% plug in NO3 input rate; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
% production_initial_driver(1, looky_plgcNH4, :)      = NH4initial_rate;	% plug in NH4 input rate; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)

% finalize initial condition calculations
[production_initial, fname_InitialProductionRates]	= f_InitialProductionRates_02012022(ODEinput, production_initial_driver, t_initial,ShowOutput);  % QQQ NH4 uptake turned OFF; initial or mean production rates (actually consumption inflow); (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE: code does NOT make transfer of bnthNH4 from surface to sub-surface boxes

% paste in initial NO3 & NH4 input rates
production_initial(looky_NO3, :)      	= NO3initial_rate;   % append initial NO3 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
  % production_initial(looky_plgcNH4, :)    = NH4initial_rate;   % append initial NH4 input rates to NO3 rows; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
  % end (case 'INITIALproduction_nutrients') ------------------------
  
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
% NOTE: >> mex mex_ECOTRANode_11272020.cpp -I/usr/local/include/ % command to compile mex function

if ShowOutput
disp('NOTE: using C++ ODE solver')
end

% step 13a: prepare ECOTRAN variables for using C++ ODE solver mex function
%           Pack parameters & drivers along proper dimensions.
if ShowOutput
disp('prepare variables for C++ ODE solver...')
end

[AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct] = f_PrepMexODE_2D_08212022(ODEinput); % pack variables to send to ODE solver QQQ needs ThorntonLessem

fname_ECOTRANode                 = 'mex_ECOTRANode_2D_08222022';	% ODE solver name
fname_PhysicalFlux_intraODE      = 'NA';	% SSS            
% -------------------------------------------------------------
  
  % step 13b: run the model in C++ ------------------------------
  if ShowOutput
disp(['Running C++ solver: ' fname_ECOTRANode])
end

tic
[output_Cpp]            = mex_ECOTRANode_2D_08222022(AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct); % run the ODE solver
time_ODE                = toc

store_T                 = output_Cpp.mat_obs_times; % (d); (vertical vector: num_t X 1)
store_ProductionRates	= output_Cpp.mat_obs_states; % (mmole N/m3/d); (2D matrix: num_t X (num_grps*num_boxes))
% -------------------------------------------------------------
  
  % step 13c: unstack result to recover variable structure of spatial boxes
re_Y    = reshape(store_ProductionRates, [num_t, num_grps, num_boxes]); % (mmole N/m3/d); (3D matrix: time X groups X num_boxes)
% end (case 'CppSolver') ------------------------------------------
  
  
  case 'MatlabSolver' % Matlab solver
% NOTE: calculation uses ode23t (trial-and-error suggests a bit better performance than ODE45)
% NOTE: size of ProductionRates is 2D matrix: time X (num_grps*num_boxes)

if ShowOutput
disp('NOTE: using MATLAB ODE solver')
end

% step 13a: solve the ODE ---------------------------------------------
  fname_ECOTRANode                = 'f_ECOTRANode_2D_08202022';         	% name of the MATLAB ODE solver
fname_PhysicalFlux_intraODE     = 'f_PhysicalFlux_intraODE_09092019';	% name of the MATLAB ODE solver sub-function

if ShowOutput
disp(['Running MATLAB solver: ' fname_ECOTRANode])
end

tic
[T, ProductionRates] = ode23t(@(t, ProductionRates_t) f_ECOTRANode_2D_08202022(t, ProductionRates_t, ODEinput), t_grid, production_initial(:)); % run MATLAB ODE solver; use this for getting soln at each time-point has corrected pb & qb usage; has transfer of bnthNH4 from surface to subsurface boxes (use for 2D only)
time_ODE = toc

store_T                                             = T;
store_ProductionRates(:, 1:(num_grps*num_boxes))	= ProductionRates; % (mmole N/m3/d)
% -------------------------------------------------------------
  
  % step 13b: unstack result to recover variable structure of spatial boxes
re_Y    = reshape(store_ProductionRates, [num_t, num_grps, num_boxes]); % (mmole N/m3/d); (3D matrix: time X groups X num_boxes)
% end (case 'MatlabSolver') ---------------------------------------
  
  end % (switch switch_ODEsolver) ---------------------------------------
  % *********************************************************************
  
  
  
  
  
  % *********************************************************************
  %% STEP 15: save run results--------------------------------------------
  
  %% step 15a: RUNlog of called functions --------------------------------
  RUNlog.fname_ECOTRANdynamic             = mfilename;
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
RUNlog.fname_E2E_MonteCarlo             = ECOTRAN_MC.fname_E2E_MonteCarlo;        	% name of this f_E2E_MonteCarlo function
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
RUNlog.fname_physiology_Q10             = 'not applied';                            % name of fname_physiology_Q10 function
RUNlog.fname_ThorntonLessem             = 'not applied';                            % name of f_ThorntonLessem function
RUNlog.fname_physiology_Kitchell      	= 'not applied';                            % name of f_physiology_Kitchell function
RUNlog.fname_ECOTRANode                 = fname_ECOTRANode;
RUNlog.fname_PhysicalFlux_intraODE      = fname_PhysicalFlux_intraODE;
RUNlog.fname_VarianceDivision           = 'f_VarianceDivision_12132018';            % SSS name of f_VarianceDivision function
RUNlog.fname_VarianceMultiplication     = 'f_VarianceMultiplication_12132018';      % SSS name of f_VarianceMultiplication function
RUNlog.time_ODE                         = time_ODE;                                 % time to run ODE (mins?)
% ---------------------------------------------------------------------
  
  
  % step 15b: save run results ------------------------------------------
  save(saveFile, 'RUNlog', 're_Y', 'store_T', 'PHYSICSinput', 'ECOTRANphysics', 'ODEinput'); % save model run
%   save(saveFile, 'RUNlog', 're_Y', 'store_T', 'PHYSICSinput', 'ECOTRANphysics', 'ECOTRANmigration', 'ODEinput', '-v7.3'); % save model run; NOTE: add option '-v7.3' for very large ODEinput cases

  % *********************************************************************
  
  end % (MonteCarlo_loop) ---------------------------------------------------
  
  
  % end m-file***************************************************************