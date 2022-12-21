function [re_Y, Q_cp, qb, pb] = f_calcQcp_result_09072022(read_file, subregion)
% calculate CONSUMPTION matrix (Q_cp) timeseries from the f_ECOTRANdynamic result (re_Y)
% revision date: 4-20-2022


% % *************************************************************************
% % STEP 1: define model run files & operating conditions--------------------
% 
% % step 1a: define run files -----------------------------------------------
% filedirectory           = '/Users/jimsebi/Documents/12_GoMex_project/6_FunctionalResponseTests/';
% run_file                = 'TestModel10_testAC_09-Mar-2022.mat';
% read_file               = [filedirectory run_file];
% % -------------------------------------------------------------------------
% 
% % step 1b: select sub-region of interest ----------------------------------
% subregion = 1;
% % *************************************************************************





% *************************************************************************
% STEP 2: load run results from run output file----------------------------

% step 2a: load current results fiel --------------------------------------
load(read_file, 're_Y', 'store_T', 'PHYSICSinput', 'ODEinput'); % load model run
% -------------------------------------------------------------------------


% step 2b: get time vectors, num_grps, num_boxes, q/b, physical fluxes ----
datestart                       = PHYSICSinput.datestart;
T                               = store_T(:,1,1);
min_t                           = T(1);
max_t                           = T(end);
num_t                           = length(T);
matlabdate                      = datevec(T + datestart - 1); % model dates as time vector; (2D matrix: num_t X [year month day hour minute second])
yearvector                      = (min(matlabdate(:, 1)):1:max(matlabdate(:, 1))); % (horizontal vector: 1 X num years)

num_grps                        = ODEinput.num_grps;
num_nutrients                   = ODEinput.num_nutrients;
num_eggs                        = ODEinput.num_eggs;
num_livingANDfleets             = ODEinput.num_livingANDfleets;
num_ANYdetritus                 = ODEinput.num_ANYdetritus;

looky_nutrients                 = ODEinput.looky_nutrients;
looky_eggs                      = ODEinput.looky_eggs;
looky_livingANDfleets           = ODEinput.looky_livingANDfleets;
looky_ANYdetritus               = ODEinput.looky_ANYdetritus;

qb                              = ODEinput.qb(:, :, subregion);         % biomass-specific consumption rate; (1/d); (2D matrix: num_t X num_grps)
pb                              = ODEinput.pb(:, :, subregion);         % biomass-specific production rate; (1/d); (2D matrix: num_t X num_grps)

FunctionalResponseParams        = ODEinput.FunctionalResponseParams(:, :, subregion);	% producer vulnerability (m_p); (2D matrix: CONSUMERS X prey group) replicated across clms (= producers)
q_TemperatureScaler             = ODEinput.q_TemperatureScaler(:, :, subregion);        % Thornton-Lessem temperature adjustment to consumption rate; value between 0 and 1; (2D matrix: num_t X num_grps)

EnergyBudget                    = ODEinput.EnergyBudget(:, :, subregion); % (2D matrix: num_grps X num_grps)

ConsumptionBudget_feces       	= ODEinput.ConsumptionBudget_feces(:, :, subregion); % (2D matrix: num_t X num_grps)
ConsumptionBudget_metabolism   	= ODEinput.ConsumptionBudget_metabolism(:, :, subregion); % (2D matrix: num_t X num_grps)
ConsumptionBudget_eggs         	= ODEinput.ConsumptionBudget_eggs(:, :, subregion); % (2D matrix: num_t X num_grps)
ConsumptionBudget_predation   	= ODEinput.ConsumptionBudget_predation(:, :, subregion); % (2D matrix: num_t X num_grps)
ConsumptionBudget_senescence  	= ODEinput.ConsumptionBudget_senescence(:, :, subregion); % (2D matrix: num_t X num_grps)
ConsumptionBudget_ba          	= ODEinput.ConsumptionBudget_ba(:, :, subregion); % (2D matrix: num_t X num_grps)
ConsumptionBudget_em          	= ODEinput.ConsumptionBudget_em(:, :, subregion); % (2D matrix: num_t X num_grps)

fate_feces                      = ODEinput.fate_feces(:, :, subregion);   	 % (2D matrix: num_ANYdetritus X num_grps)
fate_metabolism                 = ODEinput.fate_metabolism(:, :, subregion); % (2D matrix: num_nutrients X num_grps)
fate_eggs                       = ODEinput.fate_eggs(:, :, subregion);    	 % (2D matrix: num_eggs X num_grps)
fate_predation                  = ODEinput.fate_predation(:, :, subregion);  % (2D matrix: num_livingANDfleets X num_grps)
fate_senescence                 = ODEinput.fate_senescence(:, :, subregion); % (2D matrix: num_ANYdetritus X num_grps)

consumptionC_initial_repmat     = ODEinput.productionC_initial_repmat(:, :, subregion);	% initial consumer c consumption inflow rate; (mmole N/m3/d); (2D matrix: num_grps X num_grps); NOTE: replicated vertical vectors across columns

re_Y                            = re_Y(:, :, subregion); % consumption rates q of each group; (mmoles N/m3/d); (2D matrix: num_t X num_grps)
ConsumptionRates                = re_Y; % consumption rates q of each group; (mmoles N/m3/d); (2D matrix: num_t X num_grps)
% *************************************************************************





% *************************************************************************
% STEP 3: Adjust EnergyBudget to accomodate changes in ConsumptionBudget---
%         Changes in ConsumptionBudget are due to seasonal changes in physiology, migration, etc. 
%         SSS: deactivate this step if ConsumptionBudget does not change over time
%         NOTE: Box types are already accounted for
%         NOTE: Make no changes to senescence due to sinking. Senescence
%               directs biomass transfer to detritus (which is subject to
%               its own sinking rate). Sinking is an additional loss term 
%               handled as any other physical flux in dy calculation.
%         ConsumptionBudget:
%                       1) feces
%                       2) metabolism
%                       3) eggs (reproduction)
%                       4) predation
%                       5) senescence
%                       6) ba (biomass accumulation)
%                       7) em (emigration); NOTE: negative for immigration


% step 3a: replicate EnergyBudget along time-axis -------------------------
%          NOTE: time-axis lies along the 3rd dimension (easy, since we're not using a sub-region dimension here)
EnergyBudget_repmat = repmat(EnergyBudget, [1, 1, num_t]); % (3D matrix: num_grps X num_grps X num_t)
% -------------------------------------------------------------------------

% step 3b: reshape ConsumptionBudget terms -------------------------
%          NOTE: time-axis lies along the 3rd dimension
ConsumptionBudget_feces       	= reshape(ConsumptionBudget_feces',      [1, num_grps, num_t]); % (3D matrix: 1 X num_grps X num_t); NOTE transpose
ConsumptionBudget_metabolism	= reshape(ConsumptionBudget_metabolism', [1, num_grps, num_t]); % (3D matrix: 1 X num_grps X num_t); NOTE transpose
ConsumptionBudget_eggs       	= reshape(ConsumptionBudget_eggs',       [1, num_grps, num_t]); % (3D matrix: 1 X num_grps X num_t); NOTE transpose
ConsumptionBudget_predation     = reshape(ConsumptionBudget_predation',  [1, num_grps, num_t]); % (3D matrix: 1 X num_grps X num_t); NOTE transpose
ConsumptionBudget_senescence	= reshape(ConsumptionBudget_senescence', [1, num_grps, num_t]); % (3D matrix: 1 X num_grps X num_t); NOTE transpose
ConsumptionBudget_ba            = reshape(ConsumptionBudget_ba',         [1, num_grps, num_t]); % (3D matrix: 1 X num_grps X num_t); NOTE transpose
ConsumptionBudget_em            = reshape(ConsumptionBudget_em',         [1, num_grps, num_t]); % (3D matrix: 1 X num_grps X num_t); NOTE transpose
% -------------------------------------------------------------------------

% step 3c: replicate EnergyBudget along time-axis -------------------------
%          NOTE: time-axis lies along the 3rd dimension
fate_feces_repmat               = repmat(fate_feces,      [1, 1, num_t]);	% (3D matrix: num_ANYdetritus X num_grps X num_t)
fate_metabolism_repmat          = repmat(fate_metabolism, [1, 1, num_t]);	% (3D matrix: num_nutrients X num_grps X num_t)
fate_eggs_repmat                = repmat(fate_eggs,       [1, 1, num_t]);	% (3D matrix: num_eggs X num_grps X num_t)
fate_predation_repmat           = repmat(fate_predation,  [1, 1, num_t]);	% (3D matrix: num_livingANDfleets X num_grps X num_t)
fate_senescence_repmat          = repmat(fate_senescence, [1, 1, num_t]);	% (3D matrix: num_ANYdetritus X num_grps X num_t)
% -------------------------------------------------------------------------

% step 3d: ConsumptionBudget @ t: metabolism into EnergyBudget ------------
EnergyBudget_repmat(looky_nutrients, :, :)          = repmat(ConsumptionBudget_metabolism, [num_nutrients, 1, 1]) .* fate_metabolism_repmat; % excretion of NH4 and nitrification of NH4-->>NO3; (fraction); (3D matrix: num_nutrients X num_grps X num_t)
% -------------------------------------------------------------------------

% step 3e: ConsumptionBudget @ t: egg production into EnergyBudget --------
EnergyBudget_repmat(looky_eggs, :, :)               = repmat(ConsumptionBudget_eggs, [num_eggs, 1, 1]) .* fate_eggs_repmat; % egg production (should work for [] eggs or for multiple eggs); (fraction); (3D matrix: num_eggs X num_grps X num_t)
% -------------------------------------------------------------------------

% step 3f: ConsumptionBudget @ t: predation & fleet catch into EnergyBudget
EnergyBudget_repmat(looky_livingANDfleets, :, :)	= repmat(ConsumptionBudget_predation, [num_livingANDfleets, 1, 1]) .* fate_predation_repmat; % (fraction); (3D matrix: num_livingANDfleets X num_grps X num_t)
% -------------------------------------------------------------------------

% step 3g: ConsumptionBudget @ t feces & senescence into EnergyBudget -----
feces_t                                             = repmat(ConsumptionBudget_feces,      [num_ANYdetritus, 1, 1]) .* fate_feces_repmat; % (3D matrix: num_ANYdetritus X num_grps X num_t)
senescence_t                                        = repmat(ConsumptionBudget_senescence, [num_ANYdetritus, 1, 1]) .* fate_senescence_repmat; % (3D matrix: num_ANYdetritus X num_grps X num_t)
detritus_t                                          = feces_t + senescence_t;	% detritus = feces + senescence; (fraction); (3D matrix: num_ANYdetritus X num_grps X num_t)
EnergyBudget_repmat(looky_ANYdetritus, :, :)    	= detritus_t;               % detritus = feces + senescence; (fraction); (3D matrix: num_ANYdetritus X num_grps X num_t)
% *************************************************************************





% *************************************************************************
% STEP 4: calculate consumption rate matrix (Q_cp)-------------------------

% step 4a: replicate FunctionalResponseParams along the time-axis ---------
FunctionalResponseParams_repmat = repmat(FunctionalResponseParams, [1, 1, num_t]); % producer vulnerability (m_p); (3D matrix: CONSUMERS X prey group X num_t); replicated across clms (= producers)

% step 4b: reshape consumer & producer consumption rates ------------------
ConsumptionRatesC               = ConsumptionRates'; % transpose; (mmoles N/m3/d); (2D matrix: num_grps X num_t)
ConsumptionRatesC               = reshape(ConsumptionRatesC, [num_grps, 1, num_t]); % consumption rate of each consumer c group; (mmoles N/m3/d); (3D matrix: num_grps X 1 X num_t)
ConsumptionRatesP               = reshape(ConsumptionRatesC, [1, num_grps, num_t]); % consumption rate of each producer p group; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_t)

ConsumptionRatesC_repmat        = repmat(ConsumptionRatesC, [1, num_grps, 1]); % replicate consumer c consumption across clms; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_t) NOTE: VERTICAL vectors replicated across clms;
ConsumptionRatesP_repmat        = repmat(ConsumptionRatesP, [num_grps, 1, 1]);	% replicate producer p consumption down rows; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: HORIZONTAL vectors replicated down rows

consumptionC_initial_repmat     = repmat(consumptionC_initial_repmat, [1, 1, num_t]); % initial consumption inflow rate; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_t); NOTE: replicated vertical vectors across columns
% -------------------------------------------------------------------------

% step 4c: calculate FunctionalResponse consumption scaler @ each t -------
FunctionalResponse                              = ((1+FunctionalResponseParams_repmat) .* ConsumptionRatesC_repmat) ./ (FunctionalResponseParams_repmat .* consumptionC_initial_repmat + ConsumptionRatesC_repmat); % (matrix aligned with EnergyBudget); (unitless); (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_t)
FunctionalResponse(isnan(FunctionalResponse))   = 1; % catch NaNs and set to FunctionalResponse = 1 (caused by 0 biomass cells in sub-surface boxes and possibly by dummy variable NaNs) QQQ try to remove this step
% -------------------------------------------------------------------------

% step 4d: reshape and repmat q_TemperatureScaler -------------------------
%          Thornton-Lessem temperature adjustment to consumption rate; value between 0 and 1; (2D matrix: num_t X num_grps)
%          scaler gets applied to consumers (rows)
%          reshape to 3D matrix (num_grps X 1 X num_t)
%          replicate across producers (columns)
q_TemperatureScaler      	= q_TemperatureScaler'; % transpose; (2D matrix: num_grps X num_t)
q_TemperatureScaler        	= reshape(q_TemperatureScaler, [num_grps, 1, num_t]); % (3D matrix: num_grps X 1 X num_t)
q_TemperatureScaler_repmat	= repmat(q_TemperatureScaler, [1, num_grps, 1]); % (3D matrix: num_grps X num_grps X num_t)
% -------------------------------------------------------------------------

% step 4e: calculate CONSUMPTION matrix Q_cp @ each t ---------------------
Q_cp                    	= EnergyBudget_repmat .* FunctionalResponse .* q_TemperatureScaler_repmat .* ConsumptionRatesP_repmat; % (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_t)
% *************************************************************************


% end m-file***************************************************************