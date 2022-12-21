function [FluxSummary, FluxTimeseries, ODEinput] = f_ECOTRANdynamic_FluxSummary_NCC_09082022(read_file, target_year, target_qtr, target_box)
% Read the output of an individual ECOTRANdynamic run and return the mean
% biomass (q), physical fluxes, and Q_cp of a specified period of time.
%
% calls:
%       f_calcQcp_result_09072022
% 
% takes:
%       read_file
%       target_year     year(s) to use for analyses (vector)
%       target_qtr      month(s) to use for analyses (vector)
%       target_box      sub-region boxes to use for analyses (vector)
%
% returns:
%   FluxSummary
%       mean_biomass_daily           	biomass average of all hours across all days in the target_year(s) & target_qtr(s)                  (mmoles N/m3)       (2D matrix: num_grps X [epi meso bathy benthic epi+meso epi+meso+bathy total])
%       mean_biomass_midnight         	midnight biomasses average across all days in the target_year(s) & target_qtr(s)                    (mmoles N/m3)       (2D matrix: num_grps X [epi meso bathy benthic epi+meso epi+meso+bathy total])
%       mean_biomass_noon             	noon biomasses average across all days in the target_year(s) & target_qtr(s)                        (mmoles N/m3)       (2D matrix: num_grps X [epi meso bathy benthic epi+meso epi+meso+bathy total])
%
%       mean_production_daily           production average of all hours across all days in the target_year(s) & target_qtr(s)           	(mmoles N/m3/d)	(2D matrix: num_grps X [epi meso bathy benthic epi+meso epi+meso+bathy total])
%       mean_production_midnight     	midnight production average across all days in the target_year(s) & target_qtr(s)                 	(mmoles N/m3/d)	(2D matrix: num_grps X [epi meso bathy benthic epi+meso epi+meso+bathy total])
%       mean_production_noon            noon production average across all days in the target_year(s) & target_qtr(s)                     	(mmoles N/m3/d)	(2D matrix: num_grps X [epi meso bathy benthic epi+meso epi+meso+bathy total])
%
%       mean_ConsumptionRates_daily    	q rate average of all hours across all days in the target_year(s) & target_qtr(s)                   (mmoles N/m3/d)	(2D matrix: num_grps X num_boxes)
%       mean_ConsumptionRates_midnight	midnight q rate average across all days in the target_year(s) & target_qtr(s)                       (mmoles N/m3/d)	(2D matrix: num_grps X num_boxes)
%       mean_ConsumptionRates_noon    	noon q rate average across all days in the target_year(s) & target_qtr(s)                           (mmoles N/m3/d)	(2D matrix: num_grps X num_boxes)
%
%       mean_Q_cp_daily              	Q_cp matrix average of all hours across all days in the target_year(s) & target_qtr(s)              (mmoles N/m3/d)	(3D matrix: num_grps X num_grps X num_boxes)
%       mean_Q_cp_midnight            	midnight Q_cp matrix averaged across all days in the target_year(s) & target_qtr(s)                 (mmoles N/m3/d)	(3D matrix: num_grps X num_grps X num_boxes)
%       mean_Q_cp_noon                	noon Q_cp matrix averaged across all days in the target_year(s) & target_qtr(s)                     (mmoles N/m3/d)	(3D matrix: num_grps X num_grps X num_boxes)
%       mean_Q_cp_subdaily              subdaily Q_cp matrix averaged across all days in the target_year(s) & target_qtr(s)                 (mmoles N/m3/d) (4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
%
%       mean_VolumeFlux_advection_daily     volume flux average of all hours across all days in the target_year(s) & target_qtr(s) (m3/d); (horizontal vector: 1 X num_fluxes_advection)
%       mean_VolumeFlux_HorizMixing_daily	volume flux average of all hours across all days in the target_year(s) & target_qtr(s) (m3/d); (horizontal vector: 1 X num_fluxes_HorizontalMixing)
%       mean_VolumeFlux_VertMixing_daily	volume flux average of all hours across all days in the target_year(s) & target_qtr(s) (m3/d); (horizontal vector: 1 X num_fluxes_VerticalMixing)
% 
%       mean_ImportDriver_advection_daily	biomass flux average of all hours across all days in the target_year(s) & target_qtr(s) (t N/d); (2D matrix: num_grps X num_boxes DESTINATION)
%       mean_ImportDriver_HorizMixing_daily biomass flux average of all hours across all days in the target_year(s) & target_qtr(s) (t N/d); (2D matrix: num_grps X num_boxes DESTINATION)
%       mean_ImportDriver_migration_daily   biomass flux average of all hours across all days in the target_year(s) & target_qtr(s) (t N/d); (2D matrix: num_grps X num_boxes DESTINATION)
%
%       mean_advection_daily            average advecton flux over all hours across all days in the target_year(s) & target_qtr(s)          (t N/d)      	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_HorizMixing_daily          average horizontal mixing flux over all hours across all days in the target_year(s) & target_qtr(s)	(t N/d)     	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_VertMixing_daily           average vertical mixing flux over all hours across all days in the target_year(s) & target_qtr(s)	(t N/d)       	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_sinking_daily              average sinking flux over all hours across all days in the target_year(s) & target_qtr(s)           (t N/d)       	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_migration_daily            average migration flux over all hours across all days in the target_year(s) & target_qtr(s)         (t N/d)        	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
% 
%       mean_advection_midnight         average advecton flux at midnight across all days in the target_year(s) & target_qtr(s)             (t N/d)        	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_HorizMixing_midnight       average horizontal mixing flux at midnight across all days in the target_year(s) & target_qtr(s)    (t N/d)       	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_VertMixing_midnight        average vertical mixing flux at midnight across all days in the target_year(s) & target_qtr(s)      (t N/d)      	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_sinking_midnight           average sinking flux at midnight across all days in the target_year(s) & target_qtr(s)              (t N/d)        	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_migration_midnight         average migration flux at midnight across all days in the target_year(s) & target_qtr(s)            (t N/d)        	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
% 
%       mean_advection_noon             average advecton flux at noon across all days in the target_year(s) & target_qtr(s)                 (t N/d)       	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_HorizMixing_noon           average horizontal mixing flux at noon across all days in the target_year(s) & target_qtr(s)        (t N/d)        	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_VertMixing_noon            average vertical mixing flux at noon across all days in the target_year(s) & target_qtr(s)          (t N/d)       	(2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_sinking_noon               average sinking flux at noon across all days in the target_year(s) & target_qtr(s)                  (t N/d)         (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%       mean_migration_noon             average migration flux at noon across all days in the target_year(s) & target_qtr(s)                (t N/d)         (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
%
% revision date: 9-8-2022
%	9/6/2022 cleaning up code and checking units
%	9/8/2022 updated version of f_calcQcp_result_09072022



% *************************************************************************
% STEP 1: define operating conditions--------------------------------------

num_targetBoxes = length(target_box);


% step 1a: define run files -----------------------------------------------
% filedirectory	= '/Users/jimsebi/Documents/12_GoMex_project/5_calibrations_NutrientInput/';
% run_file        = 'GoMexOcn_calA3_05-Jan-2022.mat';
% -------------------------------------------------------------------------


% step 1b: define time period for analysis --------------------------------
% target_year             = 2005; % year(s) to use for analyses
% target_qtr              = [12]; % month(s) to use for analyses
% -------------------------------------------------------------------------
% *************************************************************************





% *************************************************************************
% STEP 2: load run results file and compile mean qtr biomasses ------------

% step 2a: load in physical fluxes ----------------------------------------
% load(read_file, 'store_T', 'PHYSICSinput', 'ECOTRANphysics', 'ECOTRANmigration', 'ODEinput'); % load model run
load(read_file, 'store_T', 'PHYSICSinput', 'ECOTRANphysics', 'ODEinput'); % load model run


% get time vectors, num_grps, num_boxes, q/b, physical fluxes
datestart                   = PHYSICSinput.datestart;
dt                          = PHYSICSinput.dt;
T                           = store_T(:,1,1);
num_t                       = length(T);
matlabdate                  = datevec(T + datestart - 1); % model dates as time vector; (2D matrix: num_t X [year month day hour minute second])
looky_driver            	= ODEinput.looky_driver;	% row address(es) of driver group(s) (e.g., NO3)
num_grps                 	= ODEinput.num_grps;
num_boxes                  	= ODEinput.num_boxes;       % number of spatial boxes in model
grp_row                     = 1:num_grps;

num_subdaily                = 24/(dt*24); % num_subdaily will be redefined again below (but the value should not change)

BoxVolume                   = ODEinput.BoxVolume;                           % (m3); (2D matrix: num_t X num_boxes)
BoxVolume                   = reshape(BoxVolume, [num_t, 1, num_boxes]);	% (m3); (2D matrix: num_t X 1 X num_boxes)
BoxVolume                   = repmat(BoxVolume, [1, num_grps, 1]);          % (m3); (3D matrix: num_t X num_grps X num_boxes)

ADVECTION_compact         	= ODEinput.ADVECTION_compact;           % (m3/d); (2D matrix: num_t X num_fluxes_advection)
HORIZONTALMIXING_compact	= ODEinput.HORIZONTALMIXING_compact;	% (m3/d); (2D matrix: num_t X num_fluxes_HorizontalMixing)
VERTICALMIXING_compact   	= ODEinput.VERTICALMIXING_compact;      % (m3/d); (2D matrix: num_t X num_fluxes_VerticalMixing)
SINKING_compact           	= ODEinput.SINKING_compact;             % (m3/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
MIGRATION_compact       	= ODEinput.MIGRATION_compact;           % (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)

num_fluxes_advection        = ODEinput.num_fluxes_advection;
num_fluxes_HorizontalMixing = ODEinput.num_fluxes_HorizontalMixing;
num_fluxes_VerticalMixing	= ODEinput.num_fluxes_VerticalMixing;
num_fluxes_sinking          = ODEinput.num_fluxes_sinking;
num_fluxes_migration        = ODEinput.num_fluxes_migration;
% -------------------------------------------------------------------------


% step 2b: calculate biological fluxes for 1 model & 1 sub-region ---------
for box_loop = 1:num_targetBoxes        % loop through each box
    
    currentBox              = target_box(box_loop);
    
    [re_Y_currentBox, Q_cp_currentBox, qb_currentBox, pb_currentBox]	= f_calcQcp_result_09072022(read_file, currentBox);      % CONSUMPTION matrix (Q_cp); (mmoles N/m3/d); (3D matrix: num_grps c X num_grps p X num_t)

    re_Y(:, :, box_loop)	= re_Y_currentBox;	% (3D matrix: num_t X num_grps X num_targetBoxes)
    qb(:, :, box_loop)      = qb_currentBox;	% (3D matrix: num_t X num_grps X num_targetBoxes)
    pb(:, :, box_loop)      = pb_currentBox;	% (3D matrix: num_t X num_grps X num_targetBoxes)
    Q_cp(:, :, :, box_loop)	= Q_cp_currentBox;	% (4D matrix: num_grps X num_grps X num_t X num_targetBoxes)
    
end % (box_loop)
% -------------------------------------------------------------------------


% step 2c: convert from consumption rate q (mmoles N/m3/d) to biomass density (mmoles N/m3) and to production rate (mmoles N/m3/d)
biomass_total               = re_Y ./ qb;           % use q/b value to convert q rates to standing stock biomass density; (mmole N/m3); (3D matrix: num_t X num_grps X num_targetBoxes)
production_total            = biomass_total .* pb;	% use p/b value to convert standing stock biomass density to production rates; (mmole N/m3/d); (3D matrix: num_t X num_grps X num_targetBoxes)
% *************************************************************************





% *************************************************************************
% STEP 3: select period of interest----------------------------------------
looky_year                          = find(ismember(matlabdate(:, 1), target_year));
looky_month                         = find(ismember(matlabdate(:, 2), target_qtr));
looky_period                        = find(ismember(looky_year, looky_month));
target_period                       = looky_year(looky_period);
matlabdate_select                   = matlabdate(target_period, :);

unique_years                        = unique(matlabdate_select(:, 1));
unique_months                       = unique(matlabdate_select(:, 2));
unique_dates                        = unique(matlabdate_select(:, 1:3), 'rows');

num_years                           = length(unique_years);
num_months                          = length(unique_months);
num_dates                       	= length(unique_dates); 

biomass_total_select                = biomass_total(target_period, :, :);           % (mmole N/m3);	  (3D matrix: length(target_period) X num_grps X num_targetBoxes)
production_total_select           	= production_total(target_period, :, :);     	% (mmole N/m3/d); (3D matrix: length(target_period) X num_grps X num_targetBoxes)
ConsumptionRates_select             = re_Y(target_period, :, :);                    % (mmole N/m3/d); (3D matrix: length(target_period) X num_grps X num_targetBoxes)
BoxVolume_select                    = BoxVolume(target_period, :, :);               % (m3);           (3D matrix: length(target_period) X num_grps X num_boxes)
Q_cp_select                         = Q_cp(:, :, target_period, :);                 % (mmole N/m3/d); (4D matrix: num_grps X num_grps X length(target_period) X num_targetBoxes)

ADVECTION_select                    = ADVECTION_compact(target_period, :);          % (m3/d);         (2D matrix: length(target_period) X num_fluxes_advection)
HORIZONTALMIXING_select             = HORIZONTALMIXING_compact(target_period, :);	% (m3/d);         (2D matrix: length(target_period) X num_fluxes_HorizontalMixing)
VERTICALMIXING_select               = VERTICALMIXING_compact(target_period, :);     % (m3/d);         (2D matrix: length(target_period) X num_fluxes_VerticalMixing)
SINKING_select                      = SINKING_compact(target_period, :, :);         % (m3/d);         (3D matrix: length(target_period) X num_grps X num_fluxes_sinking)
MIGRATION_select                    = MIGRATION_compact(target_period, :, :);       % (m3/d);         (3D matrix: length(target_period) X num_grps X num_fluxes_migration)
% *************************************************************************





% *************************************************************************
% STEP 4: harvest results from the period of interest----------------------
% step 4a: initialize variables -------------------------------------------
counter                                     = 0; % initialize

intraODEinput.flux_domain_import_t          = ODEinput.flux_domain_import_t;            % initialized variable, all zeros; (2D matrix: num_grps X num_boxes+1)
intraODEinput.flux_domain_export_t          = ODEinput.flux_domain_export_t;            % initialized variable, all zeros; (2D matrix: num_grps X num_boxes+1)
intraODEinput.flux_domain_import_driver_t	= ODEinput.flux_domain_import_driver_t;     % initialized variable, all zeros; (2D matrix: num_grps X num_boxes DESTINY)

build_biomass_dailymean                     = zeros(num_dates, num_grps, num_boxes);            % initialize
build_biomass_midnight                      = zeros(num_dates, num_grps, num_boxes);            % initialize
build_biomass_noon                          = zeros(num_dates, num_grps, num_boxes);            % initialize

build_production_dailymean               	= zeros(num_dates, num_grps, num_boxes);            % initialize
build_production_midnight                 	= zeros(num_dates, num_grps, num_boxes);            % initialize
build_production_noon                     	= zeros(num_dates, num_grps, num_boxes);            % initialize

build_ConsumptionRates_dailymean            = zeros(num_dates, num_grps, num_boxes);    % initialize; (3D matrix: num_dates X num_grps X num_targetBoxes)
build_ConsumptionRates_midnight             = zeros(num_dates, num_grps, num_boxes);    % initialize; (3D matrix: num_dates X num_grps X num_targetBoxes)
build_ConsumptionRates_noon                 = zeros(num_dates, num_grps, num_boxes);    % initialize; (3D matrix: num_dates X num_grps X num_targetBoxes)

build_Q_cp_dailymean                        = zeros(num_grps, num_grps, num_dates, num_boxes);     % initialize; (4D matrix: num_grps X num_grps X num_dates X num_targetBoxes)
build_Q_cp_midnight                         = zeros(num_grps, num_grps, num_dates, num_boxes);     % initialize; (4D matrix: num_grps X num_grps X num_dates X num_targetBoxes)
build_Q_cp_noon                             = zeros(num_grps, num_grps, num_dates, num_boxes);     % initialize; (4D matrix: num_grps X num_grps X num_dates X num_targetBoxes)
build_Q_cp_subdaily                      	= zeros(num_grps, num_grps, num_subdaily, num_boxes);     % QQQ not sure time dimension correct here; initialize; (4D matrix: num_grps X num_grps X num_dates X num_targetBoxes)

build_VolumeFlux_advection_dailymean        = zeros(num_dates, num_fluxes_advection);           % (m3/d);	(3D matrix: num_dates X num_fluxes_advection)
build_VolumeFlux_HorizMixing_dailymean      = zeros(num_dates, num_fluxes_HorizontalMixing);	% (m3/d);	(3D matrix: num_dates X num_fluxes_HorizontalMixing)
build_VolumeFlux_VertMixing_dailymean       = zeros(num_dates, num_fluxes_VerticalMixing);      % (m3/d);	(3D matrix: num_dates X num_fluxes_VerticalMixing)

build_ImportDriver_advection_dailymean      = zeros(num_dates, num_grps, num_boxes);	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
build_ImportDriver_HorizMixing_dailymean	= zeros(num_dates, num_grps, num_boxes);	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
build_ImportDriver_migration_dailymean      = zeros(num_dates, num_grps, num_boxes);	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)

build_advection_dailymean                   = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_HorizMixing_dailymean                 = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_VertMixing_dailymean                  = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_sinking_dailymean                     = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_migration_dailymean                   = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)

build_advection_midnight                    = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_HorizMixing_midnight                  = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_VertMixing_midnight                   = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_sinking_midnight                      = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_migration_midnight                    = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)

build_advection_noon                        = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_HorizMixing_noon                      = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_VertMixing_noon                       = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_sinking_noon                          = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
build_migration_noon                        = zeros(num_dates, num_grps, num_boxes);    % initialize; (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
% -------------------------------------------------------------------------


% step 3b: harvest data from selected year(s) -----------------------------
for year_loop = 1:num_years
    current_year                    = unique_years(year_loop);
    looky_year                      = find(matlabdate_select(:, 1) == current_year);
    current_matlabdate_year         = matlabdate_select(looky_year, :);
    
    current_biomass_year            = biomass_total_select(looky_year, :, :);       % (mmoles N/m3);   (3D matrix: length(looky_year) X num_grps X num_targetBoxes)
    current_production_year       	= production_total_select(looky_year, :, :);	% (mmoles N/m3/d); (3D matrix: length(looky_year) X num_grps X num_targetBoxes)
    current_ConsumptionRates_year	= ConsumptionRates_select(looky_year, :, :);    % (mmoles N/m3/d); (3D matrix: length(looky_year) X num_grps X num_targetBoxes)
	current_BoxVolume_year          = BoxVolume_select(looky_year, :, :);           % (m3);            (3D matrix: length(looky_year) X num_grps X num_boxes)
	current_Q_cp_year               = Q_cp_select(:, :, looky_year, :);          	% (mmoles N/m3/d); (4D matrix: num_grps X num_grps X length(looky_year) X num_targetBoxes)

    current_ADVECTION_year          = ADVECTION_select(looky_year, :);              % (m3/d);          (2D matrix: length(looky_year) X num_fluxes_advection)
    current_HORIZONTALMIXING_year	= HORIZONTALMIXING_select(looky_year, :);       % (m3/d);          (2D matrix: length(looky_year) X num_fluxes_HorizontalMixing)
    current_VERTICALMIXING_year     = VERTICALMIXING_select(looky_year, :);         % (m3/d);          (2D matrix: length(looky_year) X num_fluxes_VerticalMixing)
    current_SINKING_year            = SINKING_select(looky_year, :, :);             % (m3/d);          (3D matrix: length(looky_year) X num_grps X num_fluxes_sinking)
    current_MIGRATION_year          = MIGRATION_select(looky_year, :, :);           % (m3/d);          (3D matrix: length(looky_year) X num_grps X num_fluxes_migration)
    % ---------------------------------------------------------------------
    
    
    % step 3c: harvest results from selected month(s) ---------------------
    for month_loop = 1:num_months
        current_month                   = unique_months(month_loop);
        looky_month                     = find(current_matlabdate_year(:, 2) == current_month);
        current_matlabdate_month        = current_matlabdate_year(looky_month, :);
        
        unique_days                 	= unique(current_matlabdate_month(:, 3));
        num_days                    	= length(unique_days);
        
        current_biomass_month           = current_biomass_year(looky_month, :, :);              % (mmoles N/m3);   (3D matrix: length(looky_month) X num_grps X num_targetBoxes)
        current_production_month      	= current_production_year(looky_month, :, :);        	% (mmoles N/m3/d); (3D matrix: length(looky_month) X num_grps X num_targetBoxes)
        current_ConsumptionRates_month	= current_ConsumptionRates_year(looky_month, :, :);     % (mmoles N/m3/d); (3D matrix: length(looky_month) X num_grps X num_targetBoxes)
        current_BoxVolume_month         = current_BoxVolume_year(looky_month, :, :);            % (m3);            (3D matrix: length(looky_month) X num_grps X num_boxes)
        current_Q_cp_month              = current_Q_cp_year(:, :, looky_month, :);              % (mmoles N/m3/d); (4D matrix: num_grps X num_grps X length(looky_month) X num_targetBoxes)

        current_ADVECTION_month         = current_ADVECTION_year(looky_month, :);               % (m3/d);          (2D matrix: length(looky_month) X num_fluxes_advection)
        current_HORIZONTALMIXING_month	= current_HORIZONTALMIXING_year(looky_month, :);        % (m3/d);          (2D matrix: length(looky_month) X num_fluxes_HorizontalMixing)
        current_VERTICALMIXING_month	= current_VERTICALMIXING_year(looky_month, :);          % (m3/d);          (2D matrix: length(looky_month) X num_fluxes_VerticalMixing)
        current_SINKING_month           = current_SINKING_year(looky_month, :, :);              % (m3/d);          (3D matrix: length(looky_month) X num_grps X num_fluxes_sinking)
        current_MIGRATION_month         = current_MIGRATION_year(looky_month, :, :);            % (m3/d);          (3D matrix: length(looky_month) X num_grps X num_fluxes_migration)
        % -----------------------------------------------------------------
        
        
        % step 3d: harvest results from each day of selected month(s) -----
        for day_loop = 1:num_days
            current_day                     = unique_days(day_loop);
            looky_day                       = find(current_matlabdate_month(:, 3) == current_day);
            current_matlabdate_day          = current_matlabdate_month(looky_day, :);
            num_subdaily                    = length(looky_day);
            
            counter                         = counter + 1; % increment matrix build counter by 1
            
            current_biomass_day          	= current_biomass_month(looky_day,          :, :);        	% (mmoles N/m3);	(3D matrix: num_subdaily X num_grps X num_targetBoxes)
            current_production_day       	= current_production_month(looky_day,       :, :);     	% (mmoles N/m3/d);	(3D matrix: num_subdaily X num_grps X num_targetBoxes)
            current_ConsumptionRates_day 	= current_ConsumptionRates_month(looky_day, :, :);	% (mmoles N/m3/d);  (3D matrix: num_subdaily X num_grps X num_targetBoxes)
            current_BoxVolume_day        	= current_BoxVolume_month(looky_day,        :, :);	% (m3);             (3D matrix: num_subdaily X num_grps X num_boxes)
            current_Q_cp_day                = current_Q_cp_month(:, :, looky_day, :);           % (mmoles N/m3/d);  (4D matrix: num_grps X num_grps X length(looky_day) X num_targetBoxes)

            current_ADVECTION_day           = current_ADVECTION_month(looky_day,        :);            % (m3/d);           (2D matrix: length(looky_day) X num_fluxes_advection)
            current_HORIZONTALMIXING_day	= current_HORIZONTALMIXING_month(looky_day, :);     % (m3/d);           (2D matrix: length(looky_day) X num_fluxes_HorizontalMixing)
            current_VERTICALMIXING_day      = current_VERTICALMIXING_month(looky_day,   :);     	% (m3/d);           (2D matrix: length(looky_day) X num_fluxes_VerticalMixing)
            current_SINKING_day             = current_SINKING_month(looky_day,          :, :);        	% (m3/d);           (3D matrix: length(looky_day) X num_grps X num_fluxes_sinking)
            current_MIGRATION_day           = current_MIGRATION_month(looky_day,        :, :);       	% (m3/d);           (3D matrix: length(looky_day) X num_grps X num_fluxes_migration)
            % -------------------------------------------------------------
            
            
            % step 3e: process sub-daily results --------------------------
            if num_subdaily > 0
                
                % initialize subdaily physical flux terms -----------------
                build_VolumeFlux_advection_subdaily         = zeros(num_subdaily, num_fluxes_advection);        % (m3/d); (2D matrix: num_subdaily X num_fluxes_advection)
                build_VolumeFlux_HorizMixing_subdaily       = zeros(num_subdaily, num_fluxes_HorizontalMixing);	% (m3/d); (2D matrix: num_subdaily X num_fluxes_HorizontalMixing)
                build_VolumeFlux_VertMixing_subdaily        = zeros(num_subdaily, num_fluxes_VerticalMixing);	% (m3/d); (2D matrix: num_subdaily X num_fluxes_VerticalMixing)
                
                build_ImportDriver_advection_subdaily       = zeros(num_subdaily, num_grps, num_boxes);         % (mmoles N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
               	build_ImportDriver_HorizMixing_subdaily   	= zeros(num_subdaily, num_grps, num_boxes);         % (mmoles N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
               	build_ImportDriver_migration_subdaily       = zeros(num_subdaily, num_grps, num_boxes);         % (mmoles N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION); NOTE; that there can be a migration input driver but not a vertical mixing input driver
                
                build_flux_advection_subdaily               = zeros(num_subdaily, num_grps, num_boxes); % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_HorizMixing_subdaily             = zeros(num_subdaily, num_grps, num_boxes); % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_VertMixing_subdaily              = zeros(num_subdaily, num_grps, num_boxes); % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_sinking_subdaily                 = zeros(num_subdaily, num_grps, num_boxes); % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_migration_subdaily               = zeros(num_subdaily, num_grps, num_boxes); % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                for subdaily_loop = 1:num_subdaily
                    
                    BoxVolume_t                     = squeeze(current_BoxVolume_day(subdaily_loop, :, :));         % (m3);           (2D matrix: num_grps X num_boxes); NOTE squeeze
                    intraODEinput.repmat_BoxVolume	= BoxVolume_t;                                                 % (m3);           (2D matrix: num_grps X num_boxes)
                    
                    biomass_t                       = squeeze(current_biomass_day(subdaily_loop, :, 1:num_targetBoxes)); % (mmoles N/m3);  (2D matrix: num_grps X num_targetBoxes); NOTE squeeze
                    biomass_plus1                	= reshape(biomass_t, [num_grps, 1, num_targetBoxes]);                % (mmoles N/m3);  (3D matrix: num_grps X 1 X num_targetBoxes)
                    biomass_plus1(:, 1, (end+1))	= 0;                % add 1 layer for boundary fluxes;           (mmoles N/m3);  (3D matrix: num_grps X 1 X num_boxes+1)
                    biomass_boundary_t           	= biomass_plus1;	%                                            (mmoles N/m3);  (3D matrix: num_grps X 1 X num_boxes+1)
                    biomass_MigratorBoundary_t      = biomass_plus1;	% default case is to use same boundary conditions for migration as for physical fluxes; NOTE: deactivate this line for special definition of boundary biomass for migrators (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    
                    ConsumptionRates_t            	= squeeze(current_ConsumptionRates_day(subdaily_loop, :, :));  % (mmoles N/m3/d); (2D matrix: num_grps X num_boxes)
                    ConsumptionRates_t          	= reshape(ConsumptionRates_t, [num_grps, 1, num_boxes]);       % (mmoles N/m3/d); (3D matrix: num_grps X 1 X num_boxes)

                    ConsumptionRates_t_transpose    = current_ConsumptionRates_day(subdaily_loop, :, :); % reshape to 3D matrix ; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes); 
                    % -----------------------------------------------------
                    
                    
                    % Process physical fluxes -----------------------------
                    % ADVECTION -------------------------------------------
                    ADVECTION_compact_t                     = current_ADVECTION_day(subdaily_loop, :);	% (m3/d); (horizontal vector: 1 X num_fluxes_advection)
                    intraODEinput.compact_flux_t            = repmat(ADVECTION_compact_t, [num_grps, 1]); % replicate each flux term for all model groups; (m3/d); (2D matrix: num_grps X num_fluxes_advection)
                    intraODEinput.biomass_plus1             = biomass_plus1;                            % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.biomass_boundary          = biomass_boundary_t;                      	% boundary biomass conditions @ t; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler;                 % (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINY);
                    intraODEinput.looky_flux                = ODEinput.looky_AdvectionFlux;             % (2D matrix: num_fluxes_advection X 5-->[(DESTINY box) (SOURCE box) (DESTINY box address in compact_flux) (SOURCE box address in compact_flux) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
                    intraODEinput.looky_boundary_import     = ODEinput.looky_AdvectionBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
                    intraODEinput.looky_boundary_export     = ODEinput.looky_AdvectionBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
                    intraODEinput.looky_driver           	= looky_driver;                             % row address(es) of driver group(s) (e.g., NO3)
                    intraODEinput.repmat_looky_source          = repmat(ODEinput.looky_AdvectionFlux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_looky_destiny         = repmat(ODEinput.looky_AdvectionFlux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_grp_row          	= repmat(grp_row', [1, num_fluxes_advection]); % addressses; (2D matrix: num_grps X num_fluxes)
                    Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    %	Fluxes_t.flux_import_t	-->flux INTO each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does NOT include domain import (see flux_domain_import_t for domain input rates)
                    %	Fluxes_t.flux_export_t	-->flux OUT OF each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does include domain export
                    %	Fluxes_t.flux_net_t     -->net biomass flux into(+) or out of(-) each box; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: does include domain import & export
                    flux_net_advection_t                    = Fluxes_t.flux_net_t;                      % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    flux_domain_ImportDriver_advection_t	= Fluxes_t.flux_domain_import_driver_t;     % net biomass flux OF DRIVER GRP across domain boundary INTO each box (all other grps = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION)
                    % -----------------------------------------------------
                
                    % HORIZONTAL MIXING -----------------------------------
                    HORIZONTALMIXING_compact_t           	= current_HORIZONTALMIXING_day(subdaily_loop, :); % (m3/d); (horizontal vector: 1 X num_fluxes_HorizontalMixing)
                    intraODEinput.compact_flux_t            = repmat(HORIZONTALMIXING_compact_t, [num_grps, 1]); % replicate each flux term for all model groups; (m3/d); (2D matrix: num_grps X num_fluxes_HorizontalMixing)
                    intraODEinput.biomass_plus1             = biomass_plus1;                                    % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.biomass_boundary          = biomass_boundary_t;                              	% boundary biomass conditions @ t; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler;                         % (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
                    intraODEinput.looky_flux               	= ODEinput.looky_HorizontalMixingFlux;              % (2D matrix: num_fluxes_HorizontalMixing X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
                    intraODEinput.looky_boundary_import    	= ODEinput.looky_HorizontalMixingBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
                    intraODEinput.looky_boundary_export    	= ODEinput.looky_HorizontalMixingBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
                    intraODEinput.looky_driver           	= looky_driver;                                     % row address(es) of driver group(s) (e.g., NO3)
                    intraODEinput.repmat_looky_source       = repmat(ODEinput.looky_HorizontalMixingFlux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_looky_destiny      = repmat(ODEinput.looky_HorizontalMixingFlux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_grp_row          	= repmat(grp_row', [1, num_fluxes_HorizontalMixing]); % addressses; (2D matrix: num_grps X num_fluxes)
                    Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    %	Fluxes_t.flux_import_t	-->flux INTO each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does NOT include domain import (see flux_domain_import_t for domain input rates)
                    %	Fluxes_t.flux_export_t	-->flux OUT OF each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does include domain export
                    %	Fluxes_t.flux_net_t     -->net biomass flux into(+) or out of(-) each box; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: does include domain import & export
                    flux_net_HorizMixing_t                  = Fluxes_t.flux_net_t;                              % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    flux_domain_ImportDriver_HorizMixing_t	= Fluxes_t.flux_domain_import_driver_t;             % net biomass flux OF DRIVER GRP across domain boundary INTO each box (all other grps = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION)
                    % -----------------------------------------------------

                    % VERTICAL MIXING -------------------------------------
                    VERTICALMIXING_compact_t                = current_VERTICALMIXING_day(subdaily_loop, :); % (m3/d); (horizontal vector: 1 X num_fluxes_VerticalMixing)
                    intraODEinput.compact_flux_t            = repmat(VERTICALMIXING_compact_t, [num_grps, 1]); % replicate each flux term for all model groups; (m3/d); (2D matrix: num_grps X num_fluxes_VerticalMixing)
                    intraODEinput.biomass_plus1             = biomass_plus1;                                % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.biomass_boundary          = biomass_boundary_t;                          	% boundary biomass conditions @ t; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler;                     % (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
                    intraODEinput.looky_flux              	= ODEinput.looky_VerticalMixingFlux;           	% (2D matrix: num_fluxes_VerticalMixing X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
                    intraODEinput.looky_boundary_import   	= ODEinput.looky_VerticalMixingBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
                    intraODEinput.looky_boundary_export  	= ODEinput.looky_VerticalMixingBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
                    intraODEinput.looky_driver           	= looky_driver;                                 % row address(es) of driver group(s) (e.g., NO3)
                    intraODEinput.repmat_looky_source       = repmat(ODEinput.looky_VerticalMixingFlux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_looky_destiny      = repmat(ODEinput.looky_VerticalMixingFlux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_grp_row          	= repmat(grp_row', [1, num_fluxes_VerticalMixing]);         % addressses; (2D matrix: num_grps X num_fluxes)                    
                    Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    %	Fluxes_t.flux_import_t	-->flux INTO each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does NOT include domain import (see flux_domain_import_t for domain input rates)
                    %	Fluxes_t.flux_export_t	-->flux OUT OF each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does include domain export
                    %	Fluxes_t.flux_net_t     -->net biomass flux into(+) or out of(-) each box; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: does include domain import & export
                    flux_net_VertMixing_t                   = Fluxes_t.flux_net_t;                          % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    % -----------------------------------------------------

                    % SINKING ---------------------------------------------
                    SINKING_compact_t                       = squeeze(current_SINKING_day(subdaily_loop, :, :)); % (m3/d); (2D matrix: num_grps X num_fluxes_sinking); NOTE squeeze
                    intraODEinput.compact_flux_t            = SINKING_compact_t;                        % (m3/d); (2D matrix: num_grps X num_fluxes_sinking)
                    intraODEinput.biomass_plus1             = biomass_plus1;                            % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.biomass_boundary          = biomass_boundary_t * 0;                 	% NOTE: no boundary conditions for sinking (all grps = 0); boundary biomass conditions @ t; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler * 0;             % NOTE: apply no retention for sinking; (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
                    intraODEinput.looky_flux              	= ODEinput.looky_SinkingFlux;               % (2D matrix: num_fluxes_sinking X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
                    intraODEinput.looky_boundary_import    	= ODEinput.looky_SinkingBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
                    intraODEinput.looky_boundary_export    	= ODEinput.looky_SinkingBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
                    intraODEinput.looky_driver           	= looky_driver;                             % row address(es) of driver group(s) (e.g., NO3)
                    intraODEinput.repmat_looky_source       = repmat(ODEinput.looky_SinkingFlux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_looky_destiny      = repmat(ODEinput.looky_SinkingFlux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_grp_row          	= repmat(grp_row', [1, num_fluxes_sinking]);         % addressses; (2D matrix: num_grps X num_fluxes)                    
                    Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    %	Fluxes_t.flux_import_t	-->flux INTO each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does NOT include domain import (see flux_domain_import_t for domain input rates)
                    %	Fluxes_t.flux_export_t	-->flux OUT OF each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does include domain export
                    %	Fluxes_t.flux_net_t     -->net biomass flux into(+) or out of(-) each box; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: does include domain import & export
                    flux_net_sinking_t                      = Fluxes_t.flux_net_t;                      % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    % -----------------------------------------------------

                    % MIGRATION -------------------------------------------
                    MIGRATION_compact_t                   	= squeeze(current_MIGRATION_day(subdaily_loop, :, :)); % (m3/d); (2D matrix: num_grps X num_fluxes_migration); NOTE squeeze
                    intraODEinput.compact_flux_t            = MIGRATION_compact_t;                      % (m3/d); (2D matrix: num_grps X num_fluxes_migration)
                    intraODEinput.biomass_plus1             = biomass_plus1;                            % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.biomass_boundary          = biomass_MigratorBoundary_t;         	    % boundary biomass conditions @ t for migrators; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1)
                    intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler * 0;             % (scaler = 0 for migration); (2D matrix: num_grps X num_boxes DESTINATION);
                    intraODEinput.looky_flux                = ODEinput.looky_MigrationFlux;             % (2D matrix: num_fluxes_migration X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
                    intraODEinput.looky_boundary_import     = ODEinput.looky_MigrationBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
                    intraODEinput.looky_boundary_export     = ODEinput.looky_MigrationBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
                    intraODEinput.looky_driver           	= looky_driver;                             % row address(es) of driver group(s) (e.g., NO3)
                    intraODEinput.repmat_looky_source       = repmat(ODEinput.looky_MigrationFlux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_looky_destiny      = repmat(ODEinput.looky_MigrationFlux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
                    intraODEinput.repmat_grp_row          	= repmat(grp_row', [1, num_fluxes_migration]);                    % addressses; (2D matrix: num_grps X num_fluxes_migration)
                    Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    %	Fluxes_t.flux_import_t	-->flux INTO each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does NOT include domain import (see flux_domain_import_t for domain input rates)
                    %	Fluxes_t.flux_export_t	-->flux OUT OF each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does include domain export
                    %	Fluxes_t.flux_net_t     -->net biomass flux into(+) or out of(-) each box; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: does include domain import & export
                    flux_net_migration_t                    = Fluxes_t.flux_net_t;                      % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
                    flux_domain_ImportDriver_migration_t	= Fluxes_t.flux_domain_import_driver_t;     % net biomass flux OF DRIVER GRP across domain boundary INTO each box (all other grps = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: migration import driver could be used for specialized simulations, but deactived by default
                    % -----------------------------------------------------
                    
                    % build up matrices of subdaily physical fluxes -------
                    %       physical volume fluxes into(+) or out of(-) each box (m3/d)
                    build_VolumeFlux_advection_subdaily(subdaily_loop, :)   = ADVECTION_compact_t;          % (m3/d); (horizontal vector: 1 X num_fluxes_advection)
                    build_VolumeFlux_HorizMixing_subdaily(subdaily_loop, :)	= HORIZONTALMIXING_compact_t;	% (m3/d); (horizontal vector: 1 X num_fluxes_HorizontalMixing)
                    build_VolumeFlux_VertMixing_subdaily(subdaily_loop, :)	= VERTICALMIXING_compact_t;     % (m3/d); (horizontal vector: 1 X num_fluxes_VerticalalMixing)
                    
                    %       physical driver fluxes into(+) or out of(-) each box (mmoles N/m3/d)
                    build_ImportDriver_advection_subdaily(subdaily_loop, :, :)      = reshape(flux_domain_ImportDriver_advection_t, [1, num_grps, num_boxes]);      % (mmoles N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    build_ImportDriver_HorizMixing_subdaily(subdaily_loop, :, :)	= reshape(flux_domain_ImportDriver_HorizMixing_t, [1, num_grps, num_boxes]);	% (mmoles N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    build_ImportDriver_migration_subdaily(subdaily_loop, :, :)      = reshape(flux_domain_ImportDriver_migration_t, [1, num_grps, num_boxes]);      % (mmoles N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    
                    %       net biomass flux into(+) or out of(-) each box (mmole N/m3/d)
                    build_flux_advection_subdaily(subdaily_loop, :, :)      = flux_net_advection_t;   % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    build_flux_HorizMixing_subdaily(subdaily_loop, :, :)	= flux_net_HorizMixing_t; % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    build_flux_VertMixing_subdaily(subdaily_loop, :, :)     = flux_net_VertMixing_t;  % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    build_flux_sinking_subdaily(subdaily_loop, :, :)        = flux_net_sinking_t;     % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    build_flux_migration_subdaily(subdaily_loop, :, :)      = flux_net_migration_t;   % (mmole N/m3/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                    % -----------------------------------------------------
                    
                end % (subdaily_loop)
                
                
                % convert physical flux rates to (t N/d) ------------------
                build_ImportDriver_advection_subdaily       = build_ImportDriver_advection_subdaily .* current_BoxVolume_day;                           % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_ImportDriver_advection_subdaily       = (build_ImportDriver_advection_subdaily * ECOTRANphysics.atomic_mass_N)	/ 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                build_ImportDriver_HorizMixing_subdaily     = build_ImportDriver_HorizMixing_subdaily .* current_BoxVolume_day;                         % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_ImportDriver_HorizMixing_subdaily     = (build_ImportDriver_HorizMixing_subdaily * ECOTRANphysics.atomic_mass_N)	/ 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                build_ImportDriver_migration_subdaily       = build_ImportDriver_migration_subdaily .* current_BoxVolume_day;                           % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_ImportDriver_migration_subdaily       = (build_ImportDriver_migration_subdaily * ECOTRANphysics.atomic_mass_N)	/ 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                build_flux_advection_subdaily               = build_flux_advection_subdaily .* current_BoxVolume_day;                                   % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_advection_subdaily               = (build_flux_advection_subdaily * ECOTRANphysics.atomic_mass_N)            / 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                build_flux_HorizMixing_subdaily           	= build_flux_HorizMixing_subdaily .* current_BoxVolume_day;                                 % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_HorizMixing_subdaily          	= (build_flux_HorizMixing_subdaily * ECOTRANphysics.atomic_mass_N)          / 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                build_flux_VertMixing_subdaily            	= build_flux_VertMixing_subdaily .* current_BoxVolume_day;                                  % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_VertMixing_subdaily            	= (build_flux_VertMixing_subdaily * ECOTRANphysics.atomic_mass_N)           / 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                build_flux_sinking_subdaily                 = build_flux_sinking_subdaily .* current_BoxVolume_day;                                     % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_sinking_subdaily                 = (build_flux_sinking_subdaily * ECOTRANphysics.atomic_mass_N)              / 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                build_flux_migration_subdaily               = build_flux_migration_subdaily .* current_BoxVolume_day;                                   % (mmole N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                build_flux_migration_subdaily               = (build_flux_migration_subdaily * ECOTRANphysics.atomic_mass_N)            / 1000000000;	% (t N/d); (3D matrix: num_subdaily X num_grps X num_boxes DESTINATION)
                
                % convert consumption Q_cp rates to (t N/d) ---------------
                for box_loop = 1:num_boxes
                    currentBox                                  = target_box(box_loop);
                    BoxVolume_current                           = current_BoxVolume_day(:, :, currentBox)';                 % (m3); (2D matrix: num_grps X num_subdaily); NOTE transpose
                    BoxVolume_current                           = reshape(BoxVolume_current, [num_grps, 1, num_subdaily]);	% (m3); (3D matrix: num_grps X 1 X num_subdaily)
                    BoxVolume_current                           = repmat(BoxVolume_current, [1, num_grps, 1]);              % (m3); (3D matrix: num_grps X num_grps X num_subdaily)

                    Q_cp_subdaily(:, :, :, currentBox)          = current_Q_cp_day(:, :, :, currentBox);                    % (mmole N/m3/d); (4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
                    
% %                     % convert units (mmoles N/m3/d) to (t N/d)
% %                     Q_cp_subdaily(:, :, :, currentBox)	= Q_cp_subdaily(:, :, :, currentBox) .* BoxVolume_current;                          % (mmole N/d); (4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
% %                     Q_cp_subdaily(:, :, :, currentBox)	= (Q_cp_subdaily(:, :, :, currentBox) * ECOTRANphysics.atomic_mass_N) / 1000000000;	% (t N/d); (4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
                
                    build_Q_cp_subdaily(:, :, :, currentBox)	= build_Q_cp_subdaily(:, :, :, currentBox) + Q_cp_subdaily(:, :, :, currentBox);     % (mmole N/m3/d) or (t N/d); (4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
                
                end % (box_loop)

% %                 % convert units (mmoles N/m3) to (t N) ---------------
% %                 current_biomass_day                         = current_biomass_day           .* current_BoxVolume_day;                          	% (mmole N); (3D matrix: num_subdaily X num_grps X num_boxes)
% %                 current_biomass_day                         = (current_biomass_day          .* ECOTRANphysics.atomic_mass_N)	/ 1000000000;   % (t N); (3D matrix: num_subdaily X num_grps X num_boxes)
% %                 
% %                 current_production_day                      = current_production_day        .* current_BoxVolume_day;                           % (mmole N/10km2/d); (3D matrix: num_subdaily X num_grps X num_boxes)
% %                 current_production_day                      = (current_production_day       .* ECOTRANphysics.atomic_mass_N)    / 1000000000;   % (t N/10km2/d); (3D matrix: num_subdaily X num_grps X num_boxes)
% %                 
% %                 current_ConsumptionRates_day                = current_ConsumptionRates_day  .* current_BoxVolume_day;                           % (mmole N); (3D matrix: num_subdaily X num_grps X num_boxes)
% %                 current_ConsumptionRates_day              	= (current_ConsumptionRates_day  * ECOTRANphysics.atomic_mass_N)	/ 1000000000;	% (t N); (3D matrix: num_subdaily X num_grps X num_boxes)
% %                 % --------------------------------------------------
                

                % build up daily mean biomasses, consumption fluxes, and physical fluxes
                %       NOTE: average over all hours in the current day
                %       NOTE: unts are either (mmoles N/m3) or (t N)
                build_biomass_dailymean(counter, :, :)              = mean(current_biomass_day, 1);                 % (mmoles N/m3);	(3D matrix: num_dates X num_grps X num_boxes)
                build_production_dailymean(counter, :, :)       	= mean(current_production_day, 1);              % (mmoles N/m3/d);	(3D matrix: num_dates X num_grps X num_boxes)
                build_ConsumptionRates_dailymean(counter, :, :)     = mean(current_ConsumptionRates_day, 1);        % (mmoles N/m3/d);  (3D matrix: num_dates X num_grps X num_boxes)
                
                % daily mean consumption matrices (Q_cp), (mmoles N/m3/d) or (t N/d)
                for box_loop = 1:num_boxes
                    currentBox                  = target_box(box_loop);
                    build_Q_cp_dailymean(:, :, counter, currentBox)      = mean(Q_cp_subdaily(:, :, :, currentBox), 3);     % (mmoles N/m3/d) or (t N/d);	(4D matrix: num_grps X num_grps X num_dates X num_boxes)
                end % (box_loop)
                
                % daily mean physical fluxes -----
                build_VolumeFlux_advection_dailymean(counter, :)        = mean(build_VolumeFlux_advection_subdaily, 1);     % (m3/d);	(3D matrix: num_dates X num_fluxes_advection)
             	build_VolumeFlux_HorizMixing_dailymean(counter, :)      = mean(build_VolumeFlux_HorizMixing_subdaily, 1);	% (m3/d);	(3D matrix: num_dates X num_fluxes_HorizontalMixing)
              	build_VolumeFlux_VertMixing_dailymean(counter, :)       = mean(build_VolumeFlux_VertMixing_subdaily, 1);	% (m3/d);	(3D matrix: num_dates X num_fluxes_VerticalMixing)
                
                %   NOTE: net biomass flux into(+) or out of(-) each box (t N/d)
                build_ImportDriver_advection_dailymean(counter, :, :)	= mean(build_ImportDriver_advection_subdaily, 1);	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
             	build_ImportDriver_HorizMixing_dailymean(counter, :, :)	= mean(build_ImportDriver_HorizMixing_subdaily, 1);	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
              	build_ImportDriver_migration_dailymean(counter, :, :)	= mean(build_ImportDriver_migration_subdaily, 1);	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                
                build_advection_dailymean(counter, :, :)                = mean(build_flux_advection_subdaily, 1);           % (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_HorizMixing_dailymean(counter, :, :)              = mean(build_flux_HorizMixing_subdaily, 1);      	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_VertMixing_dailymean(counter, :, :)               = mean(build_flux_VertMixing_subdaily, 1);       	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_sinking_dailymean(counter, :, :)                  = mean(build_flux_sinking_subdaily, 1);          	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_migration_dailymean(counter, :, :)                = mean(build_flux_migration_subdaily, 1);          	% (t N/d);	(3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                % ---------------------------------------------------------
                
                
                % build up midnight & noon biomasses, consumption fluxes, and physical fluxes 
                current_hour                                    = current_matlabdate_day(:, 4);
                temp_hour                                       = abs(12 - current_hour);
                looky_midnight                                  = find(current_hour == min(current_hour));
                looky_noon                                      = find(temp_hour == min(temp_hour));
                
                % midnight & noon biomasses (t N)
                build_biomass_midnight(counter, :, :)           = current_biomass_day(looky_midnight(1), :, :);             % (mmoles N/m3) or (t N);            (3D matrix: num_dates X num_grps X num_boxes); NOTE: select the first looky_midnight in case there are hour ties (there should not be ties)
                build_biomass_noon(counter, :, :)               = current_biomass_day(looky_noon(1), :, :);                 % (mmoles N/m3) or (t N);            (3D matrix: num_dates X num_grps X num_boxes); NOTE: select the first looky_noon in case there are hour ties (there should not be ties)
                
                % midnight & noon production rates (t N/10km2/d)
                build_production_midnight(counter, :, :)      	= current_production_day(looky_midnight(1), :, :);        	% (mmoles N/m3/d) or (t N/10km2/d);	(3D matrix: num_dates X num_grps X num_boxes); NOTE: select the first looky_midnight in case there are hour ties (there should not be ties)
                build_production_noon(counter, :, :)           	= current_production_day(looky_noon(1), :, :);            	% (mmoles N/m3/d) or (t N/10km2/d);	(3D matrix: num_dates X num_grps X num_boxes); NOTE: select the first looky_noon in case there are hour ties (there should not be ties)
                
                % midnight & noon ConsumptionRates (mmoles N/m3/d)
                build_ConsumptionRates_midnight                 = current_ConsumptionRates_day(looky_midnight(1), :, :);	% (mmoles N/m3/d);  (3D matrix: num_dates X num_grps X num_boxes)
                build_ConsumptionRates_noon                     = current_ConsumptionRates_day(looky_noon(1), :, :);        % (mmoles N/m3/d);  (3D matrix: num_dates X num_grps X num_boxes)

                % midnight & noon consumption matrices (Q_cp); (mmoles N/m3/d) or (t N/d)
                for box_loop = 1:num_boxes
                    currentBox                                      = target_box(box_loop);
                    build_Q_cp_midnight(:, :, counter, currentBox)	= Q_cp_subdaily(:, :, looky_midnight(1), currentBox);        % (mmoles N/m3/d) or (t N/d);          (4D matrix: num_grps X num_grps X num_dates X num_boxes)
                    build_Q_cp_noon(:, :, counter, currentBox)    	= Q_cp_subdaily(:, :, looky_noon(1), currentBox);            % (mmoles N/m3/d) or (t N/d);          (4D matrix: num_grps X num_grps X num_dates X num_boxes)
                end % (box_loop)
                
                % midnight & noon physical fluxes; (t N/d)
                %	NOTE: net biomass flux into(+) or out of(-) each box
                %	NOTE: physical fluxes should be constant for all timepoints in any given day; these should equal _dailymean values
                build_advection_midnight(counter, :, :)         = build_flux_advection_subdaily(looky_midnight(1), :, :);   % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_HorizMixing_midnight(counter, :, :)       = build_flux_HorizMixing_subdaily(looky_midnight(1), :, :); % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_VertMixing_midnight(counter, :, :)        = build_flux_VertMixing_subdaily(looky_midnight(1), :, :);  % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_sinking_midnight(counter, :, :)           = build_flux_sinking_subdaily(looky_midnight(1), :, :);     % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_migration_midnight(counter, :, :)         = build_flux_migration_subdaily(looky_midnight(1), :, :);   % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                
                build_advection_noon(counter, :, :)             = build_flux_advection_subdaily(looky_noon(1), :, :);       % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_HorizMixing_noon(counter, :, :)           = build_flux_HorizMixing_subdaily(looky_noon(1), :, :);     % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_VertMixing_noon(counter, :, :)            = build_flux_VertMixing_subdaily(looky_noon(1), :, :);      % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_sinking_noon(counter, :, :)               = build_flux_sinking_subdaily(looky_noon(1), :, :);         % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                build_migration_noon(counter, :, :)             = build_flux_migration_subdaily(looky_noon(1), :, :);       % (t N/d);          (3D matrix: num_dates X num_grps X num_boxes DESTINATION)
                % ---------------------------------------------------------
                
            end % (if num_subdaily > 0)
            
        end % (day_loop)
    end % (month_loop)
end % (year_loop)
% *************************************************************************





% *************************************************************************
% STEP 4: calculate overall means of the target year(s) & months(s)--------
% step 4a: biomass (mmoles N/m3) or (t N) ---------------------------------
mean_biomass_daily              = squeeze(mean(build_biomass_dailymean, 1));                % (mmoles N/m3) or (t N);    (2D matrix: num_grps X num_boxes); NOTE squeeze
mean_biomass_midnight           = squeeze(mean(build_biomass_midnight, 1));                 % (mmoles N/m3) or (t N);    (2D matrix: num_grps X num_boxes); NOTE squeeze
mean_biomass_noon               = squeeze(mean(build_biomass_noon, 1));                     % (mmoles N/m3) or (t N);    (2D matrix: num_grps X num_boxes); NOTE squeeze
% -------------------------------------------------------------------------


% step 4a: production rates (t N/10km2/d) ---------------------------------
mean_production_daily              = squeeze(mean(build_production_dailymean, 1));      	% (mmoles N/m3/d) or (t N/10km2/d);	(2D matrix: num_grps X num_boxes); NOTE squeeze
mean_production_midnight           = squeeze(mean(build_production_midnight, 1));          	% (mmoles N/m3/d) or (t N/10km2/d);    (2D matrix: num_grps X num_boxes); NOTE squeeze
mean_production_noon               = squeeze(mean(build_production_noon, 1));             	% (mmoles N/m3/d) or (t N/10km2/d);    (2D matrix: num_grps X num_boxes); NOTE squeeze
% -------------------------------------------------------------------------


% step 4b: consumption rates (mmoles N/m3/d) or (t N/d) -------------------
mean_ConsumptionRates_daily     = squeeze(mean(build_ConsumptionRates_dailymean, 1));       % (mmoles N/m3/d) or (t N/d);  (2D matrix: num_grps X num_boxes); NOTE squeeze
mean_ConsumptionRates_midnight	= squeeze(mean(build_ConsumptionRates_midnight, 1));        % (mmoles N/m3/d) or (t N/d);  (2D matrix: num_grps X num_boxes); NOTE squeeze
mean_ConsumptionRates_noon      = squeeze(mean(build_ConsumptionRates_noon, 1));            % (mmoles N/m3/d) or (t N/d);  (2D matrix: num_grps X num_boxes); NOTE squeeze
% -------------------------------------------------------------------------


% step 4c: consumption matrices Q_cp (mmoles N/m3/d) or (t N/d) ------------------------------
mean_Q_cp_daily         = squeeze(mean(build_Q_cp_dailymean, 3));   	% (mmoles N/m3/d) or (t N/d);	(3D matrix: num_grps X num_grps X num_boxes)
mean_Q_cp_midnight      = squeeze(mean(build_Q_cp_midnight, 3));        % (mmoles N/m3/d) or (t N/d);	(3D matrix: num_grps X num_grps X num_boxes)
mean_Q_cp_noon          = squeeze(mean(build_Q_cp_noon, 3));            % (mmoles N/m3/d) or (t N/d);	(3D matrix: num_grps X num_grps X num_boxes)

counter_repmat          = repmat(counter, [num_grps, num_grps, num_subdaily, num_boxes]); % (4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
mean_Q_cp_subdaily      = build_Q_cp_subdaily    ./ counter_repmat;     % (mmoles N/m3/d) or (t N/d);	(4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
% -------------------------------------------------------------------------


% step 4d: physical fluxes (t N/d) ----------------------------------------
mean_VolumeFlux_advection_daily     = mean(build_VolumeFlux_advection_dailymean, 1); % (m3/d); (horizontal vector: 1 X num_fluxes_advection)
mean_VolumeFlux_HorizMixing_daily	= mean(build_VolumeFlux_HorizMixing_dailymean, 1); % (m3/d); (horizontal vector: 1 X num_fluxes_HorizontalMixing)
mean_VolumeFlux_VertMixing_daily	= mean(build_VolumeFlux_VertMixing_dailymean, 1); % (m3/d); (horizontal vector: 1 X num_fluxes_VerticalMixing)

mean_ImportDriver_advection_daily	= squeeze(mean(build_ImportDriver_advection_dailymean, 1)); % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_ImportDriver_HorizMixing_daily	= squeeze(mean(build_ImportDriver_HorizMixing_dailymean, 1)); % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_ImportDriver_migration_daily	= squeeze(mean(build_ImportDriver_migration_dailymean, 1)); % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze

mean_advection_daily            = squeeze(mean(build_advection_dailymean, 1));          % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_HorizMixing_daily          = squeeze(mean(build_HorizMixing_dailymean, 1));        % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_VertMixing_daily           = squeeze(mean(build_VertMixing_dailymean, 1));         % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_sinking_daily              = squeeze(mean(build_sinking_dailymean, 1));            % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_migration_daily            = squeeze(mean(build_migration_dailymean, 1));          % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze

mean_advection_midnight      	= squeeze(mean(build_advection_midnight, 1));           % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_HorizMixing_midnight    	= squeeze(mean(build_HorizMixing_midnight, 1));         % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_VertMixing_midnight     	= squeeze(mean(build_VertMixing_midnight, 1));          % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_sinking_midnight          	= squeeze(mean(build_sinking_midnight, 1));             % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_migration_midnight      	= squeeze(mean(build_migration_midnight, 1));           % (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze

mean_advection_noon             = squeeze(mean(build_advection_noon, 1));            	% (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_HorizMixing_noon           = squeeze(mean(build_HorizMixing_noon, 1));         	% (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_VertMixing_noon            = squeeze(mean(build_VertMixing_noon, 1));            	% (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_sinking_noon               = squeeze(mean(build_sinking_noon, 1));              	% (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
mean_migration_noon             = squeeze(mean(build_migration_noon, 1));             	% (t N/d);          (2D matrix: num_grps X num_boxes DESTINATION); NOTE squeeze
% *************************************************************************





% *************************************************************************
% STEP 5: pack results for export------------------------------------------
% step 5a: biomass (mmoles N/m3) or (t N) ---------------------------------
FluxSummary.mean_biomass_daily              = mean_biomass_daily;                   % (mmoles N/m3) or (t N);	(2D matrix: num_grps X num_boxes)
FluxSummary.mean_biomass_midnight           = mean_biomass_midnight;                % (mmoles N/m3) or (t N);	(2D matrix: num_grps X num_boxes)
FluxSummary.mean_biomass_noon               = mean_biomass_noon;                    % (mmoles N/m3) or (t N);	(2D matrix: num_grps X num_boxes)

FluxTimeseries.biomass_dailymean            = build_biomass_dailymean;              % (mmoles N/m3) or (t N);	(3D matrix: num_dates X num_grps X num_boxes)
FluxTimeseries.biomass_midnight             = build_biomass_midnight;               % (mmoles N/m3) or (t N);	(3D matrix: num_dates X num_grps X num_boxes)
FluxTimeseries.biomass_noon                 = build_biomass_noon;                   % (mmoles N/m3) or (t N);	(3D matrix: num_dates X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 5a: production rates (mmoles N/m3/d) or (t N/10km2/d) --------------
FluxSummary.mean_production_daily              = mean_production_daily;            	% (mmoles N/m3/d) or (t N/10km2/d);	(2D matrix: num_grps X num_boxes)
FluxSummary.mean_production_midnight           = mean_production_midnight;      	% (mmoles N/m3/d) or (t N/10km2/d);	(2D matrix: num_grps X num_boxes)
FluxSummary.mean_production_noon               = mean_production_noon;           	% (mmoles N/m3/d) or (t N/10km2/d);	(2D matrix: num_grps X num_boxes)

FluxTimeseries.production_dailymean            = build_production_dailymean;       	% (mmoles N/m3/d) or (t N/10km2/d);	(3D matrix: num_dates X num_grps X num_boxes)
FluxTimeseries.production_midnight             = build_production_midnight;       	% (mmoles N/m3/d) or (t N/10km2/d);	(3D matrix: num_dates X num_grps X num_boxes)
FluxTimeseries.production_noon                 = build_production_noon;           	% (mmoles N/m3/d) or (t N/10km2/d);	(3D matrix: num_dates X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 5b: consumption rates (mmoles N/m3/d) or (t N/d) -------------------
FluxSummary.mean_ConsumptionRates_daily     = mean_ConsumptionRates_daily;          % (mmoles N/m3/d) or (t N/d);  (2D matrix: num_grps X num_boxes)
FluxSummary.mean_ConsumptionRates_midnight	= mean_ConsumptionRates_midnight;       % (mmoles N/m3/d) or (t N/d);  (2D matrix: num_grps X num_boxes)
FluxSummary.mean_ConsumptionRates_noon      = mean_ConsumptionRates_noon;           % (mmoles N/m3/d) or (t N/d);  (2D matrix: num_grps X num_boxes)

FluxTimeseries.ConsumptionRates_dailymean	= build_ConsumptionRates_dailymean;     % (mmoles N/m3/d) or (t N/d);  (3D matrix: num_dates X num_grps X num_boxes)
FluxTimeseries.ConsumptionRates_midnight	= build_ConsumptionRates_midnight;      % (mmoles N/m3/d) or (t N/d);  (3D matrix: num_dates X num_grps X num_boxes)
FluxTimeseries.ConsumptionRates_noon        = build_ConsumptionRates_noon;          % (mmoles N/m3/d) or (t N/d);  (3D matrix: num_dates X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 5c: consumption matrices Q_cp (mmoles N/m3/d) or (t N/d) -----------
FluxSummary.mean_Q_cp_daily                 = mean_Q_cp_daily;                      % (mmoles N/m3/d) or (t N/d);	(2D matrix: num_grps X num_grps X num_boxes)
FluxSummary.mean_Q_cp_midnight              = mean_Q_cp_midnight;                   % (mmoles N/m3/d) or (t N/d);	(2D matrix: num_grps X num_grps X num_boxes)
FluxSummary.mean_Q_cp_noon                  = mean_Q_cp_noon;                       % (mmoles N/m3/d) or (t N/d);	(2D matrix: num_grps X num_grps X num_boxes)
FluxSummary.mean_Q_cp_subdaily              = mean_Q_cp_subdaily;                   % (mmoles N/m3/d) or (t N/d);	(4D matrix: num_grps X num_grps X num_subdaily X num_boxes)
% -------------------------------------------------------------------------


% step 5d: physical fluxes (t N/d) ----------------------------------------
FluxSummary.mean_VolumeFlux_advection_daily     = mean_VolumeFlux_advection_daily;      % (m3/d); (horizontal vector: 1 X num_fluxes_advection)
FluxSummary.mean_VolumeFlux_HorizMixing_daily	= mean_VolumeFlux_HorizMixing_daily;	% (m3/d); (horizontal vector: 1 X num_fluxes_HorizontalMixing)
FluxSummary.mean_VolumeFlux_VertMixing_daily	= mean_VolumeFlux_VertMixing_daily;     % (m3/d); (horizontal vector: 1 X num_fluxes_VerticalMixing)

FluxSummary.mean_ImportDriver_advection_daily	= mean_ImportDriver_advection_daily;	% (t N/d); (2D matrix: num_grps X num_boxes DESTINATION)
FluxSummary.mean_ImportDriver_HorizMixing_daily	= mean_ImportDriver_HorizMixing_daily;	% (t N/d); (2D matrix: num_grps X num_boxes DESTINATION)
FluxSummary.mean_ImportDriver_migration_daily	= mean_ImportDriver_migration_daily;	% (t N/d); (2D matrix: num_grps X num_boxes DESTINATION)

FluxSummary.mean_advection_daily                = mean_advection_daily;                 % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_HorizMixing_daily              = mean_HorizMixing_daily;               % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_VertMixing_daily               = mean_VertMixing_daily;                % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_sinking_daily                  = mean_sinking_daily;                   % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_migration_daily                = mean_migration_daily;                 % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);

FluxSummary.mean_advection_midnight             = mean_advection_midnight;              % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_HorizMixing_midnight           = mean_HorizMixing_midnight;            % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_VertMixing_midnight            = mean_VertMixing_midnight;             % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_sinking_midnight               = mean_sinking_midnight;                % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_migration_midnight             = mean_migration_midnight;              % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);

FluxSummary.mean_advection_noon                 = mean_advection_noon;                  % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_HorizMixing_noon               = mean_HorizMixing_noon;                % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_VertMixing_noon                = mean_VertMixing_noon;                 % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_sinking_noon                   = mean_sinking_noon;                    % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
FluxSummary.mean_migration_noon                 = mean_migration_noon;                  % (t N/d); (2D matrix: num_grps X num_boxes DESTINATION);
% *************************************************************************


% end m-file***************************************************************