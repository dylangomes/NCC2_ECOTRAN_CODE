function [ECOTRANphysics] = f_ECOTRANphysics_NCC2_upwelling_09042022(PHYSICSinput, upwelling_driver,ShowOutput,CUTI_YEARS)
% returns Advection and Mixing exchanges between boxes
% use for NCC upwelling setting
%
% calls:
  %       f_read_NWPO3data_02012022           read wind data from station NWPO3 on Newport South Jetty (no smoothing)	(m/s)                           (vertical vector)
%       f_UpwellingIndex                    calculat Bakun Upwelling Index from wind data                           (m3/s per 100m of coastline)	(vertical vector)
%       f_smooth                            mean smoothing over provided smoothing_window (in days before and after time-point)
%       f_read_ERDdata_05122021             calculate daily median upwelling intensities from ERD products (http://http://www.pfeg.noaa.gov/products/PFELData/upwell/6_hourly/upwell45N125W)
%       [deactivated] f_ekmandepth          calculate ekman depth
%       f_prep_ERD_CUTI_08052022           process ERD CUTI timeseries for upwelling flux cropped to datestart:dateend time period	(m2/s)	(vertical vector: num_t X 1)
%           f_smooth                            mean smoothing over provided smoothing_window (in days before and after time-point)
%       round2                              round to a specified number of decimal places
%       f_CompactFluxTimeSeries_11182019	compact flux time series when arranged as (3D matrix: time X source box X destiny box)
%       f_UnCompactFluxTimeSeries_12112019  UnCompact a flux time series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole
%           f_calcNetFlux_12112019              calculate net flux into and net flux out of each model box and across outer domain boundaries
%       f_EvaluateFluxBalance_11262021      examine for flux time-series for imbalances IN & OUT of invidual boxes and IN & OUT of the overall domain
%       calcur_res.mat                      (dataset with monthly mean salinity, temperature, and NO3+NO2 data)
%       f_LightIntensity_12112020           instantaneous (W/m2), daily mean averaged across 24 h (W m^-2 h^-1), & daily integrated (W m^-2 d^-1) solar raditation at ocean surface; (vertical vector: num_t X 1)
%
% takes:
  %       PHYSICSinput:
  %           datestart                       time-series start date in matlab double format (scalar)
%           dateend                         time-series end date in matlab double format (scalar)
%           dt                              time-step (days) (scalar)
%           t_grid                          time-series dates in matlab double format (vertical vector: num_t X 1)
%           smoothing_window                window for smoothing before and after time-point (e.g., 2 = a window of 5 days centered on time-point t)
%       upwelling_driver                    specify flux time-series to use ('Brink_BUI', 'NWPO3_BUI', 'ERD_BUI', 'ERD_CUTI', 'Fake_Upwelling')
%
% returns 
%       ECOTRANphysics:
  %           AdvectionMatrix         advection rate = volume transported per day; (m3/d)
  %                                       (3D matrix: time X source box X destination box)
  %                                       (values can be positive or negative)
  %           VerticalMixMatrix       vertical mixing rate = volume mixed per day; (m3/d)
  %                                       (3D matrix: time X source box X destination box)
  %                                       (values can be positive or negative; sign defines direction of net mixing)
  %
  % NOTE: What happens in cases when there is a complete advective washout of a spatial box?)
%       I added cap to advection rate to prevent >100% box washout (but do not account for river flux)
%       Added warning for days when AdvectionMatrix > BoxVolume; see ECOTRANphysics.BoxWashoutTime
%
% revision date: 9-4-2022
%           8/13/2022 carrying through SinkLink terms for benthic NH4 needs in 2D cross-shelf model
%           8/22/2022 adding null temperature time-series
%           9/4/2022 changed function name, changed directories


% *************************************************************************
  % STEP 1: basic domain parameters (NCC upwelling)--------------------------
  fname_PhysicalModel      = mfilename; % save name of this m-file to keep in saved model results
  if ShowOutput
  display(['Running: ' fname_PhysicalModel])
  end
  
  % directories (all gathered in one place)
  % NutrientFile_directory      = '/Users/jamesruzicka/Documents/10_ECOTRAN_code/7_NCC_code/NCC_physics/'; % directory: NH-Line nutrient climatology
  % NWP03winddata_directory     = '/Users/jamesruzicka/Documents/10_ECOTRAN_code/7_NCC_code/NCC_physics/nwpo3h_data_raw/'; % rdirectory: aw south jetty wind data
  % ERD_BUI_directory           = '/Users/jamesruzicka/Documents/10_ECOTRAN_code/7_NCC_code/NCC_physics/ERD_data_raw/ERD_45N125W_1967-2015.csv'; % file: upwelling from ERD BUI product
  % ERD_CUTI_directory          = '/Users/jamesruzicka/Documents/10_ECOTRAN_code/7_NCC_code/NCC_physics/ERD_CUTI_data/CUTI_daily.nc'; % file: upwelling from ERD CUTI product
  
  % Ebi laptop
  NutrientFile_directory      = '/Users/dgome/Documents/NCC2_ECOTRAN_CODE/EcoTran_Code/'; % directory: NH-Line nutrient climatology
  NWP03winddata_directory     = '/Users/jimsebi/Documents/10_ECOTRAN_code/7_NCC_code/NCC_physics/nwpo3h_data_raw/'; % directory: raw south jetty wind data
  ERD_BUI_directory           = '/Users/jimsebi/Documents/10_ECOTRAN_code/7_NCC_code/NCC_physics/ERD_data_raw/ERD_45N125W_1967-2015.csv'; % file: upwelling from ERD BUI product
  ERD_CUTI_directory          = '/Users/dgome/Documents/NCC2_ECOTRAN_CODE/EcoTran_Code/CUTI/CUTI_daily.nc'; % file: upwelling from ERD CUTI product
  % ERD_CUTI_AVG_directory      = '/Users/dgome/Documents/NCC2_ECOTRAN_CODE/EcoTran_Code/CUTI/CUTI_AVERAGE_daily.csv'; % file: upwelling from ERD CUTI product
  ERD_CONST_directory          = '/Users/dgome/Documents/NCC2_ECOTRAN_CODE/EcoTran_Code/CUTI/CUTI_daily0.01.nc'; % file: upwelling from ERD CUTI product
  
  % -------------------------------------------------------------------------
    
    
    % step 1a: physical & other operating parameters --------------------------
    latitude                    = 44.6516667;              % latitude of NH line (decimal degress north) based on LTOP transects
  rho_water                   = 1000;                    % seawater density; (kg/m3)
  rho_air                     = 1.220;                   % air density (kg/m3) (Schwing et al. 1996)
  coriollis                   = 1.025 * 10^(-4);         % coriollis factor (1/s); (NH-Line NNPPZD model: Ruzicka et al. 2011)
  midmonth_doy                = [15 45 74 105 135 166 196 227 258 288 319 349]; % mid-month day-of-year
  reference_temperture_season	= [5 9]; % start and end months to define reference temperature
  smoothing_window            = PHYSICSinput.smoothing_window; % window for smoothing before and after time-point (e.g., 2 = a window of 5 days centered on time-point t)
  target_latitudes            = PHYSICSinput.target_latitudes; % QQQ chosse only 1 latitude for now; FFF in future, choose 1 or more target latitude(s); [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
  % -------------------------------------------------------------------------
    
    
    % step 1b: time domain and vectors ----------------------------------------
  datestart    = PHYSICSinput.datestart;  % enter starting date
  dateend      = PHYSICSinput.dateend;    % enter ending date
  dt           = PHYSICSinput.dt;         % t-step (d)
  t_grid       = PHYSICSinput.t_grid;     % t_grid; (vertical vector)
  num_t        = length(t_grid);          % length of t_grid; (scaler)
  min_t        = min(t_grid);
  max_t        = max(t_grid);
  num_years    = ceil(max_t/365);
  
  % build up t_grid_midmonth to cover the full span of simulation years
  build_midmonth_doy = [];
  for year_loop = 1:num_years
  build_midmonth_doy = [build_midmonth_doy (midmonth_doy + (365 * (year_loop - 1)))];  % (horizontal vector: 1 X (12*years)
                                                                                          end % (year_loop)
                                                                                          t_grid_midmonth     = [1 build_midmonth_doy (365*num_years)]'; % day; (vertical vector: (12*num_years + 2) X 1)
% -------------------------------------------------------------------------


% step 1c: Box ID info ----------------------------------------------------
num_boxes                   = 5;        % I, II, III, IV, V
num_AdvecFluxes             = 6;        % advective fluxes: A_12, A_13, A_24, A_35, A_4offshore, A_5offshore
num_MixFluxes               = 2;        % vertical mix fluxes: w_23, w_45
num_HorizMixFluxes          = 6;        % horizontal mix fluxes: u_12 u_13 u_24 u_35 u_4ocean u_5ocean

looky_SubSurfaceBoxes       = [3, 5];
looky_SurfaceBoxes          = [2, 4];
looky_OffshoreSurfaceBox    = 4;
looky_OffshoreSubSurfaceBox = 5;
looky_RiverFluxBox          = 2;        % River input from alongshore; (Box 2 in CGoA downwelling model)
% -------------------------------------------------------------------------


% step 1d: spatial dimensions ---------------------------------------------
%          NOTE: Mixed Layer Depth (only used for light intensity calculations, not used for advection rate calculations)
%          NOTE: box geometry with slanty bottom
maxMLD                      = 15;                                                       % max MLD; (m) % alt 40
minMLD                      = 15;                                                       % min MLD; (m) % alt 10
rangeMLD                    = maxMLD - minMLD;                                          % set minimum mixed layer depth and annual range of the mixed layer depth; (m)
MLD                         = (cos(((t_grid/365)*pi)*2) * (rangeMLD * 0.5) + (minMLD + rangeMLD * 0.5)); % MLD; (m); (vertical vector: time X 1)
BoxLength                   = [10000, 20000, 20000, 20000, 20000];                      % (m); I, II, III, IV, V
BoxLength                   = repmat(BoxLength, [num_t, 1]);                            % (m); (2D matrix: time X num_boxes)
BoxWidth                    = [1, 1, 1, 1, 1];                                          % (m); I, II, III, IV, V
BoxWidth                    = repmat(BoxWidth, [num_t, 1]);                             % (m); (2D matrix: time X num_boxes)
BoxHeight_inshore           = [repmat(0,  [num_t, 1]), MLD, (30-MLD),  MLD, (140-MLD)]; % (m); (2D matrix: time X num_boxes); I, II, III, IV, V
BoxHeight_offshore          = [repmat(30, [num_t, 1]), MLD, (140-MLD), MLD, (250-MLD)]; % (m); (2D matrix: time X num_boxes); I, II, III, IV, V
BoxHeight                   = 0.5 * (BoxHeight_inshore + BoxHeight_offshore);           % mean box height; (m); (2D matrix: time X num_boxes); I, II, III, IV, V
BoxArea                     = BoxHeight .* BoxLength;                                   % (m2); (2D matrix: time X num_boxes); I, II, III, IV, V

BoxFloorArea              	= BoxLength .* BoxWidth;                % (m2); (2D matrix: time X num_boxes SOURCE)

BoxVolume                   = BoxArea .* BoxWidth;                                      % I, II, III, IV, V; (m3); (2D matrix: time X num_boxes);
% *************************************************************************





% *************************************************************************
% STEP 2: conversion factors-----------------------------------------------
atomic_mass_N             = 14.00674;   % (g N / mole N)
atomic_mass_C             = 12.0107;    % (g C / mole C)

Chla_to_N_diatoms         = 2.19;       % (mg Chla mmoles N^-1) (Dickson & Wheeler 1995 L&O 40(3):533-543)
% Chla_to_N_diatoms         = 2.5;        % 2.49 to 2.67 (mg Chla mmoles N^-1) (Chan 1980 J Phycol 16:428-432)
Chla_to_N_dinoflagellates = 0.84;       % 0.73 to 0.95 (mg Chla mmoles N^-1) (Chan 1980 J Phycol 16:428-432)

% C_to_N_phytoplankton      = 7.3;        % general phytoplankton C:N (moles C / moles N) (Geider & LaRoche 2002 Eur. J. Phycol 37:1-17)
C_to_N_phytoplankton      = (106/16);   % Redfield ratio = C:N:P = 106:16:1 ---> C:N = 106:16 = 6.625 mole C/mole N
C_to_N_zooplankton        = 5;          % general zooplankton C:N (moles C moles N^-1) for mesozooplankton (W. Peterson pers com); C:N=5.16-5.26 [Schneider 1990 (Mar Biol 106:219-225)]; C:N=3.9 [Uye 1982]

WWT_to_C                  = 8.77;       % (g WWT / g C); (Steele et al. 2007), NOTE: I could not find this value in 2007 citation, must have gotten value directly from John, value is close to other lit. values; {Strickland, 1966 #1521} says "From the foregoing, it appears that the 10-12 mg algal weight:mg C is the most suitable factor for estimating phytoplankton biomass as ?wet weight?";
mgWWT_per_mmolesN         = (C_to_N_phytoplankton * atomic_mass_C * WWT_to_C); % (mg WWT/mmoles N)

mgram_to_gram             = 1000;       % (mg/g)
gram_to_ton               = 1000000;	% (g/t)
mgram_to_ton              = mgram_to_gram * gram_to_ton; % (mg/t)

W_to_microE               = 4.6;        % convert (W m^-2) to (microE m^-2 s^-1) at 400-700nm (Brock 1981 Ecol Model 14: 1-19)
% *************************************************************************





% *************************************************************************
% STEP 3: load nutrient data (NCC)-----------------------------------------
%         The file calcur_res.mat contains the structure "calcur" with the fields
%         "model_box", "var", and "monthly_mean". To keep usage (i.e., searching)
%         simple, model_box and var use a numeric code:
% 
%           model_box:
%               1: coastal box: everything shoreward of the 60m isobath
%               2: shelf mixed layer (ML) box: the upper 15m between the 60m and 200m isobath
%               3: shelf deep box: 15m to the bottom, between the 60m and 200m isobath
%               4: shelf break ML box: the upper 15m in a 20km swatch offshore of the 200m isobath
%               5: shelf break deep box: 15m to the bottom in a 20km swatch offshore of the 200m isobath
% 
%           var:
%               1: salinity
%               2: temperature
%               3: NO3+NO2
% 
%       The field "monthly_mean" contains 12 values with the monthly means.
%       The structure "calcur" has 15 entries: one for each model box and variable.
% step 3a: load NH-Line nutrient climatology ------------------------------
readNutrientFile           = 'calcur_res.mat'; % monthly mean nutrients
readNutrientFile           = [NutrientFile_directory readNutrientFile];
load(readNutrientFile)
looky_NO3                   = find([calcur.model_box]==1 & [calcur.var] == 3);
box_I_NO3                   = calcur(looky_NO3).monthly_mean;
looky_NO3                   = find([calcur.model_box]==2 & [calcur.var] == 3);
box_II_NO3                  = calcur(looky_NO3).monthly_mean;
looky_NO3                   = find([calcur.model_box]==3 & [calcur.var] == 3);
box_III_NO3                 = calcur(looky_NO3).monthly_mean;
looky_NO3                   = find([calcur.model_box]==4 & [calcur.var] == 3);
box_IV_NO3                  = calcur(looky_NO3).monthly_mean;
looky_NO3                   = find([calcur.model_box]==5 & [calcur.var] == 3);
box_V_NO3                   = calcur(looky_NO3).monthly_mean;

dat_NO3                     = [box_I_NO3' box_II_NO3' box_III_NO3' box_IV_NO3' box_V_NO3']; % monthly mean NO3 + NO2; (micro-moles N/L) = (mmoles N/m3); (matrix: month X box; I, II, III, IV, V)
% -------------------------------------------------------------------------
  
  
  % step 3b: interpolate monthly means to daily values ----------------------
  %          allow for cyclical repeats over several years
%          (t_grid_midmonth was calculated in step 2c)
build_NO3               = [];

for year_loop = 1:num_years
build_NO3 = [build_NO3; dat_NO3];
end

NO3_midmonth            = [mean(dat_NO3([1 12], :)); build_NO3; mean(dat_NO3([1 12], :))]; % nutrient concentration; (mmoles N/m3); (2D matrix: (12*num_years + 2) X num_boxes)
NO3_doy                 = interp1(t_grid_midmonth, NO3_midmonth, t_grid, 'spline'); % (mmoles N/m3); (2D matrix: num_t X num_boxes)
NO3_doy(NO3_doy < 0)	= 0; % filter out negative NO3_doy values created during spline interpolation
% *************************************************************************
  
  
  
  
  
  % *************************************************************************
  % STEP 4: select and prepare advection & mixing parameters-----------------
  
  % step 4a: mixing & river (cross-shelf) fluxes ----------------------------
  
  % river infux (Qr) ................................................
%       (m3/s per 1m BoxWidth)
%       NOTE: all 0 for NCC; FFF may add CR influence in future
Qr      = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% vertical mixing rate ............................................
%       (m/s)
%       NOTE: no direction, so all we are positive
we      = [2.09051350944941e-05, 1.70040101641047e-05, 1.15783391649564e-05, 4.45890782591914e-06, 1.60603701399058e-06, 1.20405100592767e-06, 1.61576927579862e-06, 1.51120326652934e-06, 2.07979462147540e-06, 4.74448045109569e-06, 1.40606441256106e-05, 2.13712842043411e-05];

% horizontal mixing rate ..........................................
%       (m/s)
%       NOTE: not considered for NCC, so all are 0
ue      = [
  0 0 0 0 0 0 0 0 0 0 0 0     % u_12
  0 0 0 0 0 0 0 0 0 0 0 0     % u_13
  0 0 0 0 0 0 0 0 0 0 0 0     % u_24
  0 0 0 0 0 0 0 0 0 0 0 0     % u_35
  0 0 0 0 0 0 0 0 0 0 0 0     % u_4ocean
  0 0 0 0 0 0 0 0 0 0 0 0];   % u_5ocean

% allow for cyclical repeats over several years ...................
build_Qr           = [];
build_we           = [];
build_ue           = [];
for year_loop = 1:num_years
build_Qr           = [build_Qr Qr];  % (horizontal vector: 1 X (12*years))
build_we           = [build_we we];  % (m/s); (horizontal vector: 1 X (12*years))
build_ue           = [build_ue ue];  % (m/s); (2D matrix: num_HorizMixFluxes X (12*years))
end % (year_loop)

% interpolate monthly means to daily values .......................
Qr_midmonth         = [mean(Qr([1 12])) build_Qr mean(Qr([1 12]))]'; % river infux (Qr); (m3/s per 1m BoxWidth); (vertical vector: (12*num_years + 2) X 1)
we_midmonth         = [mean(we(:, [1 12]), 2) build_we mean(we(:, [1 12]), 2)]'; % vertical mixing rate; (m/s); (2D matrix: (12*num_years + 2) X num_VertMixFluxes)
ue_midmonth         = [mean(ue(:, [1 12]), 2) build_ue mean(ue(:, [1 12]), 2)]'; % horizontal mixing rate; (m/s); (2D matrix: (12*num_years + 2) X num_HorizMixFluxes)

Qr_doy              = interp1(t_grid_midmonth, Qr_midmonth, t_grid, 'spline');	% (kg/m/s2); (vertical vector: num_t X 1)
we_doy              = interp1(t_grid_midmonth, we_midmonth, t_grid, 'spline'); % (m/s); (2D matrix: num_t X num_VertMixFluxes)
% we_doy              = we_doy / 2; % NOTE: division by 1/2 already accounted for in original calculation of rates; vertical mixing extends 1/2 into top box & 1/2 into bottom box; (m/s); (2D matrix: num_t X num_VertMixFluxes)
ue_doy              = interp1(t_grid_midmonth, ue_midmonth, t_grid, 'spline'); % (m/s); (2D matrix: num_t X num_HorizMixFluxes)

Qr_doy(Qr_doy < 0)	= 0; % filter out negative ue values created during spline interpolation
ue_doy(ue_doy < 0)	= 0; % filter out negative ue values created during spline interpolation

        
RiverFluxRate       = Qr_doy;           	% (m3/s per 1m BoxWidth); (vertical vector: num_t X 1)
% -------------------------------------------------------------------------


% step 4b: select and prepare upwelling (advection) driver ----------------
switch upwelling_driver

    case 'Brink_BUI' % ----------------------------------------------------
        % Brink Physics raw values (NCC upwelling)
        %	NOTE: monthly means
        %	NOTE: ONshore direction is positive
        
        % wind-driven advection (Qw) ......................................
        %       Qw = tau .* (1/coriollis) .* (1/rho_water) = daily advective flux rate; (m3/s per 1m BoxWidth)
        %       Alongshore wind stress (tau) monthly mean: (N/m^2) = (kg*m/m2/s2) = (kg/m/s2)
        tau     = [0.0680748069736075, 0.0433269605437942, 0.0448631025587659, 0.0176074709095014, -0.0105705908497718, -0.0202353968838113, -0.0298904158691320, -0.0193376707154904, -0.0124022877853344, 0.0219368906632379, 0.0642634881319758, 0.0709770739741683];
        % Qw      = tau_doy .* (1/coriollis) .* (1/1000) * (-1); % daily advective flux rate; (m3/s per 1m BoxWidth); (vertical vector); (NOTE: change of sign to make offshore "POSITIVE")

        % allow for cyclical repeats over several years ...................
        build_tau          = [];
        % build_Qw           = [];
        for year_loop = 1:num_years
            build_tau          = [build_tau tau];  % (horizontal vector: 1 X (12*years))
        %     build_Qw           = [build_Qw Qw];
        end % (year_loop)

        % interpolate monthly means to daily values .......................
        tau_midmonth        = [mean(tau([1 12])) build_tau mean(tau([1 12]))]'; % Alongshore wind stress (tau) monthly mean: (N/m^2) = (kg*m/m2/s2) = (kg/m/s2); (vertical vector: (12*num_years + 2) X 1)
% Qw_midmonth         = [mean(Qw([1 12])) build_Qw mean(Qw([1 12]))]'; % daily advective flux rate; (m3/s per 1m BoxWidth); (NOTE: offshore "POSITIVE"); (vertical vector: (12*num_years + 2) X 1)

tau_doy             = interp1(t_grid_midmonth, tau_midmonth, t_grid, 'spline');	% (kg/m/s2); (vertical vector: num_t X 1)
% Qw_doy              = interp1(t_grid_midmonth, Qw_midmonth, t_grid, 'spline');	% (kg/m/s2); (vertical vector: num_t X 1)

% wind-stress driver -.............................................
Brink_BUI       = tau_doy .* (1/coriollis) .* (1/1000) * (-1); % daily advective flux rate; (m3/s per 1m BoxWidth); (vertical vector: num_t X 1); (NOTE: change of sign to make offshore "POSITIVE")

AdvecFluxRate           = Brink_BUI;            % (m3/s per 1m of BoxWidth); (vertical vector: num_t X 1); positive is upwelling, negative is downwelling
fname_UpwellingDriver	= 'Brink_BUI';          % save name of this driver to keep in saved model results
if ShowOutput
disp('   -->Ken Brink wind stress physics');
end
% end case 'Brink_BUI' ------------------------------------------------
  
  
  case 'NWPO3_BUI' % ---------------------------------------------
  % NWPO3 BUI
%	FFF this needs to be fixed for error when start dates not = Jan 1
WindParams.file_winddata        = NWP03winddata_directory; % raw south jetty wind data
WindParams.dt                   = dt;
WindParams.t_grid               = t_grid;
WindParams.datestart            = datestart;
WindParams.dateend              = dateend;

% calculate median winds at proper time-resolution ........................
%   NOTE: winds are in "ocean convention" ie. direction blowing towards
%   NOTE: median winds over time intervals, wind vectors 1 shorter than time vector
[easterlywinds, northerlywinds]	= f_read_NWPO3data_02012022(WindParams); % (m/s); (vertical vector)

% calculate upwelling index from local NWPO3 wind data ....................
UpwellingParams.northerlywinds  = northerlywinds; % (m/s); (vertical vector)
UpwellingParams.easterlywinds   = easterlywinds; % (m/s); (vertical vector)
UpwellingParams.rho_air         = rho_air; % air density (kg/m3) (Schwing et al. 1996)
UpwellingParams.rho_water       = rho_water; % seawater density; (kg/m3)
UpwellingParams.coriollis       = coriollis; % coriollis factor (1/s); (NH-Line NNPPZD model: Ruzicka et al. 2011)

[UpwellingIndex]                = f_UpwellingIndex(UpwellingParams);
NWPO3_BUI                       = UpwellingIndex.BUI; % (m3/s per 100m of coastline); (vertical vector)
NWPO3_BUI                       = NWPO3_BUI / 100; % (m3/s per 1m of coastline); (vertical vector)
NWPO3_BUI_smooth                = f_smooth(NWPO3_BUI, smoothing_window);
tau_y_NWPO3                     = UpwellingIndex.tau_y; % north-south wind stress; (kg/m/s2); (vertical vector)

AdvecFluxRate           = NWPO3_BUI_smooth;     % (m3/s per 1m of BoxWidth); (vertical vector)
fname_UpwellingDriver	= 'NWPO3_BUI';   % save name of this driver to keep in saved model results
if ShowOutput
disp('   -->NWPO3_BUI upwelling index'); 
end
% end case 'NWPO3_BUI' ------------------------------------------------
  
  
  case 'ERD_BUI' % ------------------------------------------------------
  % ERD BUI
%	get daily medians of ERD data ...........................................
%       from --->>> http://www.pfeg.noaa.gov/products/PFEL/modeled/indices/transports/transports.html
%	NOTE: prepared .csv file covers 01-Jan-1969 through 30-Sept-2015 and has 6-hr resolution
%   NOTE: code will NOT replicate time-series to fill out to dateend
ERDparams.file_ERDdata  = ERD_BUI_directory; % upwelling from ERD product
ERDparams.datestart     = datestart;
ERDparams.dateend       = dateend;
ERDparams.t_grid        = t_grid;
ERDdata                 = f_read_ERDdata_05122021(ERDparams);
% depth_ekman             = f_ekmandepth(ERDdata, 30); % Ekman depth (m)
BUI_ERD                 = ERDdata.ERD_ekman_offshore; % (m3/s per 100m of coastline); (vertical vector: num_t X 1)
BUI_ERD                 = BUI_ERD / 100; % (m3/s per 1m of coastline); (vertical vector: num_t X 1)
BUI_ERD                 = f_smooth(BUI_ERD, smoothing_window); % (m3/s per 1m of coastline); (vertical vector: num_t X 1)

AdvecFluxRate           = BUI_ERD;              % (m3/s per 1m of BoxWidth); (vertical vector)
fname_UpwellingDriver	= 'BUI_ERD';            % save name of this driver to keep in saved model results
if ShowOutput
disp('   -->ERD_BUI upwelling index');      
end
% end case 'ERD_BUI' --------------------------------------------------
  
  
  case 'ERD_CUTI' % -----------------------------------------------------      
  % ERD CUTI
%	get daily CUTI from ERD
%	NOTE: ERDcuti file covers 1-Jan-1988 through 30-June-2022 and has 1-day resolution
%   NOTE: code will replicate time-series to fill out to dateend

if ShowOutput
disp('   -->ERD_CUTI upwelling index');
end

ERD_CUTI_input.readFile_ERD_CUTI    = ERD_CUTI_directory;
ERD_CUTI_input.datestart            = datestart;
ERD_CUTI_input.dateend              = dateend;
ERD_CUTI_input.num_t                = num_t;
ERD_CUTI_input.dt                   = dt;
ERD_CUTI_input.smoothing_window     = smoothing_window; % set moving average smoothing as time points for averaging before and after current time point
ERD_CUTI_input.target_latitudes  	= target_latitudes; % QQQ chosse only 1 latitude for now; FFF in future, choose 1 or more target latitude(s); [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]


ERD_CUTI                            = f_prep_ERD_CUTI_SELECTYEAR_10132022(ERD_CUTI_input,CUTI_YEARS);

AdvecFluxRate                       = ERD_CUTI;        % (m3/s per 1m of BoxWidth); (vertical vector)
fname_UpwellingDriver               = 'ERD_CUTI';       % save name of this driver to keep in saved model results
% end case 'ERD_CUTI' -------------------------------------------------
  
  
  case 'ERD_AVG_CUTI' % -----------------------------------------------------      
  % ERD CUTI
%	get daily CUTI from ERD
%	NOTE: ERDcuti file covers 1-Jan-1988 through 30-June-2022 and has 1-day resolution
%   NOTE: code will replicate time-series to fill out to dateend

if ShowOutput
disp('   -->ERD_CUTI_AVERAGED (repeated) upwelling index');
end

ERD_CUTI_input.readFile_ERD_CUTI    = ERD_CUTI_directory;
ERD_CUTI_input.datestart            = datestart;
ERD_CUTI_input.dateend              = dateend;
ERD_CUTI_input.num_t                = num_t;
ERD_CUTI_input.dt                   = dt;
ERD_CUTI_input.smoothing_window     = smoothing_window; % set moving average smoothing as time points for averaging before and after current time point
ERD_CUTI_input.target_latitudes  	= target_latitudes; % QQQ chosse only 1 latitude for now; FFF in future, choose 1 or more target latitude(s); [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]


ERD_CUTI                            = f_prep_ERD_CUTI_AVGYEAR_10132022(ERD_CUTI_input);


AdvecFluxRate                       = ERD_CUTI;        % (m3/s per 1m of BoxWidth); (vertical vector)
fname_UpwellingDriver               = 'ERD_AVG_CUTI';       % save name of this driver to keep in saved model results
% end case 'ERD_CUTI' -------------------------------------------------
  
  case 'ERD_CONST' % -----------------------------------------------------      
  % ERD CONSTANT
%	get daily ZERO from ERD
%	NOTE: ERDZERO file covers 1-Jan-1988 through 30-June-2022 and has 1-day resolution
%   NOTE: code will replicate time-series to fill out to dateend

if ShowOutput
disp('   --> CONSTANT upwelling index');
end

ERD_CUTI_input.readFile_ERD_CUTI    = ERD_CONST_directory;
ERD_CUTI_input.datestart            = datestart;
ERD_CUTI_input.dateend              = dateend;
ERD_CUTI_input.num_t                = num_t;
ERD_CUTI_input.dt                   = dt;
ERD_CUTI_input.smoothing_window     = smoothing_window; % set moving average smoothing as time points for averaging before and after current time point
ERD_CUTI_input.target_latitudes  	= target_latitudes; % QQQ chosse only 1 latitude for now; FFF in future, choose 1 or more target latitude(s); [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]

ERD_CUTI                            = f_prep_ERD_CUTI_08052022(ERD_CUTI_input);

AdvecFluxRate                       = ERD_CUTI;        % (m3/s per 1m of BoxWidth); (vertical vector)
fname_UpwellingDriver               = 'ERD_CONST';       % save name of this driver to keep in saved model results
% end case 'ERD_CUTI' -------------------------------------------------
  
  
  case 'Fake_Upwelling' % -----------------------------------------------
  % prepare fake upwelling
maxTREND       = 0.45; % max MLD; (m)
minTREND       = -0.84; % min MLD; (m)
prdTREND       = 365.25/2;
max_dev        = 0.86;
min_dev        = -0.92;
prd_dev        = 13.33; % (15 days) 13.33 gives mean duration of 15.14 days in April-Sept.

rangeTREND     = abs(minTREND - maxTREND); % set minimum mixed layer depth and annual range of the mixed layer depth; (m)
range_dev      = abs(max_dev - min_dev);

% TREND          = (cos(((t_grid*(365/prdTREND)/365)*pi)*1) * (rangeTREND * 0.5) + (maxTREND + rangeTREND * 0.5)); % artificial MLD; (m); (simulates PFEL data for 45N 126W); (vertical vector of length time)
TREND          = -cos(((t_grid*(365/prdTREND)/365)*pi)*1); % artificial MLD; (m); (simulates PFEL data for 45N 126W); (vertical vector of length time)
TREND          = TREND * (rangeTREND * 0.5);
TREND          = TREND + (rangeTREND * 0.5) + minTREND;

dev            = cos(((t_grid*(365/prd_dev)/365)*pi)*1); % oscillate between 1 & -1 over prd of days
dev            = dev * (range_dev * 0.5); % oscillate between +range/2 & -range/2 over prd of days
dev            = dev + (range_dev * 0.5) + min_dev;

Fake_Upwelling = TREND + dev;

AdvecFluxRate           = Fake_Upwelling;       % (m3/s per 1m BoxWidth); (vertical vector: num_t X 1)
fname_UpwellingDriver	= 'Fake_Upwelling';     % save name of this driver to keep in saved model results
if ShowOutput
disp('   -->Fake Upwelling - 15 day'); 
end
% end case 'Fake_Upwelling' -------------------------------------------
  
  end % switch upwelling_driver
% *************************************************************************
  
  
  
  
  
  % *************************************************************************
  % STEP 5: put a limit on AdvecFluxRate and RiverFluxRate-------------------
  %         prevent total washout of small boxes on extreme days
%         NOTE: still possible to have box washed out if RiverFluxRate tips the balance
max_AdvecFluxRate                   = min((BoxVolume * 1/(60*60*24)), [], 2); % max absolute advection rate (m3/s) = one BoxVolume per day; (vertical vector)
looky_WashoutTime                   = find(abs(AdvecFluxRate) > max_AdvecFluxRate);
AdvecFluxRate(looky_WashoutTime)    = sign(AdvecFluxRate(looky_WashoutTime)) .* max_AdvecFluxRate(looky_WashoutTime);
looky_WashoutTime                   = find(abs(RiverFluxRate) > max_AdvecFluxRate);
RiverFluxRate(looky_WashoutTime)    = sign(RiverFluxRate(looky_WashoutTime)) .* max_AdvecFluxRate(looky_WashoutTime);
% *************************************************************************
  
  
  
  
  
  % *************************************************************************
  % STEP 6: define Advective fluxes and River fluxes to use for TimeSeries---
  % step 6a: advection ------------------------------------------------------
  %          NOTE: onshore is POSITIVE, offshore is NEGATIVE
Advection_TimeSeries        = repmat(AdvecFluxRate, [1, 6]) .* repmat([(-1) (1) (-1) (1) (-1) (1)], [num_t, 1]); % (m3/s per 1m BoxWidth); (2D matrix: time X num_AdvecFluxes: A_12 A_13 A_24 A_35 A_4ocean A_5ocean
                                                                                                                                            RiverFlux_TimeSeries        = repmat(RiverFluxRate, [1, 6]) .* repmat([(0)  (0) (-1) (0) (-1) (0)], [num_t, 1]); % (m3/s per 1m BoxWidth); (2D matrix: time X num_AdvecFluxes: A_12 A_13 A_24 A_35 A_4ocean A_5ocean
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        % step 6b: mixing ---------------------------------------------------------
                                                                                                                                                                                                                                                                                          VerticalMixing_TimeSeries   = repmat(we_doy, [1, 2])        .* repmat([1 1], [num_t, 1]); % (m/s); (2D matrix: time X num_MixFluxes) w_23, w_45
                                                                                                                                                                                                                                                                                        %          NOTE: all POSITIVE & already divided 1/2 upward and 1/2 downward
                                                                                                                                                                                                                                                                                        %          NOTE: the second term (repmat([1 1], [num_t, 1])) is not necessary
                                                                                                                                                                                                                                                                                        HorizontalMixing_TimeSeries = ue_doy; % (m3/s per 1m BoxWidth); (2D matrix: time X num_HorizMixFluxes) u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        % *************************************************************************
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          % *************************************************************************
                                                                                                                                                                                                                                                                                          % STEP 7: Physical Structure (NCC upwelling)-------------------------------
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          % step 7a: advection ------------------------------------------------------
                                                                                                                                                                                                                                                                                          %          define advective connections between boxes (3D matrix: destination boxes X source boxes X specific advection flux)
                                                                                                                                                                                                                                                                                        %          NOTE: the sign is used to identify source box depending on onshore-vs-offshore flux
                                                                                                                                                                                                                                                                                        AdvcDirectionSwitch(1:5, 1:5, 1) = [
                                                                                                                                                                                                                                                                                          -1	1	0	0	0
                                                                                                                                                                                                                                                                                          -1	1	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_12 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvcDirectionSwitch(1:5, 1:5, 2) = [
                                                                                                                                                                                                                                                                                          -1	0	1	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          -1	0	1	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_13 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvcDirectionSwitch(1:5, 1:5, 3) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	-1	0	1	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	-1	0	1	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_24 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvcDirectionSwitch(1:5, 1:5, 4) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	-1	0	1
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	-1	0	1]; % A_35 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvcDirectionSwitch(1:5, 1:5, 5) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	-1	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_4ocean (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvcDirectionSwitch(1:5, 1:5, 6) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	-1]; % A_5ocean (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        %  assign proper sign to advective flux (3D matrix: destination boxes X source boxes X specific advection flux)
                                                                                                                                                                                                                                                                                        %          NOTE: AdvectionSign refers to whether a flux is a gain or a loss to a box, it does not refer to a onshore/offshore direction
                                                                                                                                                                                                                                                                                        AdvectionSign(1:5, 1:5, 1) = [
                                                                                                                                                                                                                                                                                          1	1	0	0	0
                                                                                                                                                                                                                                                                                          -1	-1	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_12 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvectionSign(1:5, 1:5, 2) = [
                                                                                                                                                                                                                                                                                          1	0	1	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          -1	0	-1	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_13 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvectionSign(1:5, 1:5, 3) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	1	0	1	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	-1	0	-1	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_24 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvectionSign(1:5, 1:5, 4) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	1	0	1
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	-1	0	-1]; % A_35 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvectionSign(1:5, 1:5, 5) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	1	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % A_4ocean (matrix: spatial cell X spatial cell: I II III IV V); (produces offshore loss to ocean)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        AdvectionSign(1:5, 1:5, 6) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	1]; % A_5ocean (matrix: spatial cell X spatial cell: I II III IV V); (produces offshore loss to ocean)
                                                                                                                                                                                                                                                                                        % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          % step 7b: horizontal mixing ----------------------------------------------
                                                                                                                                                                                                                                                                                          %          define horizontal mixing connections between boxes (3D matrix: DestinationBoxes X SourceSoxes X mixing flux)
                                                                                                                                                                                                                                                                                        HorizontalMixingLink(1:5, 1:5, 1) = [
                                                                                                                                                                                                                                                                                          0	1	0	0	0
                                                                                                                                                                                                                                                                                          1	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % u_12 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        HorizontalMixingLink(1:5, 1:5, 2) = [
                                                                                                                                                                                                                                                                                          0	0	1	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          1	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % u_13 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        HorizontalMixingLink(1:5, 1:5, 3) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	1	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	1	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % u_24 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        HorizontalMixingLink(1:5, 1:5, 4) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	1
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	1	0	0]; % u_35 (matrix: spatial cell X spatial cell: I II III IV V)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        HorizontalMixingLink(1:5, 1:5, 5) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % u_4ocean (matrix: spatial cell X spatial cell: I II III IV V); NOTE: this does NOT define link between box 4 and ocean
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        HorizontalMixingLink(1:5, 1:5, 6) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % u_5ocean (matrix: spatial cell X spatial cell: I II III IV V); NOTE: this does NOT define link between box 5 and ocean
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        % mixing face areas for horizontal mixing fluxes (u_12 u_13 u_24 u_35 u_4ocean u_5ocean)
                                                                                                                                                                                                                                                                                        BoxMixingFace = [BoxHeight_inshore(:, 2) BoxHeight_inshore(:, 3) BoxHeight_inshore(:, 4) BoxHeight_inshore(:, 5) BoxHeight_offshore(:, 4) BoxHeight_offshore(:, 5)] .* repmat(BoxWidth(:, 1), [1, num_HorizMixFluxes]); % (m2); (2D matrix: time X num_HorizMixFluxes); u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          % step 7c: vertical mixing ------------------------------------------------
                                                                                                                                                                                                                                                                                          %          define vertical mixing connections between boxes (3D matrix: DestinationBoxes X SourceSoxes X mixing flux)
                                                                                                                                                                                                                                                                                        VerticalMixingLink(1:5, 1:5, 1) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	1	0	0
                                                                                                                                                                                                                                                                                          0	1	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % w_23
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        VerticalMixingLink(1:5, 1:5, 2) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	1
                                                                                                                                                                                                                                                                                          0	0	0	1	0]; % w_45
                                                                                                                                                                                                                                                                                        % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          % step 7d: sinking --------------------------------------------------------
                                                                                                                                                                                                                                                                                          %          define sinking connections between boxes (3D matrix: DestinationBoxes X SourceSoxes X sinking flux)
                                                                                                                                                                                                                                                                                        SinkLink(1:5, 1:5, 1) = [
                                                                                                                                                                                                                                                                                          -1	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % S_1benthos
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        SinkLink(1:5, 1:5, 2) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	-1	0	0	0
                                                                                                                                                                                                                                                                                          0	1	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % S_23
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        SinkLink(1:5, 1:5, 3) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	-1	0
                                                                                                                                                                                                                                                                                          0	0	0	1	0]; % S_45
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        SinkLink(1:5, 1:5, 4) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	-1	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0]; % S_3benthos
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        SinkLink(1:5, 1:5, 5) = [
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	0
                                                                                                                                                                                                                                                                                          0	0	0	0	-1]; % S_5benthos
                                                                                                                                                                                                                                                                                        % *************************************************************************
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          % *************************************************************************
                                                                                                                                                                                                                                                                                          % STEP 8: net advection & mixing gains/losses to/from each box for each t--
                                                                                                                                                                                                                                                                                          AdvectionMatrix                  = zeros(num_t, num_boxes, num_boxes); % initialize; (3D matrix: time X source X destination (I II III IV V))
                                                                                                                                                                                                                                                                                        RiverFluxMatrix                  = zeros(num_t, num_boxes, num_boxes); % initialize; (3D matrix: time X source X destination (I II III IV V))
                                                                                                                                                                                                                                                                                        VerticalMixMatrix                = zeros(num_t, num_boxes, num_boxes); % initialize; (3D matrix: time X source X destination (I II III IV V))
                                                                                                                                                                                                                                                                                        HorizontalMixMatrix              = zeros(num_t, num_boxes, num_boxes); % initialize; (3D matrix: time X source X destination (I II III IV V))
                                                                                                                                                                                                                                                                                        HorizontalMix_OffshoreSurface    = zeros(num_t, 1);                    % initialize; (vertical vector: time X 1); u_4ocean
                                                                                                                                                                                                                                                                                        HorizontalMix_OffshoreSubSurface = zeros(num_t, 1);                    % initialize; (vertical vector: time X 1); u_5ocean
                                                                                                                                                                                                                                                                                        OceanSourceBoxID                 = zeros(num_t, 1);                    % initialize; (vertical vector)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        for time_loop = 1:num_t % FFF can vectorize all this in future to speed things up
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        current_BoxWidth                     = BoxWidth(time_loop, :);      % (m);  (horizontal vector: 1 X num_boxes)
                                                                                                                                                                                                                                                                                        current_BoxLength                    = BoxLength(time_loop, :);     % (m);  (horizontal vector: 1 X num_boxes)
                                                                                                                                                                                                                                                                                        current_BoxMixingFace                = BoxMixingFace(time_loop, :); % (m2); (horizontal vector: 1 X num_HorizMixFluxes)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        % step 7b: advection --------------------------------------------------
                                                                                                                                                                                                                                                                                          current_advection                    = Advection_TimeSeries(time_loop, :);                          % (m3/s per 1m BoxWidth); (horizontal vector: 1 X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        current_advection                    = reshape(current_advection, 1, 1, num_AdvecFluxes);           % (m3/s per 1m BoxWidth); (vector-along-layers: 1 X 1 X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        advection_repmat                     = repmat(current_advection, [num_boxes, num_boxes, 1]);        % (m3/s per 1m BoxWidth); (3D matrix:  destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        SwitchResult                         = AdvcDirectionSwitch .* advection_repmat;                     % (m3/s per 1m BoxWidth); value is negative where flow direction is wrong @ t; (3D matrix: destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        looky_negative                       = find(SwitchResult <= 0);                                     % find zero or backwards flow entries in matrix
                                                                                                                                                                                                                                                                                        advection_repmat(looky_negative)     = 0;                                                           % set fluxes in source box clms to 0; (m3/s per 1m BoxWidth); (3D matrix:  destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        advection_repmat                     = advection_repmat .* AdvectionSign;                           % put correct sign on direction of flux (offshore is "+"); (m3/s per 1m BoxWidth); (3D matrix: destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        AdvectionMatrix_temp                 = sum(advection_repmat, 3) * (60*60*24);                       % advection rate (m3/d per 1m BoxWidth); (2D matrix: destination boxes X source boxes)
                                                                                                                                                                                                                                                                                        % sum all the individual advection fluxes into and out of each "destination" box,
                                                                                                                                                                                                                                                                                        AdvectionMatrix_temp                 = AdvectionMatrix_temp  .* repmat(current_BoxWidth', [1, num_boxes]); % (m3/d); scale by BoxWidth; (2D matrix: DestinationBox X SourceBox); this is just a formality unless BoxWdith is not 1 m
    AdvectionMatrix(time_loop, :, :)     = reshape(AdvectionMatrix_temp', [1, num_boxes, num_boxes]);   % advection rate (m3/d); (3D matrix:time X SourceBox X DestinationBox (I II III IV V))
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        % step 7c: identify upwelling vs downwelling events by source box ID --
                                                                                                                                                                                                                                                                                          Advection_DestinationSum             = sum(AdvectionMatrix_temp, 2);                                % advected volume (m3/d); (3D matrix: 1 X 1 X DestinationBoxes)
                                                                                                                                                                                                                                                                                        Advection_DestinationSum             = round2(Advection_DestinationSum, 8);                         % round off the very tiny error caused by spline fit
                                                                                                                                                                                                                                                                                        %   NOTE: all advection terms except those from external inputs (i.e. ocean or river) cancel out; 
                                                                                                                                                                                                                                                                                        looky_ReceivingBox = find(Advection_DestinationSum < 0);
                                                                                                                                                                                                                                                                                        if isempty(looky_ReceivingBox)
                                                                                                                                                                                                                                                                                        OceanSourceBoxID(time_loop) = looky_OffshoreSurfaceBox;            % default to downwelling box if all advection = 0
                                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                                          OceanSourceBoxID(time_loop)  = find(Advection_DestinationSum < 0); % box number of box receiving advection from outside
                                                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        % step 7c: riverine advection -----------------------------------------
                                                                                                                                                                                                                                                                                          current_RiverFlux                    = RiverFlux_TimeSeries(time_loop, :);                                  % (m3/s per 1m BoxWidth); (horizontal vector: 1 X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        current_RiverFlux                    = reshape(current_RiverFlux, 1, 1, num_AdvecFluxes);                   % (m3/s per 1m BoxWidth); (vector-along-layers: 1 X 1 X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        RiverFlux_repmat                     = repmat(current_RiverFlux, [num_boxes, num_boxes, 1]);                % (m3/s per 1m BoxWidth); (3D matrix:  destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        SwitchResult                         = AdvcDirectionSwitch .* RiverFlux_repmat;                             % (m3/s per 1m BoxWidth); value is negative where flow direction is wrong @ t; (3D matrix: destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        looky_negative                       = find(SwitchResult <= 0);                                             % find zero or backwards flow entries in matrix
                                                                                                                                                                                                                                                                                        RiverFlux_repmat(looky_negative)     = 0;                                                                   % set fluxes in source box clms to 0; (m3/s per 1m BoxWidth); (3D matrix:  destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        RiverFlux_repmat                     = RiverFlux_repmat .* AdvectionSign;                                   % put correct sign on direction of flux (offshore is "+"); (m3/s per 1m BoxWidth); (3D matrix: destination box X source box X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        RiverFluxMatrix_temp                 = sum(RiverFlux_repmat, 3) * (60*60*24);                               % advection rate (m3/d per 1m BoxWidth); (2D matrix: destination boxes X source boxes)
                                                                                                                                                                                                                                                                                        % sum all the individual advection fluxes into and out of each "destination" box,
                                                                                                                                                                                                                                                                                        RiverFluxMatrix_temp                 = RiverFluxMatrix_temp  .* repmat(current_BoxWidth', [1, num_boxes]);  % scale by BoxWidth;      (2D matrix: DestinationBox X SourceBox);                              (m3/d); this is just a formality unless BoxWdith is not 1 m
    RiverFluxMatrix(time_loop, :, :)     = reshape(RiverFluxMatrix_temp', [1, num_boxes, num_boxes]);           % advection rate (m3/d); (3D matrix:time X SourceBox X DestinationBox (I II III IV V))
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        % step 7d: vertical mixing --------------------------------------------
                                                                                                                                                                                                                                                                                          current_VertMixing                   = VerticalMixing_TimeSeries(time_loop, :);                         % (m/s); (horizontal vector: 1 X num_MixFluxes)
                                                                                                                                                                                                                                                                                        current_VertMixing                   = reshape(current_VertMixing, 1, 1, num_MixFluxes);                % (m/s); (vector-along-layers: 1 X 1 X num_AdvecFluxes)
                                                                                                                                                                                                                                                                                        VertMixing_repmat                    = repmat(current_VertMixing, [num_boxes, num_boxes, 1]);           % (m/s); (3D matrix: num_boxes X num_boxes X num_MixFluxes)
                                                                                                                                                                                                                                                                                        VertMixing_repmat                    = VertMixing_repmat .* VerticalMixingLink;                         % put sign on direction of flux (sign here matters in ODE when multiplied against bottom_box - top_box); (m/s); (3D matrix: num_boxes X num_boxes X num_MixFluxes)
                                                                                                                                                                                                                                                                                        VerticalMixMatrix_temp               = sum(VertMixing_repmat, 3) * (60*60*24);                          % vertical mixing rate (m/d); (2D matrix: destination boxes X source boxes)
                                                                                                                                                                                                                                                                                        VerticalMixMatrix_temp               = VerticalMixMatrix_temp .* repmat((current_BoxLength .* current_BoxWidth)', [1, num_boxes]); % scale by box size (horizontal area); mixing rate (m3/d); (2D matrix: DestinationBoxes X SourceBoxes)
    VerticalMixMatrix(time_loop, :, :)   = reshape(VerticalMixMatrix_temp', [1, num_boxes, num_boxes]);     % mixing rate (m3/d); (3D matrix: time X SourceBox X DestinationBox (I II III IV V))
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        % step 7e: horizontal mixing ------------------------------------------
                                                                                                                                                                                                                                                                                          current_HorizMixing                         = HorizontalMixing_TimeSeries(time_loop, :);                % (m/s); (horizontal vector: 1 X num_HorizMixFluxes); u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        current_HorizMixing                         = current_HorizMixing * (60*60*24);                         % (m/d); (horizontal vector: 1 X num_HorizMixFluxes); u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        current_HorizMixing                         = current_HorizMixing .* current_BoxMixingFace;             % (m3/d); (horizontal vector: 1 X num_HorizMixFluxes); u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        current_HorizMixing                         = reshape(current_HorizMixing, 1, 1, num_HorizMixFluxes);   % (m3/d); (vector-along-layers: 1 X 1 X num_HorizMixFluxes); u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        HorizMixing_repmat                          = repmat(current_HorizMixing, [num_boxes, num_boxes, 1]);   % (m3/d); (3D matrix: DestinationBox X SourceBox X num_HorizMixFluxes); u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        HorizMixing_repmat                          = HorizMixing_repmat .* HorizontalMixingLink;               % (m3/d); (3D matrix: DestinationBox X SourceBox X num_HorizMixFluxes); u_12 u_13 u_24 u_35 u_4ocean u_5ocean
                                                                                                                                                                                                                                                                                        HorizMixing_repmat                          = sum(HorizMixing_repmat, 3);                               % (m3/d); (2D matrix: DestinationBox X SourceBox);
                                                                                                                                                                                                                                                                                        HorizontalMixMatrix(time_loop, :, :)        = reshape(HorizMixing_repmat', [1, num_boxes, num_boxes]);  % horizontal mixing rate; (m3/d); (3D matrix: time X SourceBox X DestinationBox);
    HorizontalMix_OffshoreSurface(time_loop)    = current_HorizMixing(1, 1, 5);                             % horizontal mixing rate u_4ocean; (m3/d); (vertical vector: time X 1)
    HorizontalMix_OffshoreSubSurface(time_loop) = current_HorizMixing(1, 1, 6);                             % horizontal mixing rate u_5ocean; (m3/d); (vertical vector: time X 1)
    
end
% -------------------------------------------------------------------------


% step 8b: sinking fluxes to/from each box --------------------------------
%          (identical for each t; t dimension not defined)
SinkLink                                = sum(SinkLink, 3);                 % flatten SinkLink; (2D matrix: destination boxes X source boxes)
SinkLink_surface                        = zeros(num_boxes);
SinkLink_surface(:, looky_SurfaceBoxes) = SinkLink(:, looky_SurfaceBoxes);  % only consider sinking from surface to sub-surface boxes
SinkLink_benthos                        = SinkLink;
SinkLink_benthos(:, looky_SurfaceBoxes) = 0;                                % only consider sinking from sub-surface to benthos (also includes surface to benthos mixed-boxes (Box I))

SinkLink_surface                 = reshape(SinkLink_surface, [num_boxes, 1, num_boxes]);    % (3D matrix: destination boxes X 1 X source boxes)
SinkLink_benthos                 = reshape(SinkLink_benthos, [num_boxes, 1, num_boxes]);    % (3D matrix: destination boxes X 1 X source boxes)
% -------------------------------------------------------------------------


% step 8c: put flux terms into ECOTRAN-2 format ---------------------------

% ADVECTION --------
ADVECTION           = AdvectionMatrix + RiverFluxMatrix; % advection rate (m3/d); (3D matrix: time X SourceBox X DestinationBox)

ADVECTION(:, 6, 6)	= 0; % add cells for boundary fluxes
export              = sum(ADVECTION, 3); % advection OUT OF each box; (m3/d); (2D matrix: num_t X num_boxes)
import_temp         = sum(ADVECTION, 2); % advection INTO each box; (m3/d); (3D matrix: num_t X 1 X num_boxes)
ADVECTION(:, :, 6)  = export * -1;
ADVECTION(:, 6, :)  = import_temp * -1; % total advection; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))

ADVECTION(:, 1, 1)  = 0; % remove negative box loss terms on the diagonals
ADVECTION(:, 2, 2)  = 0;
ADVECTION(:, 3, 3)  = 0;
ADVECTION(:, 4, 4)  = 0;
ADVECTION(:, 5, 5)  = 0;
ADVECTION(:, 6, 6)  = 0;


% HORIZONTALALMIXING -------
HORIZONTALMIXING            = HorizontalMixMatrix; % horizontal mixing rate (m3/d); (3D matrix: time X SourceBox X DestinationBox)

HORIZONTALMIXING(:, 6, 6)	= 0; % add cells for boundary fluxes; (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
HORIZONTALMIXING(:, 4, 6)	= HorizontalMix_OffshoreSurface;
HORIZONTALMIXING(:, 6, 4)	= HorizontalMix_OffshoreSurface;
HORIZONTALMIXING(:, 5, 6)	= HorizontalMix_OffshoreSubSurface;
HORIZONTALMIXING(:, 6, 5)	= HorizontalMix_OffshoreSubSurface;

HORIZONTALMIXING(:, 1, 1)	= 0; % remove negative box loss terms on the diagonals
HORIZONTALMIXING(:, 2, 2)	= 0;
HORIZONTALMIXING(:, 3, 3)	= 0;
HORIZONTALMIXING(:, 4, 4)	= 0;
HORIZONTALMIXING(:, 5, 5)	= 0;
HORIZONTALMIXING(:, 6, 6)	= 0;


% VERTICALMIXING -------
VERTICALMIXING              = VerticalMixMatrix; % advection rate (m3/d); (3D matrix: time X SourceBox X DestinationBox)

VERTICALMIXING(:, 6, 6)     = 0; % add cells for boundary fluxes (usually 0 for vertical mixing and sinking); (m3/d); (3D matrix: num_t X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))

export                    	= sum(VERTICALMIXING, 3); % (2D matrix: num_t X num_boxes+1 SOURCE)
import_temp                 = sum(VERTICALMIXING, 2); % (3D matrix: num_t X 1 X num_boxes+1 DESTINATION)

% VERTICALMIXING(:, :, 6)    = export * -1;
% VERTICALMIXING(:, 6, :)    = import_temp * -1;


% SINKING -------
SINKING             = reshape(SinkLink', [1, num_boxes, num_boxes]); % (linkage between source to destiny boxes); (3D matrix: 1 X SourceBox X DestinationBox); NOTE transpose
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        SINKING(1, 1, 1)	= 0; % set diagonals to 0
                                                                                                                                                                                                                                                                                        SINKING(1, 2, 2)	= 0;
                                                                                                                                                                                                                                                                                        SINKING(1, 3, 3)	= 0;
                                                                                                                                                                                                                                                                                        SINKING(1, 4, 4)	= 0;
                                                                                                                                                                                                                                                                                        SINKING(1, 5, 5)	= 0;
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        SINKING(1, 6, 6)	= 0; % add cells for boundary fluxes (usually 0 for vertical mixing and sinking); (3D matrix: 1 X SourceBox+1 X DestinationBox+1)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        SINKING(1, 1, 6)	= 1; % sinking from box 1 to benthos; NOTE: DON'T FORGET THIS!
SINKING(1, 3, 6)	= 1; % sinking from box 3 to benthos; NOTE: DON'T FORGET THIS!
                                                                                                                                                                                                                                                                                          SINKING(1, 5, 6)	= 1; % sinking from box 5 to benthos; NOTE: DON'T FORGET THIS!

SINKING             = repmat(SINKING, [num_t, 1, 1]);

BoxFloor            = BoxLength .* BoxWidth;                % (m2); (2D matrix: time X num_boxes SOURCE)
BoxFloor(:, 6)      = 0; % add dimensions (0m2) of boundary import box; (m2); (2D matrix: time X num_boxes SOURCE+1)
repmat_BoxFloor     = repmat(BoxFloor, [1, 1, (num_boxes+1)]);	% (m2); (3D matrix: time X num_boxes SOURCE+1 X num_boxes DESTINY+1)

SINKING             = SINKING .* repmat_BoxFloor;           % sinking source box floor area; (m2); (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
% *************************************************************************





% *************************************************************************
% STEP 9: error checking for conservation of water-volume------------------
%       NOTE: f_CompactFluxTimeSeries_06212019 will be called again externally, these compacted flux terms are not transferred out of this physical function
%       NOTE: no check of SINKING, which is the box floor area and not a flux at this stage
%       CompactFlux
%           compact_flux            (2D matrix: num_t X num_fluxes)
%                                       NOTE: fluxes include all linked boxes +1 for external links
%           looky_flux              (2D matrix: num_fluxes X 5)
%                                       clm 1: (destiny box) = list of boxes importing water volume (+ boundary)
%                                       clm 2: (source box) = list of boxes exporting water volume (+ boundary)
%                                       clm 3: (destiny box address) = addresses of destiny boxes in compact_flux
%                                       clm 4: (source box address) = addresses of source boxes in compact_flux
%                                       clm 5: flux address in non-compacted flux 2D matrix: destiny X source
%                                       NOTE: fluxes include all linked boxes +1 for external links
%                                       NOTE: values constant independent of t
%           looky_boundary_import	(2D matrix: num_fluxes_BoundaryImport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume
%                                       clm 2: (source box) = identity of boxes exporting water volume (always the boundary flux number)
%                                       clm 3: (destiny box address) = addresses of destiny boxes in compact_flux
%           looky_boundary_export   (2D matrix: num_fluxes_BoundaryExport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume (always the boundary flux number)
%                                       clm 2: (source box) = identity of boxes exporting water volume
%                                       clm 3: (source box address) = addresses of source boxes in compact_flux
%           unique_source           (vertical vector: list of source boxes (+ boundary))
%           unique_destiny          (vertical vector: list of destiny boxes (+ boundary))
%           num_fluxes              number of realized fluxes between boxes (or boundary) over full time-series
%       UnCompactFlux
%           flux_import             flux INTO each box (includes external domain fluxes) (2D matrix: num_t X num_boxes)
%           flux_export             flux OUT OF each box (includes external domain fluxes) (2D matrix: num_t X num_boxes)
%           flux_imbalance          imbalance of fluxes for each box (includes and external domain fluxes added into each box-specific flux) (2D matrix: num_t X num_boxes)
%           flux_domain_import      total flux INTO domain from outside (vertical vector: num_t X 1)
%           flux_domain_export      total flux OUT OF domain to outside (vertical vector: num_t X 1)
%           flux_domain_imbalance	total flux imbalance across all domain boundaries (vertical vector: num_t X 1)
%           num_t                   number of time steps (scaler)
%           num_boxes               number of boxes (scaler

% step 9a: check total advection (x + y + z) ------------------------------
CompactFlux     = f_CompactFluxTimeSeries_11182019(ADVECTION,ShowOutput);
UnCompactFlux	= f_UnCompactFluxTimeSeries_12112019(CompactFlux,ShowOutput);
if ShowOutput
display('Testing ADVECTION for flux imbalances...')
end
FluxImbalance	= f_EvaluateFluxBalance_11262021(UnCompactFlux,ShowOutput); % evaluate for any flux imbalance


% step 9b: check horizontal mixing ----------------------------------------
CompactFlux     = f_CompactFluxTimeSeries_11182019(HORIZONTALMIXING,ShowOutput);
UnCompactFlux	= f_UnCompactFluxTimeSeries_12112019(CompactFlux,ShowOutput);
if ShowOutput
display('Testing HORIZONTALMIXING for flux imbalances...')
end
FluxImbalance	= f_EvaluateFluxBalance_11262021(UnCompactFlux,ShowOutput); % evaluate for any flux imbalance


% step 9c: check vertical mixing ------------------------------------------
CompactFlux     = f_CompactFluxTimeSeries_11182019(VERTICALMIXING,ShowOutput);
UnCompactFlux	= f_UnCompactFluxTimeSeries_12112019(CompactFlux,ShowOutput);
if ShowOutput
display('Testing VERTICALMIXING for flux imbalances...')
end
FluxImbalance	= f_EvaluateFluxBalance_11262021(UnCompactFlux,ShowOutput); % evaluate for any flux imbalance
% *************************************************************************





% *************************************************************************
% STEP 10: light parameters------------------------------------------------
LightParams.latitude        = latitude;   % latitude of Deep-Water Horizon oil spill (decimal degress north)
LightParams.SolarConstant	= 1366;       % Solar const. (W m^-2) recent ave. (Froehlich & Lean (1998, Geographical Research Letters 25(23): 4377-4380))

LightParams.par_frac        = 0.43;       % fraction of 40 that is PAR, clear atmos. (Fasham, Ducklow, & McKelvie (1990, J. Mar. Res. 48:591-639), p. 604)
LightParams.surf_trans      = 0.85;       % 1-reflectance (Parsons, Takahachi, & Hargrave (1984, Biological Oceanographic Processes, 3rd ed.), p. 69)
LightParams.cloud_cover     = 0.50;       % fitting daily average and weekly-smoothed light model curve to weekly average of HMSC roof data by eye
% LightParams.cloud_cover     = 0.70;       % Ave. 1990-97, COADS Standard 1-deg. for 125.5W, 45.5N
LightParams.cloud_trans     = 1.0 - (0.7 * LightParams.cloud_cover);  % (Evans and Parslow (1985, Biol. Oceanogr. 3(3):327-347), p. 344)
LightParams.t_grid          = t_grid;

SurfaceRadiation            = f_LightIntensity_12112020(LightParams);	% instantaneous (W/m2), daily mean averaged across 24 h (W m^-2 h^-1), & daily integrated (W m^-2 d^-1) solar raditation at ocean surface
daily_average_light         = SurfaceRadiation.daily_average_light;     % daily integrated solar raditation at ocean surface averaged across 24 hours (W m^-2 h^-1); (vertical vector: num_t X 1)
%                                                 NOTE: the daily average is also what was used by Spitz in her NPZD models
%                                                 NOTE: the daily integrated value avergaed over 24 hours is almost identical to the mean of the current instantaneous light intensity averaged over 24 hours (I checked)
daily_integrated_light      = SurfaceRadiation.daily_integrated_light;	% time-series of daily integrated solar raditation at ocean surface (W m^-2 d^-1)<--QQQ check units, integrated value should not be a rate but a TOTAL-->(W m^-2)

current_light               = SurfaceRadiation.current_light;           % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day
sunrise                     = SurfaceRadiation.sunrise;                  % time of sunrise (time in h from midnight; 12 = noon, set as a default)
sunset                      = SurfaceRadiation.sunset;                   % time of sunset  (time in h from midnight; 12 = noon, set as a default)

Kw                          = 0.067;    % Light attenuation_seawater (Newberger et al., 2003); (1/m)
Kp                          = 0.0095;   % Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N)
% *************************************************************************





% *************************************************************************
% STEP 11: temperature parameters------------------------------------------
temperature_timeseries      = ones(num_t, num_boxes) * 12; % default place-holder temperature time-series; (deg C); (2D matrix: num_t X num_boxes)
temperature_reference       = ones(1, num_boxes) * 12; % default place-holder reference temperatures; (deg C); (horizontal vector: 1 X num_boxes)
temperature_dates           = (datestart:dt:dateend)'; % (vertical vector: num_t X 1); QQQ has not been tested with dt ~= 1;
                                                                                                                                                                                                                                                                                        % *************************************************************************
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                          % *************************************************************************
                                                                                                                                                                                                                                                                                          % STEP 12: calcluate initial NO3 & NH4 input rates--------------------------
                                                                                                                                                                                                                                                                                          %          NOTE: using summer 150-200 of physical flux as initial flux rates
                                                                                                                                                                                                                                                                                        %          NOTE: using annual mean NO3 concentrations as initial boundary conditions
                                                                                                                                                                                                                                                                                        initial_doy                 = 150:200; % days over which to average initial rates
                                                                                                                                                                                                                                                                                        initial_BoxVolume           = mean(BoxVolume(initial_doy, :), 1); % initial BoxVolume averaged over first 50 time-steps; (m3); (2D matrix: 1 X num_boxes)
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        initial_ADVECTION           = mean(ADVECTION(initial_doy, :, :), 1); % initial advection averaged over first 50 time-steps; (m3/d); (3D matrix: 1 X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
                                                                                                                                                                                                                                                                                        initial_HORIZONTALMIXING	= mean(HORIZONTALMIXING(initial_doy, :, :), 1); % initial advection averaged over first 50 time-steps; (m3/d); (3D matrix: 1 X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
                                                                                                                                                                                                                                                                                        initial_VERTICALMIXING      = mean(VERTICALMIXING(initial_doy, :, :), 1); % initial advection averaged over first 50 time-steps; (m3/d); (3D matrix: 1 X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                        meanNO3_annual              = mean(dat_NO3, 1); % mean annual NO3 concentration; (mmole N/m3); (horizontal vector: 1 X num_boxes)
                                                                                                                                                                                                                                                                                        meanNO3_doy                	= mean(NO3_doy(initial_doy, :), 1); % (mmoles N/m3); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                        initial_NO3_doy             = meanNO3_doy; % (mmoles N/m3); (horizontal vector: 1 X num_boxes)
                                                                                                                                                                                                                                                                                        initial_NO3_doy             = repmat(initial_NO3_doy', [1, (num_boxes+1)]); % (mmoles N/m3); (2D matrix: num_boxes X num_boxes); NOTE transpose
initial_NO3_doy((end+1), 1:(end-1)) = meanNO3_doy; % add reflective external SOURCE boundary conditions; (mmoles N/m3); (2D matrix: SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
initial_NO3_doy             = reshape(initial_NO3_doy, [1, (num_boxes+1), (num_boxes+1)]); % (mmoles N/m3); (3D matrix: 1 X SOURCE (num_boxes+1) X DESTINY (num_boxes+1)); NOTE transpose

NO3initial_rate             = (initial_NO3_doy .* initial_ADVECTION) + (initial_NO3_doy .* initial_HORIZONTALMIXING) + (initial_NO3_doy .* initial_VERTICALMIXING); % (mmoles N/d); (3D matrix: 1 X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
NO3initial_rate             = squeeze(sum(NO3initial_rate, 2)); % (mmoles N/d); (vertical vector: DESTINY (num_boxes+1) X 1)
NO3initial_rate             = NO3initial_rate(1:(end-1))' ./ initial_BoxVolume; % initial NO3 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))
                                                                                                                                                                                                                                                                                                                             % *************************************************************************
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % *************************************************************************
                                                                                                                                                                                                                                                                                                                               % STEP 13: pack results for export-----------------------------------------
                                                                                                                                                                                                                                                                                                                               % step 13a: save names of this m-file and sub-functions -------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.fname_PhysicalModel              = fname_PhysicalModel;                      % name of this sub-function
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.fname_UpwellingDriver        	= fname_UpwellingDriver;                    % name of the upwelling driver technique
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.fname_LightIntensity             = SurfaceRadiation.fname_LightIntensity;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.fname_CompactFlux                = CompactFlux.fname_CompactFlux;            % file name of f_CompactFluxTimeSeries function
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.fname_EvaluateFluxBalance        = FluxImbalance.fname_EvaluateFluxBalance;	% name of this f_EvaluateFluxBalance sub-function
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.fname_UnCompactFluxTimeSeries	= UnCompactFlux.fname_UnCompactFluxTimeSeries; % name of this f_UnCompactFluxTimeSeries sub-function
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.fname_CalcNetFlux                = UnCompactFlux.fname_CalcNetFlux;          % name of this f_CalcNetFlux function
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13b: time dimensions -----------------------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.datestart                        = datestart;	% starting date
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.dateend                          = dateend;      % ending date
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.t_grid                           = t_grid;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.num_t                            = num_t;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.dt                               = dt;           % time-step (d)
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13c: space dimensions ----------------------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.num_boxes                    = num_boxes;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxLength                    = BoxLength;                 % (m); (2D matrix: time X num_boxes)
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.BoxLength_y               	  = BoxLength_y; 	% North-South (y) orientation; (m); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.BoxWidth_x               	  = BoxWidth_x;  	% East-West (x) orientation; (m); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.BoxHeight                    = BoxHeight;   	% mean box height over time; (m); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxHeight_inshore            = BoxHeight_inshore;         % I, II, III, IV, V; (m); (2D matrix: time X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxHeight_offshore           = BoxHeight_offshore;        % I, II, III, IV, V; (m); (2D matrix: time X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxHeight                    = BoxHeight;                 % I, II, III, IV, V; mean box height; (m); (2D matrix: time X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxWidth                     = BoxWidth;                  % I, II, III, IV, V; (m); (2D matrix: time X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxArea                      = BoxArea;                   % I, II, III, IV, V; (m2); (2D matrix: time X num_boxes)
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.BoxArea_y                 	  = BoxArea_y;  	% side or face area; North-South (y) orientation; (m2); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.BoxArea_x                 	  = BoxArea_x;   	% side or face area; East-West (x) orientation; (m2); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxFloorArea                 = BoxFloorArea;             % (m2); (2D matrix: time X num_boxes SOURCE)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.BoxVolume                    = BoxVolume;                 % I, II, III, IV, V; (m3); (2D matrix: time X num_boxes)
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % % % % step 13d: spatial relationships -----------------------------------------
                                                                                                                                                                                                                                                                                                                               % % % ECOTRANphysics.advection_links_x            = advection_links_x(1, :, :); % (3D matrix: 1 X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.advection_links_y            = advection_links_y(1, :, :); % (3D matrix: 1 X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.advection_links_z            = advection_links_z(1, :, :); % (3D matrix: 1 X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.HorizMixing_links            = HorizMixing_links(1, :, :); % QQQ need to differentiate x & y; (3D matrix: 1 X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.VertMixing_links             = VertMixing_links(1, :, :);  % (3D matrix: 1 X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.sink_links                   = sink_links(1, :, :);        % (3D matrix: 1 X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             % % % % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13e: physical driver time-series -----------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.AdvecFluxRate                = AdvecFluxRate;    % upwelling index time-series;  (m3/s per 1m BoxWidth); (vertical vector: num_t X 1)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.ADVECTION                    = ADVECTION;        % advection rate                (m3/d); (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.HORIZONTALMIXING             = HORIZONTALMIXING;	% horizontal mixing rate        (m3/d); (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.VERTICALMIXING               = VERTICALMIXING;	% vertical mixing rate          (m3/d); (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.SINKING                   	= SINKING;          % sinking source box floor area (m2);   (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13f: box addresses -------------------------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.looky_SubSurfaceBoxes       = looky_SubSurfaceBoxes;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.looky_SurfaceBoxes          = looky_SurfaceBoxes;
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.looky_OffshoreSurfaceBox    = looky_OffshoreSurfaceBox;
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.looky_OffshoreSubSurfaceBox = looky_OffshoreSubSurfaceBox;
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.looky_RiverFluxBox          = looky_RiverFluxBox;        % Optional river alongshore input; (Box 2 in CGoA downwelling model)
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.OceanSourceBoxID            = OceanSourceBoxID;          % ocean sourcebox ID; (vertical vector: time X 1)
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % % step 13g: raw physical driver time-series (only for NWPO3_BUI or Brink_BUI cases)-------------------------------
                                                                                                                                                                                                                                                                                                                               % ECOTRANphysics.tau                          = tau;      % Alongshore wind stress, monthly mean: (N/m^2) = (kg*m/s2/m2) = (kg/s2/m); (horizontal vector: 1 X month)
                                                                                                                                                                                                                                                                                                                             % ECOTRANphysics.tau_doy                      = tau_doy;	% Alongshore wind stress, monthly mean: (N/m^2) = (kg*m/m2/s2) = (kg/m/s2); (vertical vector: num_t X 1)
                                                                                                                                                                                                                                                                                                                             % ECOTRANphysics.we                           = we;       % vertical mixing rate,   monthly mean: (m/s); (2D matrix: num_VertMixFluxes X month)
                                                                                                                                                                                                                                                                                                                             % ECOTRANphysics.we_doy                       = we_doy;	% vertical mixing rate,   monthly mean: (m/s); (2D matrix: num_t X num_VertMixFluxes)
                                                                                                                                                                                                                                                                                                                             % ECOTRANphysics.ue                           = ue;       % horizontal mixing rate, monthly mean: (m/s); (2D matrix: num_HorizMixFluxes X month)
                                                                                                                                                                                                                                                                                                                             % ECOTRANphysics.ue_doy                       = ue_doy;	% horizontal mixing rate, monthly mean: (m/s); (2D matrix: num_t X num_HorizMixFluxes)
                                                                                                                                                                                                                                                                                                                             % % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13h: misc info -----------------------------------------------------
                                                                                                                                                                                                                                                                                                                               % ECOTRANphysics.SinkLink_surface            = SinkLink_surface;          % define sinking connections from surface to sub-surface boxes; (unitless); (2D matrix: destination boxes X source boxes)
                                                                                                                                                                                                                                                                                                                             % ECOTRANphysics.SinkLink_benthos            = SinkLink_benthos;          % define sinking connections from sub-surface boxes to benthic detritus; (unitless); (2D matrix: destination boxes X source boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.SinkLink_surface            = SinkLink_surface;          % define sinking connections from surface to sub-surface boxes; (unitless); (3D matrix: destination boxes X 1 X source boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.SinkLink_benthos            = SinkLink_benthos;          % define sinking connections from sub-surface boxes to benthic detritus; (unitless); (3D matrix: destination boxes X 1 X source boxes)
                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.num_AdvecFluxes             = num_AdvecFluxes;
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.num_MixFluxes               = num_MixFluxes;
                                                                                                                                                                                                                                                                                                                             % % % ECOTRANphysics.num_HorizMixFluxes          = num_HorizMixFluxes;
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13i: nutrient time-series ------------------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.dat_NO3                      = dat_NO3;              	% monthly mean NO3 + NO2 concentration; (mmole N/m3); (2D matrix: month (12) X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.meanNO3_annual               = meanNO3_annual;           % mean annual NO3 concentration; (mmole N/m3); (horizontal vector: 1 X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.NO3timeseries_conc           = NO3_doy;              	% daily NO3 + NO2 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.NO3initial_rate              = NO3initial_rate;          % initial NO3 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13j: light model info ----------------------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.MLD                          = MLD;                   	% time-series of Mixed Layer Depth (for light intensity calculations, NOT advection calcs); (m); (vertical vector: num_t X 1)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.Io                           = daily_average_light;    	% time-series of surface PAR light intensity; daily mean solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1) (Brock, 1981, R. Ozretich unpub. data); (vertical vector: num_t X 1)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.daily_integrated_light       = daily_integrated_light;	% time-series of daily integrated solar raditation at ocean surface (W m^-2 d^-1)<--QQQ check units, integrated value should not be a rate but a TOTAL-->(W m^-2); (vertical vector: num_t X 1)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.current_light                = current_light;            % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day; (vertical vector: num_t X 1)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.sunrise                      = sunrise;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.sunset                       = sunset;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.Kw                           = Kw;                       % Light attenuation of seawater; (scalar)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.Kp                           = Kp;                     	% Light attenuation of phytoplankton (Newberger et al., 2003); (m2/mmol N); (scalar)
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13k: temperature time-series ---------------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.temperature_timeseries       = temperature_timeseries;	% (deg C); (2D matrix: num_t X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.temperature_reference        = temperature_reference;    % (deg C); (horizontal vector: 1 X num_boxes)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.temperature_dates            = temperature_dates;      	% (vertical vector: num_t X 1)
                                                                                                                                                                                                                                                                                                                             % -------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % step 13L: conversion factors --------------------------------------------
                                                                                                                                                                                                                                                                                                                               ECOTRANphysics.atomic_mass_C               = atomic_mass_C;             % (g C / mole C)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.atomic_mass_N               = atomic_mass_N;             % (g N / mole N)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.C_to_N_phytoplankton        = C_to_N_phytoplankton;      % Redfield ratio = C:N:P = 106:16:1 ---> C:N = 106:16 = 6.625 mole C/mole N
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.C_to_N_zooplankton          = C_to_N_zooplankton;        % general zooplankton C:N (moles C moles N^-1) for mesozooplankton (W. Peterson pers com); C:N=5.16-5.26 [Schneider 1990 (Mar Biol 106:219-225)]; C:N=3.9 [Uye 1982]
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.Chla_to_N_diatoms           = Chla_to_N_diatoms;         % (mg Chla mmoles N^-1) (Dickson & Wheeler 1995 L&O 40(3):533-543)
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.Chla_to_N_dinoflagellates   = Chla_to_N_dinoflagellates;
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.mgWWT_per_mmolesN           = mgWWT_per_mmolesN;         % (mg WWT/mmoles N); for phytoplankton
                                                                                                                                                                                                                                                                                                                             ECOTRANphysics.WWT_to_C                    = WWT_to_C;                  % (g WWT / g C); (Steele et al. 2007); QQQ recheck this for phytoplankton
                                                                                                                                                                                                                                                                                                                             % *************************************************************************
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                               % end m-file***************************************************************