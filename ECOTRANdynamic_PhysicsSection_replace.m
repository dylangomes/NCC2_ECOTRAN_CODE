%% *************************************************************************
% STEP 6: prepare physical parameters--------------------------------------

% step 6a: Define time-frame and time-step info
%          NOTE: prepared BUI_ERD .csv file covers: 01-Jan-1969 through 30-Sep-2015 and has 6-hr resolution
%          NOTE: ERDcuti file covers:               28-Feb-1990 through 28-Feb-2022  and has 1-day resolution
%          NOTE: NH Line temperature covers:        05-May-1997 through 08-Jul-2021 and has 1-day resolution

switch switch_PhysicalModel
    case '2D_upwelling'
        
        disp('--->>> 2D upwelling physics')

        datestart                       = datenum('01-Jan-1998'); % SSS --> enter starting date
        dateend                         = datenum('31-Dec-2022'); % SSS --> enter ending date
        
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
        PHYSICSinput.smoothing_window	= 5; % window for smoothing before and after time-point (e.g., 2 = a window of 5 days centered on time-point t)
        
        spatial_BiomassScalers         	= [1 1 0 1 0];	% NCC scalers for estimating initial (or mean) primary producer biomasses across model domain

        % -------------------------------------------------------------------------


        % step 6b: prepare advection & mixing time-series for each model box ------
        % 2D upwelling driver
        ECOTRANphysics                  = f_ECOTRANphysics_upwelling_09052022(PHYSICSinput, 'ERD_CUTI'); % SSS specify physical flux time-series to use: 'Brink_BUI', 'NWPO3_BUI', 'ERD_BUI', 'ERD_CUTI', or 'Fake_Upwelling'

        % apply any desired change to the temperature time-series (DO NOT change temperature_reference)
        ECOTRANphysics.temperature_timeseries	= ECOTRANphysics.temperature_timeseries + temperature_scenario; % (2D matrix: num_t X num_boxes)
        % -------------------------------------------------------------------------


        % step 6c: migration of each group (shared box face areas) ----------------
        %          NOTE: code does not handle migration flux uncertainty (nor physical flux uncertainty)
        ECOTRANmigration            = f_ECOTRANmigration_NCC_2D_09042022(ECOTRANphysics);

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

        
        % Special case for 2D cross-shelf upwelling physics
        ODEinput.SinkLink_surface = ECOTRANphysics.SinkLink_surface; % only used for 5-box, 2D cross-shelf settings; (3D matrix: destination boxes X 1 X source boxes)
        
	% end (case '2D_upwelling') --------------

    
    case '3D_ROMS'
        
        disp('--->>> 3D ROMS physics')

        datestart                       = datenum('01-Jan-2005'); % SSS --> enter starting date
        dateend                         = datenum('31-Dec-2006'); % SSS --> enter ending date (default for dynamic runs tests ('31-Dec-2020'))
        dt                              = 24/24; % t-step; (days); (dt = 24/24 = 1 d; dt = 3/24 = 3 hours)
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
        PHYSICSinput.t_grid_real     	= linspace(datestart, (dateend+1-dt), ((dateend - datestart + 1)/dt))'; % QQQ NEW VERSION!!!; (vertical vector); t_grid runs from datestart to dateend and is inclusive of dateend; intervals = dt
        % -------------------------------------------------------------------------


        % step 6b: prepare advection & mixing time-series for each model box ------
        % 3D ROMS driver
        ECOTRANphysics                  = f_ECOTRANphysics_NCC_ROMS_12022022(PHYSICSinput); % use for NCC ROMS

        % p_plotECOTRANphysics_06242019(ECOTRANphysics)
        % -------------------------------------------------------------------------


        % step 6c: migration of each group (shared box face areas) ----------------
        %          NOTE: code does not handle migration flux uncertainty (nor physical flux uncertainty)
        ECOTRANmigration                = f_ECOTRANmigration_NCC_02222022(ECOTRANphysics);  % SSS; use for NCC ROMS

        % ODEinput.biomass_migrator       = ECOTRANmigration.biomass_migrator;  % SSS special definition of boundary biomasses for migrators; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes); NOTE: de-comment migrator biomass lines within ODE code
        % -------------------------------------------------------------------------

        
        CompactFlux_ADVECTION           = ECOTRANphysics.CompactFlux_ADVECTION; % (structure)
        CompactFlux_HORIZONTALMIXING	= ECOTRANphysics.CompactFlux_HORIZONTALMIXING; % (structure)
        CompactFlux_VERTICALMIXING      = ECOTRANphysics.CompactFlux_VERTICALMIXING; % (structure)
        CompactFlux_SINKING             = ECOTRANphysics.CompactFlux_SINKING; % (structure); compact SINKING as box floor areas and connectivity information; apply functional group sinking speeds in ECOTRANdynamic_ code
        CompactFlux_MIGRATION           = ECOTRANmigration.CompactFlux_MIGRATION; % (structure)

        % QQQ these 3 lines should come out of the "case" switch; make same for 2D mode
        num_boxes                       = ECOTRANphysics.num_boxes;
        num_domains                     = ECOTRANphysics.num_domains;	% number of geographic domains (does not consider depth layers)
        num_z                           = ECOTRANphysics.num_z;         % number of depth layers

        spatial_BiomassScalers         	= ones(1, num_boxes);	% NCC scalers for estimating initial (or mean) primary producer biomasses across model domain; NOTE: these values are assumed; NOTE: x2 in Box I used to compensate for 30m depth relative to 15 m depths in Boxes II & IV; FFF apply NPZD scalers here
        
        % ROMS BioGeoChemical (BGC) model info
        ROMS_temperature_initial     	= ECOTRANphysics.ROMS_temperature_initial;          % (deg C); (horizontal vector: 1 X num_boxes)
        ROMS_diatom                     = ECOTRANphysics.diatom_timeseries;                 % (mmole N/m3); (2D matrix: num_t X num_boxes)
        ROMS_nanophytoplankton          = ECOTRANphysics.nanophytoplankton_timeseries;      % (mmole N/m3); (2D matrix: num_t X num_boxes)
        ROMS_diatom_initial          	= ECOTRANphysics.ROMS_diatom_initial;               % (mmoles N/m3); (horizontal vector: 1 X num_boxes)
        ROMS_nanophytoplankton_initial	= ECOTRANphysics.ROMS_nanophytoplankton_initial;	% (mmoles N/m3); (horizontal vector: 1 X num_boxes)
        
	% end (case '3D_ROMS') --------------

end % (switch switch_PhysicalModel) ---------------------------------------
% -------------------------------------------------------------------------


% step 6d: unpack & process physics variables -----------------------------
num_boxes                    	= ECOTRANphysics.num_boxes;
BoxVolume                    	= ECOTRANphysics.BoxVolume;     	% (m3); (2D matrix: num_t X num_boxes)
% BoxLength                    	= ECOTRANphysics.BoxLength;        	% (m); (2D matrix: num_t X num_boxes)
% BoxHeight                     	= ECOTRANphysics.BoxHeight;       	% (m); (2D matrix: num_t X num_boxes)
% BoxWidth                     	= ECOTRANphysics.BoxWidth;       	% (m); (2D matrix: num_t X num_boxes)

Io                            	= ECOTRANphysics.Io;             	% time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (vertical vector: num_t X 1)
current_light                   = ECOTRANphysics.current_light;     % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day; (vertical vector: num_t X 1)
sunrise                         = ECOTRANphysics.sunrise;           % time of sunrise (time in h from midnight; 12 = noon, set as a default)
sunset                          = ECOTRANphysics.sunset;            % time of sunset  (time in h from midnight; 12 = noon, set as a default)
Kw                            	= ECOTRANphysics.Kw;              	% Light attenuation_seawater; (scalar)
Kp                             	= ECOTRANphysics.Kp;              	% Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N); (scalar)
MLD                           	= ECOTRANphysics.MLD;            	% mixed-layer depth; (m); (vertical vector: num_t X 1)
EuphoticDepth                  	= repmat(MLD, [1 num_boxes]);     	% FFF (eventually move to f_physics code); depth of euphotic zone, used when converting the vertically-integrated EwE primary producer biomass to biomass/volume; depth; (m); (2D matrix: num_t X num_boxes)

% nutrient info for drivers (e.g., NH Line nutrient drivers)
NO3timeseries_conc              = ECOTRANphysics.NO3timeseries_conc;    % NO3 + NO2 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
% NH4timeseries_conc              = ECOTRANphysics.NH4timeseries_conc;    % NH4 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
NO3initial_rate                 = ECOTRANphysics.NO3initial_rate;       % initial NO3 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))
% NH4initial_rate                 = ECOTRANphysics.NH4initial_rate;       % initial NH4 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))

WWT_to_C                      	= ECOTRANphysics.WWT_to_C;              % (scalar)
atomic_mass_C                  	= ECOTRANphysics.atomic_mass_C;         % (scalar)
C_to_N_phytoplankton          	= ECOTRANphysics.C_to_N_phytoplankton;  % (scalar)

grp_row                      	= 1:num_grps;
% -------------------------------------------------------------------------


%% step 6e: compact fluxes -------------------------------------------------
%          remove information defining non-existing box links
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
SinkingSpeed                                = zeros(num_t, num_grps, num_fluxes_sinking);	% initialze sinking speed time-series; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)

SinkingSpeed(:, rc_plgc_detritus, :)      	= repmat((10.5), [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% SinkingSpeed(:, rc_bnth_detritus, :)       	= repmat((25),   [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% SinkingSpeed(:, rc_fishery_offal, :)       	= repmat((40),   [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
disp('NOTICE: for vertically-integrated food-web technique IN 2D & 3D ROMS SHELF SETTINGS, sinking is applied ONLY FOR PELAGIC DETRITUS (change manually in code)')

% MichaelisMenten_w                               = [0.6 1.0];                                    % SSS; sinking speed; (m/d); [Sm Phytoplankton, Lg Phytoplankton]
% SinkingSpeed(:, looky_ANYPrimaryProducer, :)	= repmat(MichaelisMenten_w, [num_t, 1, num_fluxes_sinking]);	% sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% -------------------------------------------------------------------------


%% step 6g: define DIEL VERTICAL MIGRATION (DVM) speeds & duration for each group (m/d) --------
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


%% step 6h: define production loss fractions -------------------------------
%          NOTE: this used for initial conditions only (for dynamic models)
%          NOTE: the ODE accounts for physical removal (or addition) at each time-step; while the static model accounts for physical loss as a reduction in Transfer Efficiency
%                this is the way it was handled in John's original code outline and in the pubs
PhysicalLossFraction        = repmat(ProductionLossScaler', [num_t, 1]);	% (2D matrix: num_t X num_grps); NOTE transpose
PhysicalLossFraction        = PhysicalLossFraction * 0;                     % SSS set to 0 for time-dynamic runs
% -------------------------------------------------------------------------


%% step 6i: pack physics values for ODE into ODEinput ----------------------
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


