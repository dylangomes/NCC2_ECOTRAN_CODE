function [Fluxes_t] = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput)
% calculate advection, mixing, & sinking exchange rates between boxes and out of domain
% by Jim Ruzicka
%         NOTE: flux_domain_import_t value = 0 for driver group, driver group domain input is given in flux_domain_import_driver_t
%         NOTE: flux_net_t is sum of all fluxes INTO each box (including domain import) MINUS sum of all fluxes OUT OF each box (including domain export)
%         NOTE: flux_domain_export_t values are all positive even though this is a loss from the box
%         NOTE: code currently only works with 1 defined external driver group (e.g. NO3) (FFF could chnage this in the future)
% takes:
%	intraODEinput
%       flux_domain_import_t	initialized variable; (2D matrix: num_grps X num_boxes+1)
%       flux_domain_export_t	initialized variable; (2D matrix: num_grps X num_boxes+1)
%       flux_domain_import_driver_t	initialized variable; (2D matrix: num_grps X num_boxes DESTINATION) NOTE: passed through ODE as ODEinput.flux_domain_import_driver_t and is not the same as external_driver_t interpolated within the ODE
%       repmat_BoxVolume      	(m3); (2D matrix: num_grps X num_boxes)
%       compact_flux_t          volume transported per day @ t; (m3/d); (horizontal vector: 1 X num_fluxes)
%                                   NOTE: fluxes include all linked boxes +1 for external links
%       looky_flux              (2D matrix: num_fluxes X 3)
%                                   clm 1: (destiny box) = identity of boxes importing water volume (+ boundary)
%                                   clm 2: (source box) = identity of boxes exporting water volume (+ boundary)
%                                   clm 3: flux address in non-compacted flux 2D matrix: destiny X source
%                                   NOTE: fluxes include all linked boxes +1 for external links
%                                   NOTE: values constant independent of t
%       looky_boundary_import	(2D matrix: num_fluxes_BoundaryImport X 3)
%                                   clm 1: (destiny box) = identity of boxes importing water volume
%                                   clm 2: (source box) = identity of boxes exporting water volume (always the boundary flux box number)
%                                   clm 3: (import flux address) = addresses of import flux clm in compact_flux
%       looky_boundary_export   (2D matrix: num_fluxes_BoundaryExport X 3)
%                                   clm 1: (destiny box) = identity of boxes importing water volume (always the boundary flux box number)
%                                   clm 2: (source box) = identity of boxes exporting water volume
%                                   clm 3: (export flux address) = addresses of export flux clm in compact_flux
%       repmat_looky_source     source box for each flux in compact_flux_t; clm 2 of looky_flux; (2D matrix: num_grps X num_fluxes)
%       repmat_looky_destiny	destiny box for each flux in compact_flux_t; clm 1 of looky_flux; (2D matrix: num_grps X num_fluxes)
%       repmat_grp_row          addressses; list of grp numbers; (2D matrix: num_grps X num_fluxes); NOTE: num_fluxes is specific to physical flux type (e.g., advection, mixing, or sinking)
%       biomass_plus1       	biomass density @ t; (mmoles N/m3)(3D matrix: num_grps X 1 X num_boxes+1) NOTE: sinking flux units (mmoles N/m2/d)
%       biomass_boundary    	biomass conditions @ t; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: equals biomass_plus1 for reflective boundary, but may be different when using defined boundary conditions
%       looky_driver          	row address(es) of driver group(s) (e.g., NO3)
%       RetentionScaler         (scaler 0-1); (2D matrix: num_grps X num_boxes+1 DESTINATION)
%	ODEinput
%       num_grps
%       num_boxes
% returns:
%   Fluxes_t
%       flux_import_t               flux INTO each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does NOT include domain import (see flux_domain_import_t for domain input rates)
%       flux_export_t               flux OUT OF each box; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION); NOTE: does include domain export
%       flux_net_t                  net biomass flux into(+) or out of(-) each box; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: does include domain import & export
%       flux_domain_import_t        net flux across domain boundary INTO each box (driver grp = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: driver group set to 0 (see flux_domain_import_driver_t for driver import rate)
%       flux_domain_export_t        net flux across domain boundary OUT OF each box; (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION)
%       flux_domain_import_driver_t	net flux OF DRIVER GRP across domain boundary INTO each box (all other grps = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION)
% revision date: 9-9-2019


% *************************************************************************
% STEP 1: unpack variables-------------------------------------------------
compact_flux_t              = intraODEinput.compact_flux_t;         % volume transported per day; (m3/d); (2D matrix: num_grps X num_fluxes); NOTE: sinking flux units (m2)
flux_domain_import_t        = intraODEinput.flux_domain_import_t;	% initialized as all zeros; (2D matrix: num_grps X num_boxes+1)
flux_domain_export_t        = intraODEinput.flux_domain_export_t;	% initialized as all zeros; (2D matrix: num_grps X num_boxes+1)
flux_domain_import_driver_t	= intraODEinput.flux_domain_import_driver_t; % initialized as all zeros; (2D matrix: num_grps X num_boxes DESTINATION)
biomass_plus1               = intraODEinput.biomass_plus1;          % (mmoles N/m3); NOTE: sinking flux units (mmoles N/m2/d); (3D matrix: num_grps X 1 X num_boxes+1); FFF can I keep as 2D matrix?
biomass_boundary_t          = intraODEinput.biomass_boundary;       % boundary biomass conditions @ t; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: equals biomass_plus1 for reflective boundary, but may be different when using defined boundary conditions
RetentionScaler             = intraODEinput.RetentionScaler;        % (scaler 0-1); (2D matrix: num_grps X num_boxes+1 DESTINATION)
looky_flux                  = intraODEinput.looky_flux;             % (2D matrix: num_fluxes X 5-->[(DESTINY box) (SOURCE box) (DESTINY box address in compact_flux) (SOURCE box address in compact_flux) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
looky_boundary_import       = intraODEinput.looky_boundary_import;	% (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
looky_boundary_export       = intraODEinput.looky_boundary_export;  % (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
looky_driver                = intraODEinput.looky_driver;           % row address(es) of driver group(s) (e.g., NO3)
repmat_looky_source         = intraODEinput.repmat_looky_source;	% source box for each flux in compact_flux_t; clm 2 of looky_flux; (2D matrix: num_grps X num_fluxes)
repmat_looky_destiny        = intraODEinput.repmat_looky_destiny;   % destiny box for each flux in compact_flux_t; clm 1 of looky_flux; (2D matrix: num_grps X num_fluxes)
repmat_grp_row              = intraODEinput.repmat_grp_row;         % addressses; list of grp numbers; (2D matrix: num_grps X num_fluxes); NOTE: num_fluxes is specific to physical flux type (e.g., advection, mixing, sinking, or migration)
repmat_BoxVolume           	= intraODEinput.repmat_BoxVolume;     	% (m3); (2D matrix: num_grps X num_boxes)
num_grps                    = ODEinput.num_grps;                    % FFF put these into intraODEinput?
num_boxes                   = ODEinput.num_boxes;                   % FFF put these into intraODEinput?



% *************************************************************************
% STEP 2: prepare box biomass terms----------------------------------------
% step 2a: intra-domain source & export box biomasses ---------------------
compact_biomass_source      = biomass_plus1(:, looky_flux(:, 2));            % biomasses of each source box exporting water; (mmoles N/m3); (2D matrix: num_grps X num_fluxes); sinking units (mmoles N/m2/d); NOTE: matlab automatically truncates the singleton dimension in biomass_plus1
biomass_boundary_export     = biomass_plus1(:, looky_boundary_export(:, 2)); % biomasses of each source box exporting water out of domain; (mmoles N/m3); (2D matrix: num_grps X num_export_fluxes); sinking units (mmoles N/m2/d); NOTE: matlab automatically truncates the singleton dimension in biomass_plus1

% step 2b: boundary conditions for import boxes ---------------------------
biomass_boundary_import     = biomass_boundary_t(:, looky_boundary_import(:, 1)); % biomasses immediately outside the model domain of each destiny box that imports water from outside the model domain; (mmoles N/m3); (2D matrix: num_grps X num_import_fluxes); sinking units (mmoles N/m2/d); NOTE: matlab automatically truncates the singleton dimension in biomass_boundary_t



% *************************************************************************
% step 3: express fluxes in terms of biomass exchange----------------------
flux_biomass                = compact_flux_t .* compact_biomass_source;        % flux biomasses of each source box exporting water; (mmoles N/d); (2D matrix: num_grps X num_fluxes); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)
flux_domain_import          = compact_flux_t(:, looky_boundary_import(:, 3)) .* biomass_boundary_import; % flux biomasses into each destiny box that imports water from outside the model domain; (mmoles N/d); (2D matrix: num_grps X num_import_fluxes); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)
flux_domain_export          = compact_flux_t(:, looky_boundary_export(:, 3)) .* biomass_boundary_export; % flux biomasses out of each source box that exports water out of the model domain; (mmoles N/d); (2D matrix: num_grps X num_export_fluxes); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)



% *************************************************************************
% STEP 4: pool fluxes in source & destiny boxes----------------------------
flux_import_t               = accumarray([repmat_grp_row(:) repmat_looky_destiny(:)], flux_biomass(:));	% pool import to destiny boxes; (mmoles N/d); (3D matrix: num_grps X num_boxes+1); NOTE: does NOT include domain import (see flux_domain_import for domain input rates)
flux_export_t               = accumarray([repmat_grp_row(:) repmat_looky_source(:)],  flux_biomass(:));	% pool export from source boxes; (mmoles N/d); (3D matrix: num_grps X num_boxes+1); NOTE: does include domain export
flux_import_t(:, end)   	= []; % delete boundary term; (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_export_t(:, end)   	= []; % delete boundary term; (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)



% *************************************************************************
% STEP 5: prepare cross-domain boundary fluxes-----------------------------
flux_domain_import_t(:, looky_boundary_import(:, 1))	= flux_domain_import;     % (mmoles N/d); (2D matrix: num_grps X num_boxes+1)
flux_domain_export_t(:, looky_boundary_export(:, 2))	= flux_domain_export;     % (mmoles N/d); (2D matrix: num_grps X num_boxes+1)
flux_domain_import_t(:, end)    = []; % delete boundary term; (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_export_t(:, end)	= []; % delete boundary term; (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)


% step 5b: calculate driver group production (flux_domain_import_driver_t)
%          NOTE: keep production_input groups (e.g., NO3) separate from other
%                boundary condition groups as driver_input_t is treated
%                separately in dy calculation
flux_domain_import_driver_t(looky_driver, :)	= flux_domain_import_t(looky_driver, :);	% (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_import_t(looky_driver, :)           = 0;                                        % zero out driver group(s); (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)



% *************************************************************************
% STEP 6: apply retention scalers------------------------------------------
flux_import_t              	= flux_import_t               .* (1-RetentionScaler); % (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_export_t            	= flux_export_t               .* (1-RetentionScaler); % (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_import_t     	= flux_domain_import_t        .* (1-RetentionScaler); % (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_export_t      	= flux_domain_export_t        .* (1-RetentionScaler); % (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_import_driver_t	= flux_domain_import_driver_t .* (1-RetentionScaler); % (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)



% *************************************************************************
% STEP 7: biomass flux density relative to detination box volume-----------
flux_import_t             	= flux_import_t               ./ repmat_BoxVolume; % (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_export_t             	= flux_export_t               ./ repmat_BoxVolume; % (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_import_t        = flux_domain_import_t        ./ repmat_BoxVolume; % (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_export_t      	= flux_domain_export_t        ./ repmat_BoxVolume; % (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINY)
flux_domain_import_driver_t	= flux_domain_import_driver_t ./ repmat_BoxVolume; % (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINY)



% *************************************************************************
% STEP 8: calculate flux_net_t---------------------------------------------
flux_net_t                	= (flux_import_t + flux_domain_import_t) - flux_export_t; % net biomass import to each box; (mmoles N/d); (2D matrix: num_grps X num_boxes DESTINY)



% *************************************************************************
% STEP 9: reshape flux matrices--------------------------------------------
flux_import_t             	= reshape(flux_import_t,  [1, num_grps, num_boxes]);                % transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
flux_export_t             	= reshape(flux_export_t,  [1, num_grps, num_boxes]);                % transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
flux_net_t                  = reshape(flux_net_t,  [1, num_grps, num_boxes]);                   % transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
flux_domain_import_t        = reshape(flux_domain_import_t,  [1, num_grps, num_boxes]);         % transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
flux_domain_export_t        = reshape(flux_domain_export_t,  [1, num_grps, num_boxes]);         % transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
flux_domain_import_driver_t = reshape(flux_domain_import_driver_t,  [1, num_grps, num_boxes]);	% transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)



% *************************************************************************
% STEP 10: pack results----------------------------------------------------
Fluxes_t.flux_import_t             	 = flux_import_t;               % (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
Fluxes_t.flux_export_t             	 = flux_export_t;               % (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
Fluxes_t.flux_net_t                  = flux_net_t;                  % (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
Fluxes_t.flux_domain_import_t        = flux_domain_import_t;        % (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
Fluxes_t.flux_domain_export_t        = flux_domain_export_t;        % (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
Fluxes_t.flux_domain_import_driver_t = flux_domain_import_driver_t;	% (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)


% end m-file **************************************************************