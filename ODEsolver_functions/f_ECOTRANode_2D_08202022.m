function dy = f_ECOTRANode_2D_08202022(t, ProductionRates_t, ODEinput)
% generic ODE for ECOTRAN_DynamicScenario
% default is for reflective boundary, but has built-in options for defined boundary conditions
% by Jim Ruzicka
%
% calls:
%       f_PhysicalFlux_intraODE_09092019    calculate biomass fluxes into and out of each box and across domain boundaries (mmoles N/m3/d)	(3D matrix: 1 X num_grps X num_boxes DESTINATION)
%       f_MichaelisMenten_05152016          (optional) Michaelis-Menton uptake kinetics for phytoplankton
%
% takes:
%   ODEinput
%       (trophic model info)
%           EnergyBudget              consumption fates                                  (fraction)      (3D matrix: num_grps X num_grps X num_boxes)
%           EnergyBudget_BoxType             turn off specific trophic linkages in sub-surface boxes            (3D matrix: num_grps X num_grps X num_boxes)
%           ConsumptionBudget_feces         consumption fraction to feces                                      (3D matrix: num_t X num_grps X num_boxes)
%           ConsumptionBudget_metabolism    consumption fraction to metabolism                                 (3D matrix: num_t X num_grps X num_boxes)
%           ConsumptionBudget_predation     consumption fraction to predation                                  (3D matrix: num_t X num_grps X num_boxes)
%           ConsumptionBudget_eggs          consumption fraction to eggs                                       (3D matrix: num_t X num_grps X num_boxes)
%           ConsumptionBudget_senescence    consumption fraction to Senescence                                 (3D matrix: num_t X num_grps X num_boxes)
%           ConsumptionBudget_ba            consumption fraction to ba                                         (3D matrix: num_t X num_grps X num_boxes)
%           ConsumptionBudget_em            consumption fraction to em                                         (3D matrix: num_t X num_grps X num_boxes)
%           fate_feces                                                                                  (3D matrix: 2 X num_grps X num_boxes)
%           fate_senescence                                                                             (3D matrix: 2 X num_grps X num_boxes)
%           fate_metabolism                                                                                       (3D matrix: 2 X num_grps X num_boxes)
%           fate_eggs                                                                                             (3D matrix: num_eggs X num_grps X num_boxes)
%           pb                              intrinsic (weight-specific) growth rates (not used) (1/d)           (3D matrix: num_t X num_grps X num_boxes)
%           qb         intrinsic (weight-specific) consumption rates       (1/d)           (3D matrix: num_t X num_grps X num_boxes)
%           TransferEfficiency                                                                  (fraction)      (3D matrix: 1 X num_grps X num_boxes)
%       (production rates: initial conditions, external inputs)
%           production_input                daily advection rates of NO3 INTO model domain      (mmole N/m3/d)  (3D matrix: num_t X num_grps X num_boxes)
%           production_initial              initial & mean production rates (consumption)       (mmole N/m3/d)  (3D matrix of vertical vectors: num_grps X 1 X num_boxes)
%           productionC_initial_repmat      mean production rates for functional response calcs (mmole N/m3/d)  (3D matrix: num_grps X num_grps X num_boxes); NOTE: vertical vectors
%           BoxBiomass_OceanSurface         ocean NO3 conc. for horizontal mixing calcs         (mmole N/m3)    (vertical vector: num_grps X 1); NOTE: could define ALL group ocean biomasses here
%           BoxBiomass_OceanSubSurface      ocean NO3 conc. for horizontal mixing calcs         (mmole N/m3)    (vertical vector: num_grps X 1); NOTE: could define ALL group ocean biomasses here
%       (addresses and numbers in EnergyBudget)
%           looky_NO3
%           looky_plgcNH4
%           looky_bnthNH4
%           looky_nonNO3
%           looky_ANYdetritus
%           looky_terminalPLGCdetritus
%           looky_terminalBNTHdetritus
%           looky_ANYconsumer
%           looky_fleets
%           looky_eggs
%           num_eggs
%           num_bnthNH4
%           num_grps                        number of ECOTRAN producers or consumers
%       (physical model info)
%       	AdvectionMatrix                 net advection rate                                  (m3/d)          (3D matrix: time X SourceBox X DestinationBox (I II III IV V))
%       	RiverFluxMatrix                 net advection rate                                  (m3/d)          (3D matrix: time X SourceBox X DestinationBox (I II III IV V))
%       	VerticalMixMatrix               vertical mixing                                     (m3/d)          (3D matrix: groups X SourceBox X DestinationBox)
%       	HorizontalMixMatrix             horizontal mixing                                   (m3/d)          (3D matrix: groups X SourceBox X DestinationBox)
%           SinkFrac                        fraction of total production sinking OUT of source box              (3D matrix: num_t X num_grps X num_boxes)
%           t_grid                          day number                                                          (vertical vector: time x 1)
%       (addresses and numbers in physical model)
%           looky_SurfaceBoxes
%           looky_SubSurfaceBoxes
%           looky_OffshoreSurfaceBox
%           looky_OffshoreSubSurfaceBox
%           num_boxes                       number of spatial boxes in model
%       (functional response parameters)
%           FunctionalResponseParams        producer vulnerability (m_p) (3D matrix: group X PRODUCERS X num_boxes), replicated down rows
%                   m_p = 0 is "linear"
%                   m_p = 1 is "non-linear"  ECOSIM default; half-sat constant of individual groups 
%       (used only for Michaelis-Menten purposes)
%           looky_ANYPrimaryProducer
%           looky_macroalgae
%           num_PrimaryProducers
%           num_macroalgae
%           Kw                              light attenuation_seawater                                          (scaler)
%           Kp                              ight attenuation_phytoplankton                      (m2/mmol N)     (scaler)
%           MichaelisMenten_Vmax            Vmax  = maximum nutrient uptake rate                (1/d)           (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%           MichaelisMenten_KNO3            KNO3  = NO3 half-saturation constant                (mmol N/m3)     (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%           MichaelisMenten_KNH4            KNH4  = NH4 half-saturation constant                (mmol N/m3)     (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%           MichaelisMenten_alpha           alpha = initial slope of light response curve       (m2/W/d)        (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%           MichaelisMenten_psi             psi   = NO3 uptake inhibition by NH4                (m3/mmole N)    (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%   time-step & current production rate supplied by ode45
%           t                               current time-step
%       ProductionRates_t                   (AKA "y") current production rates of each group    (t WWT/km2*m/d) (vertical vectors: groups X 1 X num_boxes)
%
% returns:
%           dy                              daily rate of change of production rate              (t/km2/d2)     (vertical vector)
%
% NOTE: ProductionEfficiency .* consumption_IN converts consumption_IN from
%       consumption to predation and accounts for metabolism & non-assimilation as linear terms
% NOTE: TransferEfficiency .* consumption_IN converts consumption_IN from
%       consumption to predation and accounts for metabolism &
%       non-assimilation as linear terms. 
% NOTE: In ECOTRAN, Senescence, Non-Assimilation, & metabolism are
%       accounted for in the EnergyBudget, and are therefore
%       mathematically equivalent to functional num_grps, and flow to
%       Senescence, Non-Assimilation, & metabolism is handled the same way as
%       flow to predators.
% NOTE: Reflective ocean boundary conditions. (Just a cleaned-up version of f_ECOTRANode_09302016.m)
%
% revision date: 8-20-2022
%       8/13/2022   transfering benthic NH4 to sub-surface boxes (step 9c)
%       8/13/2022   corrections for proper use of pb & qb (qb for biomass calc; pb for dy calc)
%       8/13/2022   set FunctionalResponse value to 0 when group production @ t is negative (should not matter in matlab, but will be important in C++?)
%       8/16/2022   tested old method of moving bnthNH4 to sub-surface boxes, but results were faulty. Trying new technique with changes to ProductionRates_t, fate_metabolism, and/or EnergyBudget
%       8/20/2022   minor clean-up


% *************************************************************************
% STEP 1: unpack parameters & drivers from ODEinput------------------------

% step 1a: EnergyBudget ---------------------------------------------------
EnergyBudget                        = ODEinput.EnergyBudget;                    % (3D matrix: num_grps X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 1b: ConsumptionBudget terms ----------------------------------------
ConsumptionBudget_feces             = ODEinput.ConsumptionBudget_feces;        	% (3D matrix: num_t X num_grps X num_boxes)
ConsumptionBudget_metabolism        = ODEinput.ConsumptionBudget_metabolism;	% (3D matrix: num_t X num_grps X num_boxes)
ConsumptionBudget_eggs              = ODEinput.ConsumptionBudget_eggs;         	% (3D matrix: num_t X num_grps X num_boxes)
ConsumptionBudget_predation         = ODEinput.ConsumptionBudget_predation;    	% (3D matrix: num_t X num_grps X num_boxes)
ConsumptionBudget_senescence        = ODEinput.ConsumptionBudget_senescence;	% (3D matrix: num_t X num_grps X num_boxes)
ConsumptionBudget_ba                = ODEinput.ConsumptionBudget_ba;          	% (3D matrix: num_t X num_grps X num_boxes)
ConsumptionBudget_em                = ODEinput.ConsumptionBudget_em;          	% (3D matrix: num_t X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 1c: TransferEfficiency ---------------------------------------------
TransferEfficiency                  = ODEinput.TransferEfficiency;              % (3D matrix: 1 X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 1d: group addresses & counts ---------------------------------------
looky_driver                        = ODEinput.looky_driver;	% row address(es) of driver group(s) (e.g., NO3)
looky_nutrients                     = ODEinput.looky_nutrients;
% looky_NO3                           = ODEinput.looky_NO3;       % NOTE: use with Michaelis-Menten option
looky_plgcNH4                       = ODEinput.looky_plgcNH4;
looky_bnthNH4                       = ODEinput.looky_bnthNH4;
% looky_ANYPrimaryProducer            = ODEinput.looky_ANYPrimaryProducer; % NOTE: use with Michaelis-Menten option
looky_eggs                          = ODEinput.looky_eggs;
% looky_fleets                        = ODEinput.looky_fleets;
looky_livingANDfleets               = ODEinput.looky_livingANDfleets;
looky_ANYdetritus                   = ODEinput.looky_ANYdetritus;

num_grps                            = ODEinput.num_grps;
num_nutrients                       = ODEinput.num_nutrients;
num_eggs                            = ODEinput.num_eggs;
num_ANYdetritus                     = ODEinput.num_ANYdetritus;
num_livingANDfleets                 = ODEinput.num_livingANDfleets;
num_boxes                           = ODEinput.num_boxes;       % number of spatial boxes in model
% -------------------------------------------------------------------------


% step 1e: physics terms --------------------------------------------------
t_grid                              = ODEinput.t_grid;                      % (vertical vector: num_t x 1)
BoxVolume                           = ODEinput.BoxVolume;                   % (m3); (2D matrix: time X num_boxes)
ADVECTION_compact                   = ODEinput.ADVECTION_compact;           % (m3/d); (2D matrix: num_t X num_fluxes_advection)
HORIZONTALMIXING_compact            = ODEinput.HORIZONTALMIXING_compact;	% (m3/d); (2D matrix: num_t X num_fluxes_HorizontalMixing)
VERTICALMIXING_compact              = ODEinput.VERTICALMIXING_compact;      % (m3/d); (2D matrix: num_t X num_fluxes_VerticalMixing)

SinkLink_surface                    = ODEinput.SinkLink_surface; % (3D matrix: num_boxes DESTINY X 1 X num_boxes SOURCE)
% SinkLink_benthos                    = ODEinput.SinkLink_benthos; % (3D matrix: destination boxes X 1 X source boxes)
% -------------------------------------------------------------------------


% step 1f: initial & boundary conditions, external driver -----------------
productionC_initial_repmat          = ODEinput.productionC_initial_repmat;	% consumption inflow rate; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: replicated vertical vectors across columns
% biomass_boundary                    = ODEinput.biomass_boundary;            % NOTE: use this line only for NON-REFLECTIVE boundary option; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes+1)
external_driver                     = ODEinput.external_driver;             % external boundary biomass conditions for each box; external NO3 driver; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
%   FFF external_driver defined for all boxes now, but will try to compact to just the boxes with fluxes (3D matrix: num_t X num_drivers X num_fluxes_BoundaryImport)
%   FFF could allow for allow multiple external_driver grps?
% -------------------------------------------------------------------------


% step 1g: functional group parameters ------------------------------------
pb                                  = ODEinput.pb;                % intrinsic (weight-specific) production rates; (1/d); (3D matrix: num_t X num_grps X num_boxes)
qb                                  = ODEinput.qb;                % intrinsic (weight-specific) consumption rates; (1/d); (3D matrix: num_t X num_grps X num_boxes)
SINKING_compact                     = ODEinput.SINKING_compact;   % (m3/d); (2D matrix: num_t X num_grps X num_fluxes_sinking)
MIGRATION_compact                   = ODEinput.MIGRATION_compact; % (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
% biomass_migrator                    = ODEinput.biomass_migrator;  % special definition of boundary biomasses for migrators; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes); NOTE: deactived by default
% -------------------------------------------------------------------------


% step 1h: functional group fate parameters -------------------------------
fate_feces                          = ODEinput.fate_feces;   	% (3D matrix: num_ANYdetritus X num_grps X num_boxes)
fate_metabolism                     = ODEinput.fate_metabolism;	% (3D matrix: num_nutrients X num_grps X num_boxes)
fate_eggs                           = ODEinput.fate_eggs;    	% (3D matrix: num_eggs X num_grps X num_boxes)
fate_predation                      = ODEinput.fate_predation;  % (3D matrix: num_livingANDfleets X num_grps X num_boxes)
fate_senescence                     = ODEinput.fate_senescence;	% (3D matrix: num_ANYdetritus X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 1i: functional response parameters -------------------------------
FunctionalResponseParams            = ODEinput.FunctionalResponseParams; % producer vulnerability (m_p); (3D matrix: CONSUMERS X prey group X num_boxes) replicated across clms (= producers)
% -------------------------------------------------------------------------


% % step 1j: scenario rules -------------------------------------------------
% %          NOTE: this could be anything, using a vertical vector for scaling NO3 levels for now
% % NNN come back to this - find dynmic scenario code in other file
% Scenario_scaler                     = ODEinput.Scenario_scaler;     % scenario time-series (scaler); (3D matrix: time X num_grps X num_boxes)
% -------------------------------------------------------------------------


% % step 1k: terms used for Michaelis-Menten option -------------------------
% MichaelisMenten_Vmax                = ODEinput.MichaelisMenten_Vmax;      % Vmax  = maximum nutrient uptake rate; (1/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% MichaelisMenten_KNO3                = ODEinput.MichaelisMenten_KNO3;      % KNO3  = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% MichaelisMenten_KNH4                = ODEinput.MichaelisMenten_KNH4;      % KNH4  = NH4 half-saturation constant; (mmol N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% MichaelisMenten_alpha               = ODEinput.MichaelisMenten_alpha;     % alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% MichaelisMenten_psi                 = ODEinput.MichaelisMenten_psi;       % psi   = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_PrimaryProducers X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% MichaelisMenten_w                   = ODEinput.MichaelisMenten_w;         % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% MichaelisMenten_eta                 = ODEinput.MichaelisMenten_eta;       % eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% 
% MLD                                 = ODEinput.MLD;                       % time-series of Mixed Layer Depth (for light intensity calculations, NOT advection calcs); (m); (vertical vector: length = num_t); NOTE: use for Michaelis-Menten option
% Io                                  = ODEinput.Io;                        % time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); NOTE: the daily average is also what was used by Spitz in her NPZD models; NOTE: use for Michaelis-Menten function
% Kw                                  = ODEinput.Kw;                        % Light attenuation of seawater; (scaler); NOTE: use for Michaelis-Menten function
% Kp                                  = ODEinput.Kp;                        % Light attenuation of phytoplankton (Newberger et al., 2003); (m2/mmol N); NOTE: use for Michaelis-Menten function
% 
% ProductionFraction_macroalgae       = ODEinput.ProductionFraction_macroalgae; % tie any macroalgae production to a fixed multiple of total phytoplankton production in each box; (3D matrix: num_macroalgae X 1 X num_boxes)
% *************************************************************************





% *************************************************************************
% STEP 2: prepare ProductionRates_t & initialize results-------------------

% step 2a: "unstack" vertical vector of ProductionRates_t -----------------
ProductionRates_t           = reshape(ProductionRates_t, [num_grps, 1, num_boxes]); % (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)

% % QQQ for de-bugging
% t = 132; ProductionRates_t = ODEinput.production_initial; % (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
% -------------------------------------------------------------------------


% step 2b: filter out sub-minimum (negative) Production @ t ---------------
%          NOTE: this should only be necessary to capture tiny errors from ODE solver
ProductionRates_t(ProductionRates_t < 0)	= 0; % (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
% *************************************************************************





% *************************************************************************
% STEP 3: interpolate parameters @ t---------------------------------------

% step 3a: ConsumptionBudget terms @ t ------------------------------------
%          NOTE: this allows for defined seasonal changes in ConsumptionBudget terms
ConsumptionBudget_feces_t           = interp1(t_grid, ConsumptionBudget_feces,       t); % ConsumptionBudget_feces @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
ConsumptionBudget_metabolism_t      = interp1(t_grid, ConsumptionBudget_metabolism,  t); % ConsumptionBudget_metabolism @ t; (fraction); (3D matrix: 1 X num_grps X num_boxes)
ConsumptionBudget_eggs_t            = interp1(t_grid, ConsumptionBudget_eggs,        t); % ConsumptionBudget_eggs @ t;       (fraction); (3D matrix: 1 X num_grps X num_boxes)
ConsumptionBudget_predation_t       = interp1(t_grid, ConsumptionBudget_predation,   t); % ConsumptionBudget_predation @ t;  (fraction); (3D matrix: 1 X num_grps X num_boxes)
ConsumptionBudget_senescence_t      = interp1(t_grid, ConsumptionBudget_senescence,  t); % ConsumptionBudget_senescence @ t; (fraction); (3D matrix: 1 X num_grps X num_boxes)
ConsumptionBudget_ba_t              = interp1(t_grid, ConsumptionBudget_ba,          t); % ConsumptionBudget_ba @ t;         (fraction); (3D matrix: 1 X num_grps X num_boxes)
ConsumptionBudget_em_t              = interp1(t_grid, ConsumptionBudget_em,          t); % ConsumptionBudget_em @ t;         (fraction); (3D matrix: 1 X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 3b: physical flux rates @ t ----------------------------------------
ADVECTION_compact_t                 = interp1(t_grid, ADVECTION_compact, t);        % (m3/d); (horizontal vector: 1 X num_fluxes_advection); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
HORIZONTALMIXING_compact_t          = interp1(t_grid, HORIZONTALMIXING_compact, t);	% (m3/d); (horizontal vector: 1 X num_fluxes_HorizontalMixing); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
VERTICALMIXING_compact_t            = interp1(t_grid, VERTICALMIXING_compact, t);	% (m3/d); (horizontal vector: 1 X num_fluxes_VerticalMixing); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
% -------------------------------------------------------------------------


% step 3c: box dimensions @ t ---------------------------------------------
BoxVolume_t                         = interp1(t_grid, BoxVolume, t);                % (m3); (horizontal vector: 1 X num_boxes)
% -------------------------------------------------------------------------


% % step 3d: boundary biomass @ t -----------------------------------------
% %          NOTE: use only with NON-REFLECTIVE boundary conditions
% %          FFF: simplify to 2D matrix and define only for boxes with defined boundary physical fluxes
% biomass_boundary_t                  = interp1(t_grid, biomass_boundary, t);                      % (mmoles N/m3); (3D matrix: 1 X num_grps X num_boxes+1)
% biomass_boundary_t                  = reshape(biomass_boundary_t, [num_grps, 1, (num_boxes+1)]); % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
% -------------------------------------------------------------------------


% step 3e: external driver biomass @ t ------------------------------------
%          FFF: defining for all boxes now, but will try to trim down to just the boxes with fluxes (horizontal vector: 1 X num_fluxes_BoundaryImport)
external_driver_t                   = interp1(t_grid, external_driver, t);	% external driver (e.g., NO3) boundary biomass conditions for each box @ t; (mmole N/m3); (3D matrix: 1 X num_drivers X num_boxes+1)
% -------------------------------------------------------------------------


% step 3f: weight-specific production & consumption rates @ t -------------
pb_t                                = interp1(t_grid, pb, t);               % weight-specific production rates(1/d); (3D matrix: 1 X num_grps X num_boxes)
qb_t                                = interp1(t_grid, qb, t);               % weight-specific consumption rates(1/d); (3D matrix: 1 X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 3g: sinking & migration @ t ----------------------------------------
SINKING_compact_t                   = interp1(t_grid, SINKING_compact, t);	 % sinking flux OUT of SOURCE box; (m3/d); (3D matrix: 1 X num_grps X num_fluxes_sinking); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
SINKING_compact_t                   = squeeze(SINKING_compact_t);            % sinking flux OUT of SOURCE box; (m3/d); (2D matrix: num_grps X num_fluxes_sinking)
MIGRATION_compact_t                 = interp1(t_grid, MIGRATION_compact, t); % migration flux OUT of SOURCE box; (m3/d); (3D matrix: 1 X num_grps X num_fluxes_migration)
MIGRATION_compact_t                 = squeeze(MIGRATION_compact_t);          % migration flux OUT of SOURCE box; (m3/d); (2D matrix: num_grps X num_fluxes_migration)
% biomass_migrator_t                  = interp1(t_grid, biomass_migrator, t);  % NOTE: use only for special definition of boundary biomass for migrators; boundary biomass for migrators @ t; (mmoles N/m3); (3D matrix: 1 X num_grps X num_boxes); NOTE: deactived by default
% biomass_migrator_t                  = reshape(biomass_migrator_t, [num_grps, 1, num_boxes]); % NOTE: use only for special definition of boundary biomass for migrators; transpose biomass_migrator_t to vertical vectors; (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes); NOTE: deactived by default
% -------------------------------------------------------------------------


% % step 3h: scenario rules @ t ---------------------------------------------
% Scenario_scaler_t         = interp1(t_grid, Scenario_scaler, t) - 1; % QQQ remove -1 ; production_input @ t; external input to model domain; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
% % -------------------------------------------------------------------------


% % step 3i: light terms @ t ----------------------------------------------
% %          NOTE: use for Michaelis-Menten option
% MLD_t                     = interp1(t_grid, MLD, t); % Mixed Layer Depth @ t; (m); (scaler)
% Io_t                      = interp1(t_grid, Io, t); % surface PAR light intensity @ t; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (scaler); 
% *************************************************************************





% *************************************************************************
% STEP 4: Adjust EnergyBudget to accomodate changes in ConsumptionBudget---
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
% QQQ any need for ULTIMATEdetritus correction?
% QQQ ba is implied through effects on other budget terms (allowing sum of
%     EnergyBudget to be >< 1?? TEST THIS TO MAKE SURE I DONT NEED TO
%     EXPLICITLY DEFINE THIS IN dy CALC
% QQQ em should be included in dy calc and handled as with physics

% step 4a: ConsumptionBudget @ t: metabolism into EnergyBudget ------------
EnergyBudget(looky_nutrients, :, :)         = repmat(ConsumptionBudget_metabolism_t, [num_nutrients, 1, 1]) .* fate_metabolism; % excretion of NH4 and nitrification of NH4-->>NO3; (fraction); (3D matrix: num_nutrients X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 4b: ConsumptionBudget @ t: egg production into EnergyBudget --------
EnergyBudget(looky_eggs, :, :)              = repmat(ConsumptionBudget_eggs_t, [num_eggs, 1, 1]) .* fate_eggs; % egg production (should work for [] eggs or for multiple eggs); (fraction); (3D matrix: num_eggs X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 4c: ConsumptionBudget @ t: predation & fleet catch into EnergyBudget
EnergyBudget(looky_livingANDfleets, :, :)	= repmat(ConsumptionBudget_predation_t, [num_livingANDfleets, 1, 1]) .* fate_predation; % (fraction); (3D matrix: num_livingANDfleets X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 4d: ConsumptionBudget @ t feces & senescence into EnergyBudget -----
feces_t                                     = repmat(ConsumptionBudget_feces_t, [num_ANYdetritus, 1, 1]) .* fate_feces; % (3D matrix: num_ANYdetritus X num_grps X num_boxes)
senescence_t                                = repmat(ConsumptionBudget_senescence_t, [num_ANYdetritus, 1, 1]) .* fate_senescence; % (3D matrix: num_ANYdetritus X num_grps X num_boxes)
detritus_t                                  = feces_t + senescence_t;	% detritus = feces + senescence; (fraction); (3D matrix: num_ANYdetritus X num_grps X num_boxes)
EnergyBudget(looky_ANYdetritus, :, :)       = detritus_t;               % detritus = feces + senescence; (fraction); (3D matrix: num_ANYdetritus X num_grps X num_boxes)

% % % QQQ scaling of predation
% EnergyBudget(looky_livingANDfleets, 23, :)	= 10 * ((repmat(ConsumptionBudget_predation_t(1, 23, :), [24, 1, 1])) .* fate_predation(:, 23, :)); % (fraction); (3D matrix: num_livingANDfleets X num_grps X num_boxes)
% *************************************************************************





% *************************************************************************
% STEP 5: prepare Production repmat matrices for array multiplication------

% step 5a: apply Scenario_scaler_t to ProductionRates_t -------------------
%          NNN find dynamic scenario work done in other file
% % % Scenario_scaler_t           = reshape(Scenario_scaler_t, [num_grps 1 num_boxes]); % (3D matrix: num_grps X 1 X num_boxes) QQQ check that reshaping is correct
% % % % ProductionRates_t           = ProductionRates_t .* Scenario_scaler_t; % deactivated ... ; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
% % % 
% % % % ConsumptionBudget_em_t = ConsumptionBudget_em_t + (-1*Scenario_scaler_t); % QQQ an attempt to get NO3 biomass doubled in box 5
% -------------------------------------------------------------------------


% step 5b: reshape and repmat ProductionRates_t ---------------------------
ProductionRatesC_t_repmat   = repmat(ProductionRates_t, [1, num_grps, 1]);              % (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes) NOTE: VERTICAL vectors replicated across clms;
ProductionRates_t_transpose = reshape(ProductionRates_t, [1, num_grps, num_boxes]);     % reshape to 3D matrix ; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes); 
ProductionRatesP_t_repmat   = repmat(ProductionRates_t_transpose, [num_grps, 1, 1]);	% transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: HORIZONTAL vectors replicated down rows
% *************************************************************************





% *************************************************************************
% STEP 6: biomasses, nutrient concentrations, sinking flux, & boundary conditions @ t
%         FFF Keep a running biomass time-series

% step 6a: biomass @ t = q / (q/b) ----------------------------------------
biomass_t                       = ProductionRates_t ./ reshape(qb_t, [num_grps, 1, num_boxes]); % use qb_t value to convert rates to standing stock biomass; (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes); NOTE transpose of qb_t to vertical vectors
% biomass_PrimaryProducer_t       = biomass_t(looky_ANYPrimaryProducer, 1, :); % primary producer biomasses; (mmole N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes); NOTE: use with Michaelis-Menten option

biomass_plus1                	= biomass_t;        % (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes)
biomass_plus1(:, 1, (end+1))	= 0;                % add 1 clm for boundary fluxes; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)

% biomass_migrator_plus1        	= biomass_migrator_t; % special definition of boundary biomass for migrators; (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes); NOTE: deactived by default
% biomass_migrator_plus1(:, 1, (end+1))	= 0;          % add 1 clm for boundary fluxes; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: deactived by default
% -------------------------------------------------------------------------


% step 6b: boundary conditions & external driver conditions ---------------
%	OPTION 1: reflective boundary
%             NOTE: biomass is imported into each boundary box from external environment of same biomass densities
%             NOTE: boundary conditions are defined for each domain box
%                   whether that box is on the edge of the domain or not 
%                   (if not, the boundary values are not used)
    biomass_boundary_t                   	= biomass_plus1;        % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
    biomass_boundary_t(looky_driver, 1, :)	= external_driver_t;	% paste in external (boundary) driver conditions (NO3) @ t; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: matlab automatically makes transpose of rows & clms during assignment
    
    biomass_MigratorBoundary_t              = biomass_plus1;        % default case is to use same boundary conditions for migration as for physical fluxes; NOTE: deactivate this line for special definition of boundary biomass for migrators (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
%     biomass_MigratorBoundary_t              = biomass_migrator_plus1; % NOTE: reactivate this line when using special definition of boundary biomass for migrators; (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: deactived by default

    
% %   OPTION 2: defined boundary conditions
% %             NOTE: use for defined boundary conditions
% %             FFF: try to trim to 2D matrix
% 	biomass_boundary_t(looky_driver, 1, :)	= external_driver_t;	  % paste in external (boundary) driver conditions (NO3) @ t; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: matlab automatically makes transpose of rows & clms during assignment
% *************************************************************************





% *************************************************************************
% STEP 7: calculate advection, mixing, & sinking exchange rates between boxes and out of domain
%         NOTE: flux_domain_import_t value = 0 for driver group, driver group boundary input is given in external_driver_t
%         NOTE: NetFlux_t includes the flux_domain_export_t
%         NOTE: flux_domain_export_t values are all positive even though this is a loss from the box
%         NOTE: code currently only works with 1 defined external driver group (e.g. NO3) (FFF could chnage this in the future)

% step 7a: define initialized variables -----------------------------------
intraODEinput.flux_domain_import_t          = ODEinput.flux_domain_import_t;       	% initialized variable, all zeros; (2D matrix: num_grps X num_boxes+1)
intraODEinput.flux_domain_export_t          = ODEinput.flux_domain_export_t;      	% initialized variable, all zeros; (2D matrix: num_grps X num_boxes+1)
intraODEinput.flux_domain_import_driver_t	= ODEinput.flux_domain_import_driver_t;	% initialized variable, all zeros; (2D matrix: num_grps X num_boxes DESTINY)
intraODEinput.repmat_BoxVolume              = repmat(BoxVolume_t, [num_grps, 1]);  	% (m3); (2D matrix: num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 7b: ADVECTION ------------------------------------------------------
intraODEinput.compact_flux_t            = repmat(ADVECTION_compact_t, [num_grps, 1]); % replicate each flux term for all model groups; (m3/d); (2D matrix: num_grps X num_fluxes_advection)
intraODEinput.biomass_plus1             = biomass_plus1;                            % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.biomass_boundary          = biomass_boundary_t;                      	% boundary biomass conditions @ t; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler;                 % (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINY);
intraODEinput.looky_flux                = ODEinput.looky_AdvectionFlux;             % (2D matrix: num_fluxes_advection X 5-->[(DESTINY box) (SOURCE box) (DESTINY box address in compact_flux) (SOURCE box address in compact_flux) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
intraODEinput.looky_boundary_import     = ODEinput.looky_AdvectionBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
intraODEinput.looky_boundary_export     = ODEinput.looky_AdvectionBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
intraODEinput.looky_driver           	= looky_driver;                             % row address(es) of driver group(s) (e.g., NO3)
intraODEinput.repmat_looky_source       = ODEinput.repmat_looky_ADVECTION_source;	% (2D matrix: num_grps X num_fluxes_advection); NOTE: clm 2 in looky_flux
intraODEinput.repmat_looky_destiny      = ODEinput.repmat_looky_ADVECTION_destiny;	% (2D matrix: num_grps X num_fluxes_advection); NOTE: clm 1 in looky_flux
intraODEinput.repmat_grp_row          	= ODEinput.repmat_GrpRow_ADVECTION;         % addressses; (2D matrix: num_grps X num_fluxes)
Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
flux_net_advection_t                    = Fluxes_t.flux_net_t;                      % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
flux_domain_import_driver_advection_t	= Fluxes_t.flux_domain_import_driver_t;     % net biomass flux OF DRIVER GRP across domain boundary INTO each box (all other grps = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION)
% -------------------------------------------------------------------------


% step 7c: HORIZONTAL MIXING ----------------------------------------------
intraODEinput.compact_flux_t            = repmat(HORIZONTALMIXING_compact_t, [num_grps, 1]); % replicate each flux term for all model groups; (m3/d); (2D matrix: num_grps X num_fluxes_HorizontalMixing)
intraODEinput.biomass_plus1             = biomass_plus1;                                    % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.biomass_boundary          = biomass_boundary_t;                              	% boundary biomass conditions @ t; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler;                         % (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
intraODEinput.looky_flux               	= ODEinput.looky_HorizontalMixingFlux;              % (2D matrix: num_fluxes_HorizontalMixing X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
intraODEinput.looky_boundary_import    	= ODEinput.looky_HorizontalMixingBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
intraODEinput.looky_boundary_export    	= ODEinput.looky_HorizontalMixingBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
intraODEinput.looky_driver           	= looky_driver;                                     % row address(es) of driver group(s) (e.g., NO3)
intraODEinput.repmat_looky_source     	= ODEinput.repmat_looky_HORIZONTALMIXING_source;	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE: clm 2 in looky_flux
intraODEinput.repmat_looky_destiny     	= ODEinput.repmat_looky_HORIZONTALMIXING_destiny;	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE: clm 1 in looky_flux
intraODEinput.repmat_grp_row          	= ODEinput.repmat_GrpRow_HORIZONTALMIXING;          % addressses; (2D matrix: num_grps X num_fluxes_HorizontalMixing)
Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
flux_net_HorizMixing_t                  = Fluxes_t.flux_net_t;                              % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
flux_domain_import_driver_HorizMixing_t	= Fluxes_t.flux_domain_import_driver_t;             % net biomass flux OF DRIVER GRP across domain boundary INTO each box (all other grps = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION)
% -------------------------------------------------------------------------


% step 7d: VERTICAL MIXING ------------------------------------------------
intraODEinput.compact_flux_t            = repmat(VERTICALMIXING_compact_t, [num_grps, 1]); % replicate each flux term for all model groups; (m3/d); (2D matrix: num_grps X num_fluxes_VerticalMixing)
intraODEinput.biomass_plus1             = biomass_plus1;                                % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.biomass_boundary          = biomass_boundary_t;                          	% boundary biomass conditions @ t; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler;                     % (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
intraODEinput.looky_flux              	= ODEinput.looky_VerticalMixingFlux;           	% (2D matrix: num_fluxes_VerticalMixing X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
intraODEinput.looky_boundary_import   	= ODEinput.looky_VerticalMixingBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
intraODEinput.looky_boundary_export  	= ODEinput.looky_VerticalMixingBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
intraODEinput.looky_driver           	= looky_driver;                                 % row address(es) of driver group(s) (e.g., NO3)
intraODEinput.repmat_looky_source     	= ODEinput.repmat_looky_VERTICALMIXING_source;	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE: clm 2 in looky_flux
intraODEinput.repmat_looky_destiny    	= ODEinput.repmat_looky_VERTICALMIXING_destiny;	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE: clm 1 in looky_flux
intraODEinput.repmat_grp_row          	= ODEinput.repmat_GrpRow_VERTICALMIXING;        % addressses; (2D matrix: num_grps X num_fluxes_VerticalMixing)
Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
flux_net_VertMixing_t                   = Fluxes_t.flux_net_t;                          % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
% -------------------------------------------------------------------------


% step 7e: SINKING --------------------------------------------------------
%       QQQ NOTE: up until 2019, senescence mortality in EnergyBudget was
%       adjusted to accomodate sinking losses. I don't think this is
%       necessary. However, sinking of living particles to the benthos
%       box does require a complete conversion of group senescence to detritus.
%       Still need to make this adjustment in individual box budgets (outside of the ODE).
intraODEinput.compact_flux_t            = SINKING_compact_t;                        % (m3/d); (2D matrix: num_grps X num_fluxes_sinking)
intraODEinput.biomass_plus1             = biomass_plus1;                            % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.biomass_boundary          = biomass_boundary_t * 0;                 	% NOTE: no boundary conditions for sinking (all grps = 0); boundary biomass conditions @ t; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler * 0;             % NOTE: apply no retention for sinking; (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
intraODEinput.looky_flux              	= ODEinput.looky_SinkingFlux;               % (2D matrix: num_fluxes_sinking X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
intraODEinput.looky_boundary_import    	= ODEinput.looky_SinkingBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
intraODEinput.looky_boundary_export    	= ODEinput.looky_SinkingBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
intraODEinput.looky_driver           	= looky_driver;                             % row address(es) of driver group(s) (e.g., NO3)
intraODEinput.repmat_looky_source     	= ODEinput.repmat_looky_SINKING_source;     % (2D matrix: num_grps X num_fluxes_sinking); NOTE: clm 2 in looky_flux
intraODEinput.repmat_looky_destiny    	= ODEinput.repmat_looky_SINKING_destiny;	% (2D matrix: num_grps X num_fluxes_sinking); NOTE: clm 1 in looky_flux
intraODEinput.repmat_grp_row          	= ODEinput.repmat_GrpRow_SINKING;           % addressses; (2D matrix: num_grps X num_fluxes_sinking)
Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
flux_net_sinking_t                      = Fluxes_t.flux_net_t;                      % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
% -------------------------------------------------------------------------


% step 7f: MIGRATION ------------------------------------------------------
intraODEinput.compact_flux_t            = MIGRATION_compact_t;                      % (m3/d); (2D matrix: num_grps X num_fluxes_migration)
intraODEinput.biomass_plus1             = biomass_plus1;                            % (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.biomass_boundary          = biomass_MigratorBoundary_t;         	    % boundary biomass conditions @ t for migrators; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1)
intraODEinput.RetentionScaler         	= ODEinput.RetentionScaler * 0;             % (scaler = 0 for migration); (2D matrix: num_grps X num_boxes DESTINATION);
intraODEinput.looky_flux                = ODEinput.looky_MigrationFlux;             % (2D matrix: num_fluxes_migration X [(destiny box) (source box) (destiny box address) (source box address) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
intraODEinput.looky_boundary_import     = ODEinput.looky_MigrationBoundary_import;	% (2D matrix: num_fluxes_BoundaryImport X [(destiny box) (source box) (destiny box address)]);
intraODEinput.looky_boundary_export     = ODEinput.looky_MigrationBoundary_export;	% (2D matrix: num_fluxes_BoundaryExport X [(destiny box) (source box) (source box address)]);
intraODEinput.looky_driver           	= looky_driver;                             % row address(es) of driver group(s) (e.g., NO3)
intraODEinput.repmat_looky_source       = ODEinput.repmat_looky_MIGRATION_source;	% (2D matrix: num_grps X num_fluxes_migration); NOTE: clm 2 in looky_flux
intraODEinput.repmat_looky_destiny      = ODEinput.repmat_looky_MIGRATION_destiny;	% (2D matrix: num_grps X num_fluxes_migration); NOTE: clm 1 in looky_flux
intraODEinput.repmat_grp_row          	= ODEinput.repmat_GrpRow_MIGRATION;         % addressses; (2D matrix: num_grps X num_fluxes_migration)
Fluxes_t                                = f_PhysicalFlux_intraODE_09092019(intraODEinput, ODEinput); % calculate biomass fluxes into and out of each box and across domain boundaries; (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
flux_net_migration_t                    = Fluxes_t.flux_net_t;                      % net biomass flux into(+) or out of(-) each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
% flux_domain_import_driver_migration_t	= Fluxes_t.flux_domain_import_driver_t;     % net biomass flux OF DRIVER GRP across domain boundary INTO each box (all other grps = 0); (mmoles N/m3/d); (2D matrix: num_grps X num_boxes DESTINATION); NOTE: migration import driver could be used for specialized simulations, but deactived by default
% -------------------------------------------------------------------------


% step 7g: calculate net driver input from advection & horizontal mixing --
%          NOTE: vertical mixing & sinking do not provide external driver (e.g. NO3) to any box (FFF could change this in the future)
%          NOTE: migration does not provide an external driver (theoretically the migration term could be added here, but only in specialized simulations)
flux_domain_import_driver_total_t         = flux_domain_import_driver_advection_t + flux_domain_import_driver_HorizMixing_t;	% (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINATION)
% *************************************************************************





% *************************************************************************
% STEP 8: calculate consumption rate matrix Q_cp---------------------------

% step 8a: arena functional responses -------------------------------------
%   form A: Eqtn. 11 in Steele & Ruzicka 2011 -----------------------------
    FunctionalResponse = ((1+FunctionalResponseParams) .* ProductionRatesC_t_repmat) ./ (FunctionalResponseParams .* productionC_initial_repmat + ProductionRatesC_t_repmat); % (matrix aligned with EnergyBudget); (unitless); (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
    FunctionalResponse(isnan(FunctionalResponse)) = 1; % catch NaNs and set to FunctionalResponse = 1 (caused by 0 biomass cells in sub-surface boxes and possibly by dummy variable NaNs); NOTE: only happens if FunctionalResponseParams = 0 and group consumption @ t = 0, OR if initial group consumption = 0 and group consumption @ t = 0
    
%     % QQQ 8/13/2022 when ProductionRatesC_t_repmat is negative, set FunctionalResponse to 0
%     looky_ProductionNegative = find(ProductionRatesC_t_repmat < 0);
%     FunctionalResponse(looky_ProductionNegative) = 0;
    
    
    
    
%   form B: a purely quadratic form ---------------------------------------
% 	FunctionalRelation_GROUP    = ProductionGROUP_t    / production_initial(current_GROUP);
% 	FunctionalRelation_PREDATOR = ProductionPREDATOR_t / production_initial(current_PREDATOR);

% QQQ turn off specific flow pathways for debugging
% FunctionalResponse(looky_ANYconsumer, :, :)                                           = 0; % QQQ turn off grazing for de-bugging
% FunctionalResponse([looky_terminalPLGCdetritus; looky_terminalBNTHdetritus], :, :)	= 0; % QQQ turn off scenescence for de-bugging
% FunctionalResponse([looky_plgcNH4; looky_bnthNH4], :, :)                              = 0; % QQQ turn off metabolism flow to nutrients
% FunctionalResponse([looky_eggs], :, :)                                                = 0; % QQQ turn off eggs for de-bugging
% FunctionalResponse([looky_fleets], :, :)                                              = 0; % QQQ turn off fisheries for de-bugging
% FunctionalResponse(looky_ANYconsumer, looky_ANYPrimaryProducer, :)                    = 0; % QQQ turn off grazing on phyto for de-bugging
% FunctionalResponse(looky_ANYconsumer, looky_ANYPrimaryProducer, :)                    = FunctionalResponse(looky_ANYconsumer, looky_ANYPrimaryProducer, :) * 0.8; % QQQ turn down grazing on phyto for de-bugging
% FunctionalResponse([looky_plgcNH4; looky_bnthNH4], [1:101], :)                        = 0; % QQQ turn off metabolism flow to nutrients EXCEPT FOR DETRITUS
% FunctionalResponse([looky_terminalPLGCdetritus; looky_terminalBNTHdetritus], [6:end], :) = 0; % QQQ turn off scenescence (except phyto) for de-bugging
% EnergyBudget([looky_plgcNH4 looky_bnthNH4], :, :)                                     = 0; % QQQ turn off all NH4
% -------------------------------------------------------------------------


% step 8b: standard calculation using ECOTRAN EnergyBudget (Acp) ----------
%          NOTE: Q_cp is functional on its own, with or without the Michaelis-Menton option of step 7c
Q_cp        = EnergyBudget .* FunctionalResponse .* ProductionRatesP_t_repmat; % (mmole N/m3/d); (3D matrix: num_grps CONSUMER X num_grps PRODUCER X num_boxes)
% -------------------------------------------------------------------------


% % step 8c: Michaelis-Menton uptake kinetics for phytoplankton -------------
% %          NOTE: For use when driving primary production with Michaelis-Menten
% %          NOTE: deactivate all step 8c lines when driving primary production using EnergyBudget
% %          NOTE: PB (& QB) for nutrients = 1, so these rates (ProductionRatesP_t_repmat) are same as concentrations
% %          NOTE: these nutrient concentration variables are used only in Michaelis-Menten functional response calcs
% NO3_t               = ProductionRatesP_t_repmat(looky_ANYPrimaryProducer, looky_NO3,     1:num_boxes); % (mmole N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes) repeated down rows;
% plgcNH4_t           = ProductionRatesP_t_repmat(looky_ANYPrimaryProducer, looky_plgcNH4, 1:num_boxes); % (mmole N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes) repeated down rows;
% bnthNH4_t           = ProductionRatesP_t_repmat(looky_ANYPrimaryProducer, looky_bnthNH4, 1:num_boxes); % (mmole N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes) repeated down rows;
% NH4_t               = plgcNH4_t + bnthNH4_t;                                % (mmole N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes) repeated down rows;
% NH4fraction_bnth    = bnthNH4_t ./ NH4_t;                                   % value between 0 - 1; (3D matrix: num_PrimaryProducers X 1 X num_boxes) repeated down rows;
% NH4fraction_bnth(isnan(NH4fraction_bnth)) = 0;                              % change div/0 errors to 0
% 
% [Q_cp_NO3, Q_cp_NH4, Q_cp_plgcNH4, Q_cp_bnthNH4, PB_MichaelisMenten] = ...
%     f_MichaelisMenten_05152016(ODEinput, biomass_PrimaryProducer_t, NO3_t, NH4_t, NH4fraction_bnth, Io_t, MLD_t);
% 
% Q_cp(looky_ANYPrimaryProducer, looky_NO3, :)     = Q_cp_NO3;                  % paste in Michaelis-Menten         NO3 uptake by primary producers; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes)
% Q_cp(looky_ANYPrimaryProducer, looky_plgcNH4, :) = Q_cp_plgcNH4;              % paste in Michaelis-Menten pelagic NH4 uptake by primary producers; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes)
% Q_cp(looky_ANYPrimaryProducer, looky_bnthNH4, :) = Q_cp_bnthNH4;              % paste in Michaelis-Menten benthic NH4 uptake by primary producers; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes)
% 
% % if macroalgae are present, scale their producion to the total phytoplankton production rate
% %          NOTE: deactivate these step 8c lines when no macroalgae are present in the model
% Q_cp(looky_macroalgae, looky_NO3, :)          = repmat(sum(Q_cp_NO3), [num_macroalgae, 1])     .* ProductionFraction_macroalgae; % (3D matrix: num_grps X num_grps X num_boxes)
% Q_cp(looky_macroalgae, looky_plgcNH4, :)      = repmat(sum(Q_cp_plgcNH4), [num_macroalgae, 1]) .* ProductionFraction_macroalgae; % (3D matrix: num_grps X num_grps X num_boxes)
% Q_cp(looky_macroalgae, looky_bnthNH4, :)      = repmat(sum(Q_cp_bnthNH4), [num_macroalgae, 1]) .* ProductionFraction_macroalgae; % (3D matrix: num_grps X num_grps X num_boxes)
% *************************************************************************





% *************************************************************************
% STEP 9: calculate consumption_IN & predation_OUT for each group @ t------ 
%         (t/km2/y2) or (mmoles N/m3/d)

% step 9a: calculate consumption_IN for each group @ t --------------------
%          NOTE: predation_OUT = ProductionRates_t' if each clm of EnergyBudget sums to 1 (when FunctionalResponse = 1)
%          NOTE: predation_OUT includes metabolism, eggs, senescence, & feces as well as predation & fleet catch
consumption_IN                      = sum(Q_cp, 2);                                     	% total consumption flow INTO each group, nutrient pool, egg pool, or detritus pool C; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
consumption_IN                    	= reshape(consumption_IN,  [1, num_grps, num_boxes]);   % transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 9b: calculate predation_OUT predation_OUT for each group @ t -------
predation_OUT                     	= sum(Q_cp, 1);                                     	% total predation flow out from each group, nutrient pool, egg pool, or detritus pool P; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 9c: QQQ for 2D ONLY; move production of benthic NH4 from surface boxes (II, IV) to sub-surface boxes (III, V)
%          NOTE: this only allows for 1 benthic NH4 group
NH4_IN                              = consumption_IN(:, looky_bnthNH4, :);                  % (mmole N/m3/d); (3D matrix: 1 X num_bnthNH4 (1) X num_boxes SOURCE)
NH4_BoxVolume_t                    	= reshape(BoxVolume_t, [1, 1, num_boxes]);              % volume of source box; (m3); (3D matrix: 1 X 1 X num_boxes SOURCE); NOTE: this only allows for 1 benthic NH4 pool
NH4_IN                              = NH4_IN .* NH4_BoxVolume_t;                         	% absolute amount of NH4 in source box; (mmole N/d); (3D matrix: 1 X num_bnthNH4 (1) X num_boxes SOURCE)
NH4_IN_repmat                       = repmat(NH4_IN, [num_boxes, 1, 1]);                	% (mmole N/d); (3D matrix: num_boxes DESTINY X num_bnthNH4 (1) X num_boxes SOURCE)
NH4_IN_repmat                       = NH4_IN_repmat .* SinkLink_surface;                    % (mmole N/d); (3D matrix: num_boxes DESTINY X num_bnthNH4 (1) X num_boxes SOURCE)
NH4_IN_repmat                       = sum(NH4_IN_repmat, 3);                              	% net transfer to destination boxes; (mmole N/d); (2D matrix: num_boxes DESTINY X 1)
NH4_BoxVolume_t                   	= squeeze(BoxVolume_t);                                 % volume of destination box; (m3); (2D matrix: 1 X num_boxes DESTINY); (if error because of multiple benthic NH4 pools, add transpose ' to BoxVolume_t([looky_bnthNH4], 1, :)')
NH4_BoxVolume_t                     = NH4_BoxVolume_t';                                     % volume of destination box; (m3); (2D matrix: num_boxes DESTINY X 1);
NH4_transfer                     	= NH4_IN_repmat ./ NH4_BoxVolume_t;                   	% normalize by destination box volume; NH4 transfer from surface to sub-surface boxes; (mmole N/m3/d); (2D matrix: num_boxes DESTINY X 1)
NH4_transfer                      	= reshape(NH4_transfer', [1, 1, num_boxes]);        	% (mmole N/m3/d); (3D matrix: 1 X 1 X num_boxes SOURCE)
consumption_IN(1, looky_bnthNH4, :) = consumption_IN(1, looky_bnthNH4, :) + NH4_transfer;	% add NH4_transfer to Consumption_INflow to transfer benthic NH4 from surface to sub-surface boxes;  (mmole N/m3/d); (3D matrix: 1 X num_grps X source box)


% step 9b: pool pelagic NH4 & benthic NH4 in all boxes --------------------
%          NOTE: ECOTRAN model defines separate NH4 pool fates, but in this
%                spatially-resolved format, the NH4 pool definition should be
%                defined by the individual spatial box being considered.
consumption_IN(1, looky_plgcNH4, :)	= consumption_IN(1, looky_plgcNH4, :) + consumption_IN(1, looky_bnthNH4, :); % (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
consumption_IN(1, looky_bnthNH4, :) = 0;        % zero out benthic NH4 pool after combining with the pelagic NH4 pool; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
% *************************************************************************





% *************************************************************************
% STEP 10: calculate dy----------------------------------------------------                 

% step 10a: ---------------------------------------------------------------
%           NOTE: dy units (mmole N/m3/d2); (3D matrix: 1 X num_grps X num_boxes)
%           NOTE: use pb_t (weight-specific production rate)
%           NOTE: physical loss terms here correspond to John's original format
%                     where I use the term: NetAdvection_t + NetVerticalMix_t + NetHorizontalMix_t + NetSinking_t
%                     John used the term: PhysicalLossFraction * currentProductionRate
%                     These are the same in all essentials except that the
%                     new code also allows for physical GAIN as well as LOSS

dy      = pb_t .* ...
          ((TransferEfficiency .* (consumption_IN + flux_domain_import_driver_total_t)) - ...
          predation_OUT - ...
          ((ConsumptionBudget_em_t - ConsumptionBudget_ba_t) .* ProductionRates_t_transpose) + ...
          flux_net_advection_t + flux_net_HorizMixing_t + flux_net_VertMixing_t + flux_net_sinking_t + flux_net_migration_t);

% dy_evaluate = sum(sum(squeeze(dy))); % QQQ for debugging
% -------------------------------------------------------------------------


% step 10b: stack dy into 1 column vector ---------------------------------
%          NOTE: dy MUST be returned as a column vector, because that's the way MATLAB wants it;
dy                  = dy(:);                    % (mmole N/m3/d2), (vertical vector: (num_grps * num_boxes) X 1)
% dy                  = round2(dy, dy_precision); % clean out any tiny trace rounding error; NOTE: deactive by default
% *************************************************************************


% end ODE function*********************************************************