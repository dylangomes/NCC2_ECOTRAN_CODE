function [Q_cp_NO3, Q_cp_NH4, Q_cp_plgcNH4, Q_cp_bnthNH4, PB_MichaelisMenten] = f_MichaelisMenten_05152016(ODEinput, Biomass_PrimaryProducer_t, NO3_t, NH4_t, NH4fraction_bnth, Io_t, MLD_t)
% phytoplankton uptake rate of NO3 & NH4 (mmole N/m3/d) & PB @ t 
% as calculated from MichaelisMenten uptake (1/y); (2D matrix: num_PrimaryProducers X num_boxes)
% this function must go into MichaelisMenton loop
% takes:
%       ODEinput
%           Kw                          Light attenuation_seawater                              (scaler)
%           Kp                          Light attenuation_phytoplankton     (m2/mmol N)         (scaler)
%           E2E_rows
%           looky_PrimaryProducer
%           num_PrimaryProducers
%       Biomass_PrimaryProducer_t                                           (mmole N/m3)        (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%       NO3_t                             NO3 concentration @ t             (mmoles NO3/m3)     (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%       NH4_t                             NH4 concentration @ t             (mmoles NH4/m3)     (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%       NH4fraction_bnth                                                    (0 - 1)             (3D matrix: num_PrimaryProducers X 1 X num_boxes)
%       Io_t                              surface PAR light intensity @ t; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (scaler)
%       MLD_t                             mixed layer depth @ t             (m)                 (scaler)



% returns:
%       Q_cp_NO3                                                            (mmole N/m3/d)      (2D matrix: num_PrimaryProducers X num_boxes)
%       Q_cp_NH4                                                            (mmole N/m3/d)      (2D matrix: num_PrimaryProducers X num_boxes)
%       PB_MichaelisMenten                                                  (1/y)
% calls:
%       none
% revision date: 5-15-2016


% *************************************************************************
% STEP 1: unpack parameters------------------------------------------------
% step 1a: ODEinput parameters --------------------------------------------
Kw                          = ODEinput.Kw; % Light attenuation_seawater; (scaler)
Kp                          = ODEinput.Kp; % Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N)
looky_PrimaryProducer       = ODEinput.looky_PrimaryProducer;
num_PrimaryProducers        = ODEinput.num_PrimaryProducers;

% step 1b: Michaelis-Menton uptake kinetics for phytoplankton -------------
MichaelisMenten_Vmax        = ODEinput.MichaelisMenten_Vmax;  % Vmax  = maximum nutrient uptake rate;          (1/d);        (3D matrix: num_PrimaryProducers X 1 X num_boxes)
MichaelisMenten_KNO3        = ODEinput.MichaelisMenten_KNO3;  % KNO3  = NO3 half-saturation constant;          (mmol N/m3);  (3D matrix: num_PrimaryProducers X 1 X num_boxes)
MichaelisMenten_KNH4        = ODEinput.MichaelisMenten_KNH4;  % KNH4  = NH4 half-saturation constant;          (mmol N/m3);  (3D matrix: num_PrimaryProducers X 1 X num_boxes)
MichaelisMenten_alpha       = ODEinput.MichaelisMenten_alpha; % alpha = initial slope of light response curve; (m2/W/d);     (3D matrix: num_PrimaryProducers X 1 X num_boxes)
MichaelisMenten_psi         = ODEinput.MichaelisMenten_psi;   % psi   = NO3 uptake inhibition by NH4;          (m3/mmole N); (3D matrix: num_PrimaryProducers X 1 X num_boxes)

% % select only the primary producer rows
% MichaelisMenten_Vmax        = MichaelisMenten_Vmax(looky_PrimaryProducer, 1, :);  % Vmax  = maximum nutrient uptake rate;          (1/d);        (3D matrix: num_PrimaryProducers X 1 X num_boxes)
% MichaelisMenten_KNO3        = MichaelisMenten_KNO3(looky_PrimaryProducer, 1, :);  % KNO3  = NO3 half-saturation constant;          (mmol N/m3);  (3D matrix: num_PrimaryProducers X 1 X num_boxes)
% MichaelisMenten_KNH4        = MichaelisMenten_KNH4(looky_PrimaryProducer, 1, :);  % KNH4  = NH4 half-saturation constant;          (mmol N/m3);  (3D matrix: num_PrimaryProducers X 1 X num_boxes)
% MichaelisMenten_alpha       = MichaelisMenten_alpha(looky_PrimaryProducer, 1, :); % alpha = initial slope of light response curve; (m2/W/d);     (3D matrix: num_PrimaryProducers X 1 X num_boxes)
% MichaelisMenten_psi         = MichaelisMenten_psi(looky_PrimaryProducer, 1, :);   % psi   = NO3 uptake inhibition by NH4;          (m3/mmole N); (3D matrix: num_PrimaryProducers X 1 X num_boxes)

% % squeeze layers into clms to form into a 2D matrix
% MichaelisMenten_Vmax        = squeeze(MichaelisMenten_Vmax);                      % Vmax  = maximum nutrient uptake rate;          (1/d);        (2D matrix: num_PrimaryProducers X num_boxes)
% MichaelisMenten_KNO3        = squeeze(MichaelisMenten_KNO3);                      % KNO3  = NO3 half-saturation constant;          (mmol N/m3);  (2D matrix: num_PrimaryProducers X num_boxes)
% MichaelisMenten_KNH4        = squeeze(MichaelisMenten_KNH4);                      % KNH4  = NH4 half-saturation constant;          (mmol N/m3);  (2D matrix: num_PrimaryProducers X num_boxes)
% MichaelisMenten_alpha       = squeeze(MichaelisMenten_alpha);                     % alpha = initial slope of light response curve; (m2/W/d);     (2D matrix: num_PrimaryProducers X num_boxes)
% MichaelisMenten_psi         = squeeze(MichaelisMenten_psi);                       % psi   = NO3 uptake inhibition by NH4;          (m3/mmole N); (2D matrix: num_PrimaryProducers X num_boxes)

% step 1c: ----------------------------------------------------------------
NH4fraction_plgc            = 1 - NH4fraction_bnth;                                 % (0 - 1); (3D matrix: num_PrimaryProducers X 1 X num_boxes)





% *************************************************************************
% STEP 2: nutrient uptake rates--------------------------------------------
% step 2a: light attenuation of seawater + plankton (+ macroalgae) --------
Ktot                        = (Kw + (Kp .* sum(Biomass_PrimaryProducer_t, 1))) * MLD_t;     % (unitless) units cancel; (3D matrix: 1 X 1 X num_boxes)
Ktot_repmat                 = repmat(Ktot, [num_PrimaryProducers, 1]);                      % replicate down rows;  (unitless); (3D matrix: num_PrimaryProducers X 1 X num_boxes)

% step 2b: ----------------------------------------------------------------
C                           = MichaelisMenten_Vmax ./ (MichaelisMenten_alpha * Io_t);       %                                                        (unitless); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
t1                          = 1 + sqrt(1 + C.^2);                                           %                                                        (unitless); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
t2                          = exp(-Ktot_repmat) + sqrt(exp(-2 * Ktot_repmat) + C.^2);       %                                                        (unitless); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
G                           = (MichaelisMenten_Vmax ./ Ktot_repmat) .* log(t1 ./ t2);       % depth-averaged phyto. growth rate;                     (1/d);      (3D matrix: num_PrimaryProducers X 1 X num_boxes)

% step 2c: ----------------------------------------------------------------
MichaelisMenten_NO3         = (NO3_t ./ (MichaelisMenten_KNO3 + NO3_t) .* exp(-MichaelisMenten_psi .* NH4_t)); % (unitless); (3D matrix: num_E2Egroups X 1 X num_boxes)
MichaelisMenten_NH4         = (NH4_t ./ (MichaelisMenten_KNH4 + NH4_t));                                       % (unitless); (3D matrix: num_E2Egroups X 1 X num_boxes)

% step 2d: ----------------------------------------------------------------
FunctionalResponse_NO3      = G .* MichaelisMenten_NO3;                         % (1/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
FunctionalResponse_NH4      = G .* MichaelisMenten_NH4;                         % (1/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
FunctionalResponse_NO3(isnan(FunctionalResponse_NO3)) = 0;                      % set NaNs to 0 (in sub-surface boxes where Michaelis-Menten parameters are defined as = 0); (1/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
FunctionalResponse_NH4(isnan(FunctionalResponse_NH4)) = 0;                      % set NaNs to 0 (in sub-surface boxes where Michaelis-Menten parameters are defined as = 0); (1/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)

% step 2e: ----------------------------------------------------------------
%          NOTE that this is a function of consumer (PrimaryProducer) biomass
Q_cp_NO3                    = FunctionalResponse_NO3 .* Biomass_PrimaryProducer_t;  % (mmole N/m3/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
Q_cp_NH4                    = FunctionalResponse_NH4 .* Biomass_PrimaryProducer_t;  % (mmole N/m3/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
Q_cp_plgcNH4                = Q_cp_NH4 .* NH4fraction_plgc;                         % divide Q_cp_NH4 between benthic and pelagic NH4 pools; (mmole N/m3/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)
Q_cp_bnthNH4                = Q_cp_NH4 .* NH4fraction_bnth;                         % divide Q_cp_NH4 between benthic and pelagic NH4 pools; (mmole N/m3/d); (3D matrix: num_PrimaryProducers X 1 X num_boxes)



% Q_cp_NO3(Macroalgae_rowclm, 1, :) = sum(Q_cp_NO3([LargePhyto_rowclm; SmlPhyto_rowclm], 1, :)) * ProductionFraction_macroalgae; % QQQ special case for CGoA (nutrient uptake of macroalgae is a set proporation of total phytoplankton uptake)
% Q_cp_NH4(Macroalgae_rowclm, 1, :) = sum(Q_cp_NH4([LargePhyto_rowclm; SmlPhyto_rowclm], 1, :)) * ProductionFraction_macroalgae; % QQQ special case for CGoA (nutrient uptake of macroalgae is a set proporation of total phytoplankton uptake) 


% *************************************************************************
% STEP 3: calculate PB for checking----------------------------------------
P_phyto                                       = Q_cp_NO3 + Q_cp_NH4; % phytoplankton production rate; (mmole N/m3/d); (2D matrix: num_PrimaryProducers X num_boxes)
PB_MichaelisMenten                            = (P_phyto ./ Biomass_PrimaryProducer_t) * 365; % (1/y)
PB_MichaelisMenten(isnan(PB_MichaelisMenten)) = 0; % filter out NaNs from times & spaces with 0 primary producer biomass

% end m-file***************************************************************