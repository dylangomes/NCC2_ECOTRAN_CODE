function SUBREGION = NCC2_SubregionReadIn_12072022(ECOTRAN,SR) % remove this for non-function version
% force parameters for subregional model
% parameters from file: DepthZoneSubmodels_CNP_11192021.xlsx

FILE = strcat("C:/Users/dgome/Documents/NCC2_ECOTRAN_CODE/NCC_SubRegional/SubregionSubmodels_NCC_10252022_AS_csv/zone ",num2str(SR),".csv");
SubRegion = readmatrix(FILE);

SUBREGION.DepthFoodWeb                        = FILE;
% display(strcat("Read in: ",SUBREGION.DepthFoodWeb))

num_grps = ECOTRAN.num_grps;

SUBREGION.CONSUMPTION               = zeros(num_grps, num_grps);
SUBREGION.DIET                     	= zeros(num_grps, num_grps);
SUBREGION.import_diet              	= zeros(1, num_grps);

SUBREGION.DIET_NoCannibalism       	= [];
SUBREGION.CONSUMPTION_NoCannibalism	= [];

SUBREGION.ConsumptionBudget = SubRegion(235:241,6:107); 
SUBREGION.EnergyBudget      = SubRegion(249:350,6:107); 
        
SUBREGION.fate_metabolism = SubRegion(122:124,6:107); 
SUBREGION.fate_feces      = SubRegion(126:128,6:107);
SUBREGION.fate_senescence = SubRegion(130:132,6:107);    
SUBREGION.fate_eggs       = SubRegion(134:135,6:107);
SUBREGION.fate_predation  = SubRegion(137:230,6:107);

SUBREGION.ee            = SubRegion(360,6:107);
SUBREGION.ee_predation	= SubRegion(356,6:107);

% P/B terms
SUBREGION.pb              = ECOTRAN.pb'; % (1/y); pb for nutrients & detritus is always 1 regardless of time-frame; NOTE transpose
SUBREGION.pb              = SUBREGION.pb * (1/365); % (1/d); pb for nutrients & detritus is always 1 regardless of time-frame; NOTE transpose

% grab nutrients (groups ~1) and detritus (groups ~4)
cond = (ECOTRAN.GroupType<2) | (ECOTRAN.GroupType>=4 & ECOTRAN.GroupType<5);
SUBREGION.pb(cond)=1;% (1/d); pb for nutrients, eggs, & detritus is always 1 regardless of time-frame; NOTE transpose


SUBREGION.qb              =  ECOTRAN.qb'; % (1/y); qb for nutrients, detritus, & fleets is always 1 regardless of time-frame; NOTE transpose
SUBREGION.qb              = SUBREGION.qb * (1/365); % (1/d); qb for nutrients, detritus, & fleets is always 1 regardless of time-frame; NOTE transpose
% grab nutrients (groups ~1) and detritus (groups ~4) and fisheries (group == 5)
cond2 = or(ECOTRAN.GroupType==5 , cond);
SUBREGION.qb(cond2)	= 1; % (1/d); qb for nutrients, eggs, detritus, & fleets is always 1 regardless of time-frame; NOTE transpose


% % production_initial (really consumption_initial)
% %       initial q; (t WWT/km2/d); 
% %       from file DepthZoneSubmodels_CNP_11192021.xlsx (cell D166); 
% %       NOTE transpose to vertical vector
% SUBREGION.production_initial = [1	1650.63159	0	0.55200001	0.98399999	1.176000024	39.38399849	230.5081428	95.25261044	13.13930872	7.620174876	3.569868841	3.259684549	8.049150146	6.846105169	14.37798461	4.153639017	18.3979413	0.823847364	7.999457608	0	28.27951054	15.50170572	0.67731766	1.438113264	1.032414162	0	1.901769705	0.150480123	0.014965675	0.016680316	0	0	0.051134229	0.065009223	0	0	0	0	0	0	0.002097335	0.002523914	0.001950419	0	0	0.277611253	0.21225364	1169.45562	4.24787502	3400.332521	1276.41184	50.46195097	0.001875117	0	0	9.09937E-05	0.000374533]'; 

%       initial q; (mmoles N/m3/d); 
%       from INITIALproduction_nutrients for EPIpelagic (static, before time-dynamic run); see file "Initial_q_rates_06122022.xlsx" tab "Sheet 4"
SUBREGION.production_initial = ones(1,num_grps);

% SUBREGION.production_initial = SUBREGION.production_initial * (1/(1000*1000));	% area to volumetric conversion; (t/m2/y)
% SUBREGION.production_initial = SUBREGION.production_initial * (1/BoxHeight);      % area to volumetric conversion; (t/m3/y)
% SUBREGION.production_initial = SUBREGION.production_initial * 1000000;            % (g/m3/y)
% SUBREGION.production_initial = SUBREGION.production_initial * 1000;               % (mg WWT/m3/y)
% SUBREGION.production_initial = SUBREGION.production_initial * (1/WWT_to_C);       % (mg C/m3/y)
% SUBREGION.production_initial = SUBREGION.production_initial * (1/atomic_mass_C);	% (mmole C/m3/y)
% SUBREGION.production_initial = SUBREGION.production_initial * (1/C_to_N_phytoplankton);	% (mmole N/m3/y)
% SUBREGION.production_initial = SUBREGION.production_initial * (1/365);            % (mmole N/m3/d)

% calculate biomass
SUBREGION.biomass = SUBREGION.production_initial .* (1./SUBREGION.qb); % (mmoles N/m3); initial biomasses (only primary producer biomasses really used)

SUBREGION;

% save(test.mat,SUBREGION)
% end m-file***************************************************************