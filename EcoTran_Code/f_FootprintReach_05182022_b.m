function [DIET_MC, FOOTPRINT_array_MC, REACH_array_MC, FOOTPRINT_trace_MC, REACH_trace_MC, FOOTPRINT_system_MC, REACH_system_MC] = ...
    f_FootprintReach_05182022_b(ECOTRAN, ECOTRAN_MC, PhysicalLossFraction, TraceGroup)
%
% Re-calculate the DIET matrix from the ECOTRAN EnergyBudget & ConsumptionBudget
%
% Calculate "FOOTPRINT" metrics of each functional group.
% The FOOTPRINT metrics are the fraction of each PRODUCER group's production
%       flowing to CONSUMER = TraceGroup. Code calculates FOOTPRINT for
%       all functional groups as TraceGroup.
%
% Calculate the "REACH" metrics of each functional group.
% The REACH metrics are the fraction of each CONSUMER group's production ultimately 
%       originating from PRODUCER = TraceGroup. Code calculates REACH for
%       all functional groups as TraceGroup.
%
% Calculate the FOOTPRINT and REACH traces for each trophic linkage for 
%       one (1) specific functional group of interest = TraceGroup.
% FOOTPRINT_trace is the fraction of each trophic link ultimately contributing
%       to production of CONSUMER = TraceGroup
%       It is the relative contribution of each linkage to the production of 
%       CONSUMER = TraceGroup.
% REACH_trace is the fraction of each trophic link ultimately originating
%       from PRODUCER = TraceGroup.
%       It is the fraction of energy within each linkge ultimately 
%       originating from PRODUCER = TraceGroup.
%
% by Jim Ruzicka
%
% calls:
%   f_WebProductivity_03272019      - calculate consumption rate (q) for 
%                                     each group given a group driver 
%                                     production rate
%   f_DietTrace_08032020            - use DIET matrix to calculate the 
%                                     "REACH" metrics of each functional
%                                     group, and calculate the FOOTPRINT & 
%                                     REACH traces for each trophic linkage
%                                     for one (1) specific functional group
%                                     of interest = TraceGroup.
%    
% takes:
%       ECOTRAN
%           num_grps
%           GroupType
%           ee
%           GroupTypeDef_
%       ECOTRAN_MC
%           num_MC;
%           EnergyBudget_MC             (3D matrix: num_grps (consumers) X num_grps (producers) X num_MC)
%           ConsumptionBudget_MC        (3D matrix: 7 X num_grps (producers) X num_MC)
%           TransferEfficiency_MC       (2D matrix: num_MC X num_grps)
%       PhysicalLossFraction
%       TraceGroup
%
% returns:
%       DIET_MC         	- diet of each consumer
%                       	- (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%
%       FOOTPRINT_array_MC	- "FOOTPRINT" of each TraceGroup = consumer (rows) upon each producer (columns)
%                        	- the fraction of each producer group's production flowing to each consumer = TraceGroup
%                         	- (3D matrix: num_grps (consumer) X num_grps (producer) X num_MC)
%                         	- NOTE: use for web plotting footprint box colors relative to TraceGroup
%                         	- NOTE: does NOT include import diet.
%                                   There is no footprint calculated for producer = import
%                                   because the value of import production can be arbitrary.
%
%       REACH_array_MC      - "REACH" of TraceGroup = producer (rows) to each consumer (columns)
%                           - the fraction of each consumer group's production ultimately originating from producer = TraceGroup                 
%                           - (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%                           - NOTE: use to plot REACH box colors in food web diagram
%                           - NOTE: formerly called "TraceFraction_upward"
%                           - NOTE: DOES include import diet
%
%       FOOTPRINT_trace_MC	- Fraction of each trophic link ultimately contributing to production of CONSUMER = TraceGroup
%                           - the relative contribution of each linkage to production of CONSUMER = TraceGroup
%                           - (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%                           - NOTE: use to plot FOOTPRINT arrow colors in food web diagram
%                           - NOTE: formerly called "DietTrace_downward"
%                           - NOTE: DOES include import diet
%
%       REACH_trace_MC      - Fraction of each trophic link ultimately originating from PRODUCER = TraceGroup
%                           - It is the fraction of energy within each linkge ultimately originating from PRODUCER = TraceGroup.
%                           - (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC)
%                           - NOTE: use to plot REACH arrow colors in food web diagram
%                           - NOTE: formerly called "DietTrace_upward"
%                           - NOTE: DOES include import diet
%
%       FOOTPRINT_system_MC	- SYSTEM-LEVEL footprint = ratio of total TraceGroup footprint on all PRODUCERS over total production of all CONSUMER groups
%                           - total consumer production of ecosystem excludes microzooplankton
%                           - (vertical vector: num_MC X 1)
%
%       REACH_system_MC     - SYSTEM-LEVEL reach = ratio of total TraceGroup production going to all consumers over total production of all CONSUMER groups
%                           - total consumer production of ecosystem excludes microzooplankton
%                           - (vertical vector: num_MC X 1)
%
% revision date: 5-18-2022


% *************************************************************************
% STEP 1: unpack ECOTRAN info----------------------------------------------
% step 1a: ----------------------------------------------------------------
num_grps                = ECOTRAN.num_grps;
GroupType               = ECOTRAN.GroupType;
ee                      = ECOTRAN.ee;
num_MC                  = ECOTRAN_MC.num_MC;
EnergyBudget_MC         = ECOTRAN_MC.EnergyBudget_MC; % (3D matrix: num_grps (consumers) X num_grps (producers) X num_MC)
ConsumptionBudget_MC	= ECOTRAN_MC.ConsumptionBudget_MC; % (3D matrix: 7 X num_grps (producers) X num_MC)
TransferEfficiency_MC	= ECOTRAN_MC.TransferEfficiency_MC; % (2D matrix: num_MC X num_grps)
% -------------------------------------------------------------------------


% step 2b: find group addresses -------------------------------------------
looky_NO3                   = find(GroupType        == ECOTRAN.GroupTypeDef_NO3);
looky_NH4                   = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4 | GroupType == ECOTRAN.GroupTypeDef_bnthNH4);
looky_nutrients             = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYNitroNutr);
looky_detritus              = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYDetritus);
looky_terminalBnthDetr      = find(GroupType == ECOTRAN.GroupTypeDef_terminalBnthDetr); % address of the single terminal benthic detritus group
looky_ANYPrimaryProducer	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYPrimaryProd);
looky_ANYconsumer           = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYConsumer);
looky_fleets                = find(floor(GroupType) == ECOTRAN.GroupTypeDef_fleet);
looky_consumerANDfleets     = [looky_ANYconsumer; looky_fleets];
looky_micrograzers          = find(GroupType        == ECOTRAN.GroupTypeDef_micrograzers);
looky_bacteria              = find(GroupType        == ECOTRAN.GroupTypeDef_bacteria);
% -------------------------------------------------------------------------


% step 2c: define Transfer Efficiency of terminal benthic detritus --------
%          use the "type" model ee of terminalBnthDetr for TransferEfficiency
%          FFF in future, could calculate ee of individual MC models?
ee_terminalBnthDetr     = ee(looky_terminalBnthDetr);
% *************************************************************************





% *************************************************************************
% *************************************************************************
% process each Monte Carlo model

% initialize result variables ---------------------------------------------
DIET_MC                 = zeros((num_grps+1), (num_grps+1), num_MC); % initialize; (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC); NOTE: DOES include import diet
FOOTPRINT_array_MC      = zeros(num_grps, num_grps, num_MC); % initialize; (3D matrix: num_grps (consumers) X num_grps (producers) X num_MC); NOTE: does NOT include import diet
REACH_array_MC          = zeros((num_grps+1), (num_grps+1), num_MC); % (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC); NOTE: DOES include import diet
FOOTPRINT_trace_MC      = zeros((num_grps+1), (num_grps+1), num_MC); % (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC); NOTE: DOES include import diet
REACH_trace_MC          = zeros((num_grps+1), (num_grps+1), num_MC); % (3D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import) X num_MC); NOTE: DOES include import diet
FOOTPRINT_system_MC     = zeros(num_MC, 1); % initialize; (vertical vector: num_MC X 1); NOTE: does NOT include import diet
REACH_system_MC         = zeros(num_MC, 1); % initialize; (vertical vector: num_MC X 1); NOTE: DOES include import diet
% -------------------------------------------------------------------------


for MC_loop = 1:num_MC

%     display(['Processing MC model ' num2str(MC_loop) ' of ' num2str(num_MC)])

    % *********************************************************************
    % STEP 3: select current model-----------------------------------------
    current_EnergyBudget        = EnergyBudget_MC(:, :, MC_loop); % (2D matrix: num_grps (consumers) X num_grps (producers))
    current_EnergyBudget(isnan(current_EnergyBudget)) = 0; % QQQ temp patch for ECOBASE models, I need to filter out NaN in dummy terminal detritus groups
    current_ConsumptionBudget	= ConsumptionBudget_MC(:, :, MC_loop); % (2D matrix: 7 X num_grps (producers))
    current_TransferEfficiency	= TransferEfficiency_MC(MC_loop, :); % (horizontal vector: 1 X num_grps)
    current_cb_em               = current_ConsumptionBudget(7, :); % current fraction of consumption budget due to emigration; (horizontal vector: 1 X num_grps)
    current_pq                  = sum(current_ConsumptionBudget([3 4 5], :)); % current fraction of consumption budget going to production (cb_eggs + cb_predation + cb_senescence; NOTE: not using cb_em); (horizontal vector: 1 X num_grps)
    current_pq(looky_fleets)	= current_ConsumptionBudget(7, looky_fleets); % add in fleet landings fraction (ee_em) as fleet qb value
    % *********************************************************************
   
   
    
    
    
	% *********************************************************************
    % STEP 4: calculate footprint for all groups in current MC model-------
    %         FOOTPRINT_array	= footprint of each consumer upon each producer
    %                           = the fraction of each producer consumption (q) supporting (flowing to) each consumer
    %                           (2D matrix: num_grps (consumers) X num_grps (producers)
    %         NOTE: there is no footprint calculated for producer = import
    %               because the value of import production can be arbitrary
    %         FFF in future also provide as fraction of each producer's production (in addition to consumption)
    
    % step 4a: shut off recycling & boost TE on terminal benthic detritus -
    %          zero-out nutrients & detritus as consumer rows in the current_EnergyBudget
    %          NOTE: inclusion of nutrient and detritus recycling will increase the footprints substantially as a greater fraction of each group's production gets used and re-used to support the food web.
    %          NOTE: footprints for nutrients & detritus are still calculated, even with recycling deactivated, because we are driving the model for each group separately by providing a defined producion rate that is independent of input from the rest of the food web.
	working_EnergyBudget                                        = current_EnergyBudget;
    working_EnergyBudget([looky_nutrients; looky_detritus], :)	= 0;
	current_TransferEfficiency(looky_terminalBnthDetr)          = 0.9999; % footprints on terminal detritus will be underestimated by 1-TE = 0.01%.
    % ---------------------------------------------------------------------
    
    
	% step 4b: calculate FOOTPRINT_array ----------------------------------
    %          Drive the model with 1 unit of each functional group, individually
	current_FOOTPRINT_array             = zeros(num_grps, num_grps); % initialize current_FOOTPRINT_array; (2D matrix: num_grps (consumer) X num_grps (producer))
    
    for grp_loop = 1:num_grps
        current_grp                             = grp_loop; % current group as producer
        InputProductionVector                   = zeros(1, num_grps); % re-initialize for each loop iteration; (horizontal vector: 1 X num_grps)
        InputProductionVector(current_grp)      = 1; % drive model with 1 unit of current_grp; (horizontal vector: 1 X num_grps)
        q_current_grp                           = f_WebProductivity_03272019(current_TransferEfficiency, working_EnergyBudget, InputProductionVector, PhysicalLossFraction); % production of each model group (rows) when driven by each current_grp (clms); (2D matrix: num_grps X num_grps)
        q_current_grp(current_grp)           	= 0; % set footprint of current_grp on self to ZERO
        current_FOOTPRINT_array(:, grp_loop)	= q_current_grp; % (2D matrix: num_grps (consumers) X num_grps (producers))
    end % (grp_loop)
    
    FOOTPRINT_array_MC(:, :, MC_loop)	= current_FOOTPRINT_array; % fraction of each producer group's production flowing to each consumer; (3D matrix: num_grps (consumer) X num_grps (producer) X num_MC)
    % *********************************************************************

    

    
    
	% *********************************************************************
    % STEP 5: calculate ECOTRAN diet matrix for current MC model-----------
    %         NOTE: this is not exactly the same as the EwE diet because:
    %               1) cannibalism has been removed,
    %               2) ECOTRAN modifies EwE detritus fates and separates feces & senescence
    %         NOTE: largest difference from EwE diet will most likely be due to detritus contributions to group diets
	
    % step 5a: shut off nutrient recycling --------------------------------
    %          zero-out NH4 as producer columns in the current_EnergyBudget
    %          NOTE: inclusion of nutrient recycling will increase the footprints substantially as a greater fraction of each group's production gets used and re-used to support the food web.
    %          NOTE: footprints for NH4 pools are still calculated, even with recycling deactivated, because we are driving the model for each group separately by providing a defined producion rate that is independent of input from the rest of the food web.
	working_EnergyBudget                     	= current_EnergyBudget;
    working_EnergyBudget(:, [looky_NH4])     	= 0; % turn off NH4 uptake
    working_EnergyBudget([looky_NH4], [looky_detritus])	= 0; % turn off implicit bacterial metabolism of detritus
    sum_detritus                                = sum(working_EnergyBudget(:, [looky_detritus]));
    sum_detritus                                = repmat(sum_detritus, [num_grps, 1]);
    working_EnergyBudget(:, [looky_detritus])	= working_EnergyBudget(:, [looky_detritus]) ./ sum_detritus; % re-normalize detritus columns in current_EnergyBudget to sum to 1
    working_EnergyBudget(isnan(working_EnergyBudget)) = 0; % correct for div/0 errors
    % ---------------------------------------------------------------------
    
    
	% step 5b: define TE of terminal benthic detritus ---------------------
    %          NOTE: Using "type" model ee of terminalBnthDetr for TransferEfficiency for all MC models
    %          FFF could calculate ee of individual MC models?
    if ee_terminalBnthDetr == 0
        ee_terminalBnthDetr = 0.001; % added for SPF ECOBASE models that have dummey terminal benthic detritus groups
    end
    
    current_TransferEfficiency(looky_terminalBnthDetr)	= ee_terminalBnthDetr;
    % ---------------------------------------------------------------------

    
	% step 5c: calculate domestic consumption (q) of each group when driven by NO3 -----
    InputProductionVector               = zeros(1, num_grps); % re-initialize for each loop iteration; (horizontal vector: 1 X num_grps)
    InputProductionVector(looky_NO3)	= 100; % 100 units of NO3 to drive the model
    q_domestic                         	= f_WebProductivity_03272019(current_TransferEfficiency, working_EnergyBudget, InputProductionVector, PhysicalLossFraction); % total consumption q by each model group; (mmoles N/m3/d); (vertical vector: num_grps X 1)
	% ---------------------------------------------------------------------

    
	% step 5d: calculate import & total consumption of each group ---------
    q_total       = q_domestic ./ (1 - ((-1) * (current_cb_em'))); % q_total = q_domestic + q_import; (vertical vector: num_grps X 1); NOTE (-1) multiplication; NOTE transpose
    q_import      = q_total - q_domestic;
    % ---------------------------------------------------------------------
    
    
    % step 5e: calculate CONSUMPTION_matrix (Q) ---------------------------
    q_total_repmat      = repmat(q_total', [num_grps, 1]); % replicate q_total down rows; (2D matrix: num_grps X num_grps (producer)); NOTE transpose of q_total
    CONSUMPTION_matrix	= working_EnergyBudget .* q_total_repmat; % consumption rate of each producer by each consumer; (mmoles N/m3/d); (2D matrix: num_grps (consumer) X num_grps (producer));
    % ---------------------------------------------------------------------

    
    % step 5f: calculate and add import column of CONSUMPTION_matrix ------
    num_grps_plus_import	= num_grps + 1; % add import diet row & columns
    looky_import            = num_grps_plus_import;
	CONSUMPTION_matrix(:, looky_import)             = q_import; % append q_import term (import consumption) as column to end of CONSUMPTION_matrix
    CONSUMPTION_matrix(looky_fleets, looky_import)	= 0; % zero out negative import fate to fleets
    
    CONSUMPTION_matrix(looky_import, 1:num_grps_plus_import)	= 0; % add import row to CONSUMPTION_matrix
    CONSUMPTION_matrix(looky_import, looky_fleets)	= (-1) * q_import(looky_fleets, 1)'; % append landings as flow to "import"; NOTE transpose; NOTE (-1) multiplier
    
    q_total(looky_import, 1) = sum(q_import(looky_fleets, 1)); % append q_import to end of q; q(looky_import, 1) is only the fleet landings, it represents exports from system; FFF could add any non-negative q_import to the sum to include any theoretical export
    % ---------------------------------------------------------------------

    
    % step 5g: calculate diet matrix of current MC model ------------------
    %          divide CONSUMPTION_matrix by consumer's total consumption (q)
    %          NOTE: re-calculating consumption rates from CONSUMPTION_matrix instead of using q_total
    sum_consumption         = sum(CONSUMPTION_matrix, 2); % We could have used q here, but using the horizontal sum of CONSUMPTION_matrix accounts for the TransferEfficiency of terminalBnthDetr
    sum_consumption_repmat	= repmat(sum_consumption, [1, num_grps_plus_import]);       % replicate q across columns; (2D matrix: num_grps X num_grps)
    current_DietArray       = CONSUMPTION_matrix ./ sum_consumption_repmat; % fraction of each producer flowing directly to each consumer; (2D matrix: num_grps (consumers) X num_grps (producers))
	current_DietArray(isnan(current_DietArray)) = 0; % set div/0 errors to 0
    current_DietArray       = current_DietArray';    % (2D matrix: num_grps (producers) X num_grps (consumers))
	DIET_MC(:, :, MC_loop)	= current_DietArray;     % (3D matrix: num_grps (producers) X num_grps (consumers) X num_MC)
    % *********************************************************************
    
    
    
    
    
	% *********************************************************************
    % STEP 6: calculate the "REACH" metrics of each functional group, and
    %         calculate the FOOTPRINT & REACH traces for every trophic linkage
    %         for one (1) specific functional group of interest = TraceGroup.
    %       REACH_vector    - "REACH" of TraceGroup = producer to each consumer (columns)
    %                    	- the fraction of each consumer group's production ultimately originating from producer = TraceGroup                 
    %                    	- (horizontal vector: 1 X num_grps (consumers + 1 import))
    %                   	- NOTE: use to plot REACH box colors in food web diagram
    %                    	- NOTE: formerly called "TraceFraction_upward"
    %                   	- NOTE: DOES include import diet
    %
    %       FOOTPRINT_trace	- Fraction of each trophic link ultimately contributing to production of CONSUMER = TraceGroup
    %                       - the relative contribution of each linkage to production of CONSUMER = TraceGroup
    %                       - (2D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import))
    %                    	- NOTE: use to plot FOOTPRINT arrow colors in food web diagram
    %                     	- NOTE: formerly called "DietTrace_downward"
    %                     	- NOTE: DOES include import diet
    %
    %       REACH_trace     - Fraction of each trophic link ultimately originating from PRODUCER = TraceGroup
    %                       - It is the fraction of energy within each linkge ultimately originating from PRODUCER = TraceGroup.
    %                    	- (2D matrix: num_grps (producers + 1 import) X num_grps (consumers + 1 import))
    %                     	- NOTE: use to plot REACH arrow colors in food web diagram
    %                    	- NOTE: formerly called "DietTrace_upward"
    %                    	- NOTE: DOES include import diet

    % step 6a: calculate REACH_vector values for each model group ---------
    for grp_loop = 1:num_grps_plus_import
        
        current_grp                             = grp_loop;
        
        [current_REACH_vector, ~, ~]            = ...
            f_DietTrace_08032020(ECOTRAN, current_DietArray, current_grp);

        REACH_array_MC(current_grp, :, MC_loop)	= current_REACH_vector;
    
    end % (grp_loop)
	% ---------------------------------------------------------------------

    
    % step 6b: calculate FOOTPRINT_trace & REACH_trace for the single, desired TraceGroup
    %          NOTE: These trace arrays are calculated for each group in
    %                step 6a, but we only need to retain the values for 
    %                plotting food web diagrams centered about the TraceGroup.
    [~, current_FOOTPRINT_trace, current_REACH_trace] = ...
         f_DietTrace_08032020(ECOTRAN, current_DietArray, TraceGroup);

    FOOTPRINT_trace_MC(:, :, MC_loop)	= current_FOOTPRINT_trace;
    REACH_trace_MC(:, :, MC_loop)       = current_REACH_trace;
    % *********************************************************************


    
    
    
    % *********************************************************************
    % STEP 7: calculate SYSTEM-LEVEL FOOTPRINT of CONSUMER = TraceGroup & REACH of PRODUCER = TraceGroup

    % step 7a: calculate production rates of the current MC model ---------
    %          zero-out NH4 as producer columns in the working_EnergyBudget to shut off nutrient recycling
    %          NOTE: inclusion of nutrient recycling will increase the footprints substantially as a greater fraction of each group's production gets used and re-used to support the food web.
	working_EnergyBudget                        = current_EnergyBudget;
    working_EnergyBudget(:, looky_NH4)          = 0; % turn off NH4 uptake by setting NH4 column to ZERO
    working_EnergyBudget(looky_NH4, looky_detritus)	= 0; % turn off implicit bacterial metabolism of detritus
    sum_detritus                                = sum(working_EnergyBudget(:, looky_detritus));
    sum_detritus                                = repmat(sum_detritus, [num_grps, 1]);
    working_EnergyBudget(:, looky_detritus)     = working_EnergyBudget(:, looky_detritus) ./ sum_detritus; % re-normalize detritus columns in working_EnergyBudget to sum to 1
    working_EnergyBudget(isnan(working_EnergyBudget)) = 0; % correct for div/0 errors

    if ee_terminalBnthDetr == 0
        ee_terminalBnthDetr = 0.001; % added for SPF ECOBASE models that have dummy terminal benthic detritus groups
    end
    
    current_TransferEfficiency(looky_terminalBnthDetr)	= ee_terminalBnthDetr;

	InputProductionVector               = zeros(1, num_grps); % initialize; (horizontal vector: 1 X num_grps)
    InputProductionVector(looky_NO3)	= 100; % 100 units of NO3 to drive the model
    q_domestic                         	= f_WebProductivity_03272019(current_TransferEfficiency, working_EnergyBudget, InputProductionVector, PhysicalLossFraction); % total consumption q by each model group; (mmoles N/m3/d); (vertical vector: num_grps X 1)
    p_domestic                          = q_domestic' .* current_pq; % production rates of all groups when driven with 100 units of NO3; (horizontal vector: 1 X num_grps); NOTE transpose

    % production values for FOOTPRINT
    %       NOTE: total_production_PRODUCER excludes phytoplankton, micrograzers, & bacteria production
	production_PRODUCER                 = p_domestic; % (horizontal vector: 1 X num_grps)
	production_PRODUCER([looky_ANYPrimaryProducer; looky_micrograzers; looky_bacteria]) = 0; % exclude phytoplankton, micrograzer & bacteria production by setting to 0; (horizontal vector: 1 X num_grps)
    production_PRODUCER([looky_nutrients; looky_detritus]) = []; % remove nutrient & detritus columns; (horizontal vector: 1 X num_grps (- nutrients & detritus)
    total_production_PRODUCER           = sum(production_PRODUCER); % total production of all PRODUCER groups; (scaler)
    
    % production values for REACH
	total_production_CONSUMER           = p_domestic; % (horizontal vector: 1 X num_grps)
	total_production_CONSUMER(looky_micrograzers)	= 0; % QQQ still need a formal micrograzer definition, need to do this manually for now; zero-out micrograzers QQQ
    total_production_CONSUMER(looky_import)         = 0; % QQQ add import
    total_production_CONSUMER           = total_production_CONSUMER(looky_consumerANDfleets); % only consider consumer production
    total_production_CONSUMER           = sum(total_production_CONSUMER); % total production of all PRODUCER groups (scaler)
    
 	production_CONSUMER                 = p_domestic;
	production_CONSUMER(looky_import)	= 0; % QQQ add import

%     production_CONSUMER([looky_nutrients; looky_detritus]) = []; % remove non-producer group rows (nutrients, detritus, import)
    % ---------------------------------------------------------------------
    
    
    % step 7b: calculate SYSTEM-LEVEL FOOTPRINT of CONSUMER = TraceGroup --
    FOOTPRINT_array_PRODUCER            = current_FOOTPRINT_array; % (2D matrix: num_grps (consumers) X num_grps (producers))
    FOOTPRINT_array_PRODUCER(:, [looky_nutrients; looky_detritus]) = []; % remove nutrient & detritus columns; (2D matrix: num_grps (consumers) X num_grps (producers, minus nutrients & detritus clms))
    
    FOOTPRINT_vector_PRODUCER           = FOOTPRINT_array_PRODUCER(TraceGroup, :); % (horizontal vector: 1 X num_grps (minus nutrients & detritus))
    
    TraceGroup_FootprintProduction      = FOOTPRINT_vector_PRODUCER .* production_PRODUCER; % footprint of CONSUMER = TraceGroup upon each producer in terms of production; FOOTPRINT_vector_PRODUCER * production_PRODUCER = (net or gross production of each PRODUCER going to CONSUMER = TraceGroup); (horizontal vector: 1 X num_grps (minus nutrients & detritus))
    TraceGroup_FootprintProduction_total	= sum(TraceGroup_FootprintProduction); % total CONSUMER = TraceGroup footprint on all PRODUCERS; (scaler)

    current_SystemFootprint             = TraceGroup_FootprintProduction_total ./ total_production_PRODUCER; % CONSUMER = TraceGroup footprint fraction across all producer groups; (scaler)
    FOOTPRINT_system_MC(MC_loop)        = current_SystemFootprint; % (vertical vector: num_MC X 1); NOTE: does NOT include import diet
    % *********************************************************************
    
    
    
	% *********************************************************************
    % STEP 8: calculate SYSTEM-LEVEL REACH of PRODUCER = TraceGroup--------
    
    REACH_vector_TraceGroup             = REACH_array_MC(TraceGroup, :, MC_loop); % (horizontal vector: 1 X num_grps (+ 1 import))
    REACH_vector_TraceGroup(TraceGroup) = 0; % set reach of TraceGroup to self to ZERO
    
    TraceGroup_ReachProduction          = REACH_vector_TraceGroup .* production_CONSUMER; % reach of PRODUCER = TraceGroup to each consumer in terms of production; trace_fraction * (CONSUMER production) = (net PRODUCER production going to each CONSUMER); (horizontal vector: 1 X num_grps (+ 1 import))
    TraceGroup_ReachProduction_total	= sum(TraceGroup_ReachProduction); % total PRODUCER = TraceGroup production going to all consumers; (scaler)
    current_SystemReach                 = TraceGroup_ReachProduction_total ./ total_production_CONSUMER; % PRODUCER = TraceGroup reach fraction across all consumer groups; (scaler)
    REACH_system_MC(MC_loop)            = current_SystemReach; % (vertical vector: num_MC X 1); NOTE: DOES include import diet
	% *********************************************************************

    
end % (MC_loop)
% *************************************************************************
% *************************************************************************


% end m-file***************************************************************