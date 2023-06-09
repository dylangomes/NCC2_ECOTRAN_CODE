function StaticScenario_results = f_ScenarioGenerator_10212020(ScenarioConditions, ECOTRAN)
% by Jim Ruzicka
% generate a modified scenario matrix (or set of matrices)
% calculate productivity vectors for scenario and ratios compared to un-modified matrices
%
% calls:
%       f_WebProductivity_03272019	calculate production rates of each functional group (vertical vector: num_grps X 1)
%
% takes:
%       ScenarioConditions
%           modify_consumer         row number(s) of consumer group(s) to force-modify
%           ScaleFactor             factor by which to re-scale production to modify_consumer group(s); NOTE: value > 1 means increase flow to modify_consumer (and reduced flow to offset_consumer); < 1 means decrease flow to modify_consumer; NOTE: keep ScaleFactor >= 0 (negative value makes no sense)
%           offset_consumer         row number(s) of consumer group(s) to modify as offset to force-modified grp(s)NOTE: change in modify_consumer grp(s); can be [] if offset is to be distributed among ALL consumer groups
%           trgt_producer           column number(s) of producer group(s) to modify
%               EnergyBudget_base           (NOTE: NOT USED); MonteCarlo set of EnergyBudget_base models; (3D matrix: num_grps X num_grps X num_MC); NOTE: no changes made directly to EnergyBudget_base, everything is done via the fates
%           ConsumptionBudget_base	MonteCarlo set of ConsumptionBudget_base models; (3D matrix: 7 X num_grps X num_MC)
%                                       1) feces
%                                       2) metabolism
%                                       3) eggs (reproduction)
%                                       4) predation
%                                       5) senescence
%                                       6) ba (biomass accumulation)
%                                       7) em (emigration); NOTE: negative for immigration
%           fate_metabolism_base    fate of metabolism in base model; (3D matrix: num_nutrients X num_grps X num_MC)
%           fate_eggs_base          fate of eggs (reproduction) in base model; (3D matrix: num_eggs X num_grps X num_MC)
%           fate_feces_base         fate of feces detritus in base model; (3D matrix: num_ANYdetritus X num_grps X num_MC)
%           fate_predation_base     fate of production among all predators in base model; (3D matrix: num_livingANDfleets X num_grps X num_MC)
%           fate_senescence_base	fate of senescence detritus in base model; (3D matrix: num_ANYdetritus X num_grps X num_MC)
%           DiscardFraction_base	fractions of catch discarded for each group and fleet; (3D matrix: num_grps X num_fleets X num_MC)
%           TransferEfficiency          (2D matrix: num_MC X num_grps)
%           InputProductionVector       (t/km2/y); (horizontal vector: 1 X num_grps)
%           PhysicalLossFraction	(horizontal vector: 1 X num_grps)
%       ECOTRAN
%           GroupType
%           GroupTypeDef_           code number definitions for several groups
%           num_grps                number of functional groups
%           num_MC                  number of Monte Carlo models
%           num_nutrients
%           num_fleets
%           num_livingANDfleets
%           num_eggs
%           num_ANYdetritus
%
% returns:
%       StaticScenario_results
%           EnergyBudget_scenario           (proportions); (3D matrix: num_grps X num_grps X num_MC)
%           ConsumptionBudget_scenario      (proportions); (3D matrix: 7 X num_grps X num_MC); NOTE: FFF ConsumptionBudget is not changed for scenarios, this is for the future
%           DiscardFraction_scenario        (proportions); (3D matrix: num_grps X num_fleets X num_MC)
%           Production_scenario             (biomass/km2/time); (2D matrix: num_grps X num_MC)
%           fate_metabolism_scenario        (proportions); (3D matrix: num_nutrients X num_grps X num_MC)
%           fate_eggs_scenario              (proportions); (3D matrix: num_nutrients X num_grps X num_MC)
%           fate_feces_scenario             (proportions); (3D matrix: num_ANYdetritus X num_grps X num_MC)
%           fate_predation_scenario         (proportions); (3D matrix: num_livingANDfleets X num_grps X num_MC)
%           fate_senescence_scenario        (proportions); (3D matrix: num_ANYdetritus X num_grps X num_MC)
%           BuildGroup                      group code numbers; (2D matrix: num_grps X num_MC)
%           max_ScaleFactor                 maximum allowed scale factor for each group; (2D matrix: num_MC X num_grps)
%           LANDINGS_scenario               corrected landings of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
%           DISCARDS_scenario               corrected discards of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
%           CATCH_scenario                  catch of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
%           DiscardFraction_scenario        corrected discard fraction; (3D matrix: num_grps X num_fleets X num_MC)
%           catch_TotalByFleet_scenario     total catch for each fleet; (3D matrix: 1 X num_fleets X num_MC)
%           landings_TotalByFleet_scenario	total landings for each fleet; (3D matrix: 1 X num_fleets X num_MC)
%           discards_TotalByFleet_scenario	total discards for each fleet; (3D matrix: 1 X num_fleets X num_MC)
%           catch_TotalByGroup_scenario     total catch by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)
%           landings_TotalByGroup_scenario	total landings by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)
%           discards_TotalByGroup_scenario	total discards by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)
%           fname_ScenarioGenerator         name of f_ScenarioGenerator function
%
% revision date: 10-21-2020


% *************************************************************************
% STEP 1: operating parameters---------------------------------------------
fname_ScenarioGenerator     = mfilename; % name of this f_ScenarioGenerator function
display(['Running: ' fname_ScenarioGenerator])

minimum_residual            = 0.001; % minimum residual scaling factor for offset consumers (largest factor by which an offset consumer's take can be scaled down)
% *************************************************************************



% *************************************************************************
% STEP 2: unpack variables (preserve values within original structures)----
modify_consumer         = ScenarioConditions.modify_consumer;           % row number(s) of consumer group(s) to force-modify
ScaleFactor             = ScenarioConditions.ScaleFactor;               % factor by which to re-scale production to modify_consumer group(s); NOTE: value > 1 means increase flow to modify_consumer (and reduced flow to offset_consumer); < 1 means decrease flow to modify_consumer; NOTE: keep ScaleFactor >= 0 (negative value makes no sense)
offset_consumer         = ScenarioConditions.offset_consumer;           % row number(s) of consumer group(s) to modify as offset to force-modified grp(s)NOTE: change in modify_consumer grp(s); can be [] if offset is to be distributed among ALL consumer groups
trgt_producer           = ScenarioConditions.trgt_producer;             % column number(s) of producer group(s) to modify
% EnergyBudget_base       = ScenarioConditions.EnergyBudget_base;         % MonteCarlo set of EnergyBudget_base models; (3D matrix: num_grps X num_grps X num_MC); NOTE: no changes made directly to EnergyBudget_base, everything is done via the fates
ConsumptionBudget_base	= ScenarioConditions.ConsumptionBudget_base;	% MonteCarlo set of ConsumptionBudget_base models; (3D matrix: 7 X num_grps X num_MC)
fate_metabolism_base   	= ScenarioConditions.fate_metabolism_base;    	% fate of metabolism in base model; (3D matrix: num_nutrients X num_grps X num_MC)
fate_eggs_base          = ScenarioConditions.fate_eggs_base;           	% fate of eggs (reproduction) in base model; (3D matrix: num_eggs X num_grps X num_MC)
fate_feces_base      	= ScenarioConditions.fate_feces_base;          	% fate of feces detritus in base model; (3D matrix: num_ANYdetritus X num_grps X num_MC)
fate_senescence_base	= ScenarioConditions.fate_senescence_base;    	% fate of senescence detritus in base model; (3D matrix: num_ANYdetritus X num_grps X num_MC)
fate_predation_base     = ScenarioConditions.fate_predation_base;     	% fate of production among all predators in base model; (3D matrix: num_livingANDfleets X num_grps X num_MC)
DiscardFraction_base	= ScenarioConditions.DiscardFraction_base;      % Monte Carlo set of fractions of catch discarded for each group and fleet; (3D matrix: num_grps X num_fleets X num_MC)

TransferEfficiency      = ScenarioConditions.TransferEfficiency;        % (2D matrix: num_MC X num_grps)
InputProductionVector	= ScenarioConditions.InputProductionVector;     % (t/km2/y); (horizontal vector: 1 X num_grps)
PhysicalLossFraction	= ScenarioConditions.PhysicalLossFraction;      % (horizontal vector: 1 X num_grps)

if ScaleFactor < 0; error('Negative scalers make no sense. Use a value equal to or greater than ZERO.'); end

GroupType            	= ECOTRAN.GroupType;
num_grps              	= ECOTRAN.num_grps;             % number of functional groups
num_MC                	= ECOTRAN.num_MC;               % number of Monte Carlo models
num_nutrients           = ECOTRAN.num_nutrients;
num_fleets              = ECOTRAN.num_fleets;
num_livingANDfleets     = ECOTRAN.num_livingANDfleets;
num_eggs                = ECOTRAN.num_eggs;
num_ANYdetritus         = ECOTRAN.num_ANYdetritus;
% *************************************************************************



% *************************************************************************
% STEP 3: prep variables---------------------------------------------------
EnergyBudget_scenario           = zeros(num_grps, num_grps, num_MC);	% initialize; (3D matrix: num_grps X num_grps X num_MC)
ConsumptionBudget_scenario      = zeros(7, num_grps, num_MC);           % initialize; (3D matrix: 7 X num_grps X num_MC)
DiscardFraction_scenario        = zeros(num_grps, num_fleets, num_MC);  % initialize; (3D matrix: num_grps X num_fleets X num_MC);

Production_scenario             = zeros(num_grps, num_MC);              % initialize; (2D matrix: num_grps X num_MC)
max_ScaleFactor                 = zeros(num_MC, num_grps) * NaN;        % initialize; (2D matrix: num_MC X num_grps)
UseScaleFactor                  = zeros(1, num_grps);              	    % initialize; (horizontal vector: 1 X num_grps)
BuildGroup                      = zeros(num_grps, num_MC);              % initialize; (2D matrix: num_grps X num_MC)

fate_metabolism_scenario        = zeros(num_nutrients, num_grps, num_MC);       % initialize; (3D matrix: num_nutrients X num_grps X num_MC)
fate_eggs_scenario              = zeros(num_eggs, num_grps, num_MC);            % initialize; (3D matrix: num_nutrients X num_grps X num_MC)
fate_feces_scenario             = zeros(num_ANYdetritus, num_grps, num_MC);     % initialize; (3D matrix: num_ANYdetritus X num_grps X num_MC)
fate_predation_scenario         = zeros(num_livingANDfleets, num_grps, num_MC);	% initialize; (3D matrix: num_livingANDfleets X num_grps X num_MC)
fate_senescence_scenario        = zeros(num_ANYdetritus, num_grps, num_MC);     % initialize; (3D matrix: num_ANYdetritus X num_grps X num_MC)

LANDINGS_scenario               = zeros(num_grps, num_fleets, num_MC);	% initialize landings of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
DISCARDS_scenario               = zeros(num_grps, num_fleets, num_MC);	% initialize discards of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
CATCH_scenario                  = zeros(num_grps, num_fleets, num_MC);	% initialize of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
catch_TotalByFleet_scenario     = zeros(1, num_fleets, num_MC);         % initialize total catch for each fleet; (3D matrix: 1 X num_fleets X num_MC)
landings_TotalByFleet_scenario	= zeros(1, num_fleets, num_MC);         % initialize total landings for each fleet; (3D matrix: 1 X num_fleets X num_MC)
discards_TotalByFleet_scenario	= zeros(1, num_fleets, num_MC);         % initialize total discards for each fleet; (3D matrix: 1 X num_fleets X num_MC)
catch_TotalByGroup_scenario     = zeros(num_grps, 1, num_MC);           % initialize total catch by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)
landings_TotalByGroup_scenario	= zeros(num_grps, 1, num_MC);           % initialize total landings by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)
discards_TotalByGroup_scenario	= zeros(num_grps, 1, num_MC);           % initialize total discards by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)

% remove modify_consumer from offset_consumer list
looky_modify                    = find(ismember(offset_consumer, modify_consumer));
offset_consumer(looky_modify)	= [];
% *************************************************************************



% *************************************************************************
% STEP 4: find detritus, nutrient, & fleet groups--------------------------
looky_nutrients                 = find(floor(GroupType)	== ECOTRAN.GroupTypeDef_ANYNitroNutr);	% row addresses of nutrients
looky_ANYPrimaryProducer        = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYPrimaryProd);
looky_ANYconsumer               = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYConsumer);
% looky_ANYbacteria               = find(GroupType        == ECOTRAN.GroupTypeDef_plgcBacteria | GroupType == ECOTRAN.GroupTypeDef_bnthBacteria);
looky_ANYbacteria               = find(floor(GroupType) == ECOTRAN.GroupTypeDef_bacteria);
looky_fleets                    = find(floor(GroupType) == ECOTRAN.GroupTypeDef_fleet);
looky_livingANDfleets           = [looky_ANYPrimaryProducer; looky_ANYconsumer; looky_ANYbacteria; looky_fleets]; % includes primary producers & bacteria
looky_eggs                      = find(GroupType        == ECOTRAN.GroupTypeDef_eggs);
looky_terminalBNTHdetritus      = find(GroupType        == ECOTRAN.GroupTypeDef_terminalBnthDetr);
looky_ANYdetritus               = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYDetritus);

% consumer row addresses in fate_predation
looky_modify_consumer           = find(ismember(looky_livingANDfleets, modify_consumer));
looky_offset_consumer           = find(ismember(looky_livingANDfleets, offset_consumer));

num_trgt_producer               = length(trgt_producer);
num_modify_consumer             = length(looky_modify_consumer);
num_offset_consumer             = length(looky_offset_consumer);
% *************************************************************************



% *************************************************************************
% STEP 5: use "type model" ConsumptionBudget for fleets in all Monte Carlo models. 
%         NOTE: Fleet ConsumptionBudgets are revised below based on discard
%               rates of individual groups.
%         FFF: in future, make fleet correction based upon individual group discard rates at the f_E2E_MonteCarlo step
ConsumptionBudget_fleets                    = ConsumptionBudget_base(:, looky_fleets, 1); % (2D matrix: 7 X num_fleets)
ConsumptionBudget_base(:, looky_fleets, :)	= repmat(ConsumptionBudget_fleets, [1, 1, num_MC]); % (2D matrix: 7 X num_grps X num_MC)
% *************************************************************************



% *************************************************************************
% process each Monte Carlo (or "type" model)
for MonteCarlo_loop = 1:num_MC
    
    % *********************************************************************
    % STEP 6: get current model--------------------------------------------
    % current_EnergyBudget            = EnergyBudget_base(:, :, MonteCarlo_loop);         % (2D matrix: num_grps X num_grps); NOTE: no changes made directly to current_EnergyBudget, everything is done via the fates
    current_ConsumptionBudget       = ConsumptionBudget_base(:, :, MonteCarlo_loop);	% (2D matrix: 7 X num_grps)
    current_DiscardFraction_base	= DiscardFraction_base(:, :, MonteCarlo_loop);      % fractions of catch discarded for each group and fleet; (2D matrix: num_grps X num_fleets)
    current_fate_metabolism         = fate_metabolism_base(:, :, MonteCarlo_loop);      % (2D matrix: num_nutrients X num_grps)
    current_fate_eggs               = fate_eggs_base(:, :, MonteCarlo_loop);            % (2D matrix: num_eggs X num_grps)
    current_fate_feces              = fate_feces_base(:, :, MonteCarlo_loop);           % (2D matrix: num_ANYdetritus X num_grps)
    current_fate_senescence         = fate_senescence_base(:, :, MonteCarlo_loop);      % (2D matrix: num_ANYdetritus X num_grps)
    current_fate_predation          = fate_predation_base(:, :, MonteCarlo_loop);       % (2D matrix: num_livingANDfleets X num_grps)
    current_TransferEfficiency      = TransferEfficiency(MonteCarlo_loop, :);           % (2D matrix: 1 X num_grps)
    
    UseScaleFactor(1, 1:num_grps)	= ScaleFactor;                                  % initialize for each round of MonteCarlo_loop; (horizontal vector: 1 X num_grps)
    % *********************************************************************
    
    
    
    % *********************************************************************
    % STEP 7: changes to ConsumptionBudget---------------------------------
    %         NOTE: ba & em terms are NOT corrected for in the EnergyBudget (not even in base model) - that is, EnergyBudget sums do not equal 1 if there are non-zero ba or em terms
    %         FFF in future, allow for scaled changes to be directly made to the ConsumptionBudget
    modify_ConsumptionBudget                            = current_ConsumptionBudget;	% (2D matrix: 7 X num_grps)
    
    ConsumptionBudget_scenario(:, :, MonteCarlo_loop)	= modify_ConsumptionBudget;     % save modified ConsumptionBudget for each MonteCarlo_loop; NOTE: FFF ConsumptionBudget is not changed for scenarios, this is for the future
    % *********************************************************************
    
    
    
    % *********************************************************************
    % STEP 8: changes to DiscardFraction-----------------------------------
    %         FFF in future, apply forced changes to discard fractions here
    modify_DiscardFraction                              = current_DiscardFraction_base;	% fractions of catch discarded for each group and fleet; (2D matrix: num_grps X num_fleets)
    
   	DiscardFraction_scenario(:, :, MonteCarlo_loop)     = modify_DiscardFraction;	% (3D matrix: num_grps X num_fleets X num_MC)
    % *********************************************************************
    
    
    
    % *********************************************************************
    % STEP 9: changes to (non-predation) fates-----------------------------
    %         FFF in future, apply changes to (non-predation) fates here
    modify_fate_metabolism          = current_fate_metabolism;	% (2D matrix: num_nutrients X num_grps)
    modify_fate_eggs                = current_fate_eggs;        % (2D matrix: num_nutrients X num_grps)
    modify_fate_feces               = current_fate_feces;       % (2D matrix: num_ANYdetritus X num_grps)
    modify_fate_senescence          = current_fate_senescence;	% (2D matrix: num_ANYdetritus X num_grps)
    
   	fate_metabolism_scenario(:, :, MonteCarlo_loop)	= modify_fate_metabolism;	% (2D matrix: num_nutrients X num_grps)
    fate_eggs_scenario(:, :, MonteCarlo_loop)       = modify_fate_eggs;         % (2D matrix: num_nutrients X num_grps)
    fate_feces_scenario(:, :, MonteCarlo_loop)      = modify_fate_feces;        % (2D matrix: num_ANYdetritus X num_grps)
    fate_senescence_scenario(:, :, MonteCarlo_loop)	= modify_fate_senescence;	% (2D matrix: num_ANYdetritus X num_grps)
    % *********************************************************************
    
    
    
    % *********************************************************************
    % STEP 10: changes to predation fate-----------------------------------
    %         estimate maximum allowed scaling factor upon each producer group
    offset_consumer_matrix           	= current_fate_predation(looky_offset_consumer, :);	% predation fates to offset consumers; (2D matrix: num_offset_consumer X num_grps)
    min_offset_consumer_matrix          = offset_consumer_matrix .* minimum_residual;       % apply minimum trace on each element; (2D matrix: num_offset_consumer X num_grps)
    offset_consumer_freedom             = offset_consumer_matrix - min_offset_consumer_matrix; % maximum amount that can be taken from offset consumers; (2D matrix: num_offset_consumer X num_grps)
    sum_offset_consumer_freedom         = sum(offset_consumer_freedom, 1);                  % for each producer, the maximum amount that can be moved to the modify consumer; (horizontal vector: 1 X num_grps)
    
    modify_consumer_matrix              = current_fate_predation(looky_modify_consumer, :); % the original, un-scaled modify consumer take; (2D matrix: num_modify_consumer X num_grps)
    sum_modify_consumer_matrix          = sum(modify_consumer_matrix, 1);                   % for each producer, the total un-scaled take of all modify consumers; (horizontal vector: 1 X num_grps)
    
    max_scaler                                      = 1 + (sum_offset_consumer_freedom ./ sum_modify_consumer_matrix); % how many more times can the modify consumers be scaled up; (horizontal vector: 1 X num_grps)
    max_scaler(sum_modify_consumer_matrix == 0)     = ScaleFactor; % (horizontal vector: 1 X num_grps)
    max_ScaleFactor(MonteCarlo_loop, 1:num_grps)	= max_scaler; % build up a record; (2D matrix: num_MC X num_grps)

    if ScaleFactor > 1
        looky_ExceededLimits                    = find(max_scaler < ScaleFactor);
        if ~isempty(looky_ExceededLimits)
            display(['-->WARNING in MC-' num2str(MonteCarlo_loop) ': maximum scaler exceeded for target producer(s): ' num2str(looky_ExceededLimits)])
            UseScaleFactor(1, looky_ExceededLimits)	= max_scaler(looky_ExceededLimits); % replace values in UseScaleFactor where ScaleFactor exceedes max_scaler with max_scaler
        end % (~isempty(looky_ExceededLimits))
    end % (ScaleFactor > 1)
    
    % apply scaling factor to modify consumer(s) in fate_predation matrix
    modify_fate_predation                           = current_fate_predation; % (2D matrix: num_livingANDfleets X num_grps)
    modify_fate_predation(looky_modify_consumer, :)	= modify_fate_predation(looky_modify_consumer, :) .* repmat(UseScaleFactor, [num_modify_consumer, 1]);
    
    % change offset consumer(s); distribute in proportion to original values
    delta_matrix            = current_fate_predation - modify_fate_predation;	% (2D matrix: num_livingANDfleets X num_grps)
    sum_delta               = sum(delta_matrix, 1);                    % (horizontal vector: 1 X num_grps)
    
    if ~isempty(offset_consumer) % QQQ when would this ever be empty?
        
        offset_vectors          = current_fate_predation(looky_offset_consumer, :);	% (2D matrix: num_offset_consumer X num_grps)
        sum_offset              = sum(offset_vectors, 1);                           % (horizontal vector: 1 X num_grps)
        offset_fraction         = offset_vectors ./ repmat(sum_offset, [num_offset_consumer, 1]); % normalize each clm to 1; (2D matrix: num_offset_consumer X num_grps)
        offset_fraction(isnan(offset_fraction)) = 0; % remove NaN caused by div/0
        offset_fraction         = offset_fraction .* repmat(sum_delta, [num_offset_consumer, 1]); % scale offset_fraction by sum_delta; (2D matrix: num_offset_consumer X num_grps)
        offset_vectors          = offset_vectors + offset_fraction; % add offset_fraction to offset_vectors; (2D matrix: num_offset_consumer X num_grps)
        modify_fate_predation(looky_offset_consumer, :) = offset_vectors; % replace offset consumer values in modify_matrix; (2D matrix: num_offset_consumer X num_grps)

    end % (~isempty(offset_consumer))
    
	fate_predation_scenario(:, :, MonteCarlo_loop)	= modify_fate_predation;	% (3D matrix: looky_livingANDfleets X num_grps X num_MC)
    % *********************************************************************

    
    
    % *********************************************************************
    % STEP 11: build EnergyBudget_scenario---------------------------------
    %          replace terms only for the trgt_producer columns
        
    % paste in metabolism terms
    EnergyBudget_scenario(looky_nutrients, trgt_producer, MonteCarlo_loop) = ...
        modify_fate_metabolism(:, trgt_producer) .* repmat(modify_ConsumptionBudget(2, trgt_producer), [num_nutrients, 1]);
    
    % paste in egg terms
    EnergyBudget_scenario(looky_eggs, trgt_producer, MonteCarlo_loop) = ...
        modify_fate_eggs(:, trgt_producer) .* repmat(modify_ConsumptionBudget(3, trgt_producer), [num_eggs, 1]);

    % paste in predation terms
    EnergyBudget_scenario(looky_livingANDfleets, trgt_producer, MonteCarlo_loop) = ...
        modify_fate_predation(:, trgt_producer) .* repmat(modify_ConsumptionBudget(4, trgt_producer), [num_livingANDfleets, 1]);
    
    % NOTE: special case for terminal benthic detritus group(s) where predator rows are 
    %       scaled by predation + senescence so that their column sum(s) add to 1
    if ismember(looky_terminalBNTHdetritus, trgt_producer) % test if terminalBNTHdetritus is a trgt_producer
        EnergyBudget_scenario(looky_livingANDfleets, looky_terminalBNTHdetritus, MonteCarlo_loop) = ...
            EnergyBudget_scenario(looky_livingANDfleets, looky_terminalBNTHdetritus, MonteCarlo_loop) + ...
            modify_fate_predation(:, looky_terminalBNTHdetritus) .* repmat(modify_ConsumptionBudget(5, looky_terminalBNTHdetritus), [num_livingANDfleets, 1]);
    end
    
    % paste in feces terms
    EnergyBudget_scenario(looky_ANYdetritus, trgt_producer, MonteCarlo_loop) = ...
        modify_fate_feces(:, trgt_producer) .* repmat(modify_ConsumptionBudget(1, trgt_producer), [num_ANYdetritus, 1]);
    
    % ADD in senescence terms (add to feces terms just pasted into detritus)
    EnergyBudget_scenario(looky_ANYdetritus, trgt_producer, MonteCarlo_loop) = ...
        EnergyBudget_scenario(looky_ANYdetritus, trgt_producer, MonteCarlo_loop) + ...
        (modify_fate_senescence(:, trgt_producer) .* repmat(modify_ConsumptionBudget(5, trgt_producer), [num_ANYdetritus, 1]));
    % *********************************************************************

    
    
    % *********************************************************************
    % STEP 12: distinguish landings & discards in fleet production---------
    %          Adjusting fleet columns in EnergyBudget & ConsumptionBudget
	%
    %           NOTE: Recalculating fleet em (landings) & feces (discard) terms based
    %           on discard rates of individual functional groups at levels of production
    %           calculated via f_WebProductivity function. This will lead
    %           to new rates of detritus (offal) production and new rates
    %           of production by all groups eating offal -- requiring
    %           revised calculations of landings and discard rates of
    %           fleets, new calculations of offal production, and repeat...
    %           Rather than recalculate revised total fleet discard rates
    %           interatively (and I have no certainty that this will
    %           converge to an approximately terminal value, though I
    %           suspect it should), I will make this adjustment just once.
    %           I am using the "type model" fleet EnergyBudget & ConsumptionBudget
    %           terms for all Monte Carlo models, thus fleet total discard
    %           rates are tied to Monte Carlo discard rates of individual
    %           groups rather than the previously generated Monte Carlo estimate 
    %           of total fleet discard.
    %
    %           NOTE: The total fleet discard rate will not exactly match
    %           that implied by the discard rates of individual groups (as
    %           described above), but because offal discard is generally a
    %           small part of detritivore diets this error should be
    %           trivial in most cases. This will NOT cause a conflict among
    %           the subsequent calculations of productivity of total fleet
    %           offal as all fleet discard contributions back into the
    %           environment are based upon the fleet EnergyBudget and
    %           ConsumptionBudget terms and not the discard rates of
    %           individual groups.

    % step 12a: get current Monte Carlo model -----------------------------
    EnergyBudget_current	= EnergyBudget_scenario(:, :, MonteCarlo_loop); % QQQ call this modify_EnergyBudget?? (2D matrix: num_grps X num_grps)

    % step 12b: calculate current group production rates ------------------
    production_current      = f_WebProductivity_03272019(current_TransferEfficiency, EnergyBudget_current, InputProductionVector, PhysicalLossFraction); % production (actually CONSUMPTION) rates; (mmole N/m3/d); (vertical vector: num_grps X 1)
    production_current      = production_current'; % transpose; (mmole N/m3/d); (horizontal vector: 1 X num_grps)
    EnergyBudget_fleet      = EnergyBudget_current(looky_fleets, :); % fraction of flow to each fleet; (2D matrix: num_fleets X num_grps)
    production_repmat       = repmat(production_current, [num_fleets, 1]); % group production rates; (t/km2/y); (2D matrix: num_fleets X num_grps)
    
    % step 12c: calculate current catch & discard rates for fleets & individual groups
    CATCH_current                   = EnergyBudget_fleet .* production_repmat; % catch rates for each group; (t/km2/y); (2D matrix: num_fleets X num_grps)
    DISCARDS_current                = CATCH_current .* modify_DiscardFraction'; % discard rates for each group; (t/km2/y); (2D matrix: num_fleets X num_grps); NOTE transpose
	LANDINGS_current                = CATCH_current .* (1 - modify_DiscardFraction'); % landing rates for each group; (t/km2/y); (2D matrix: num_fleets X num_grps); NOTE transpose

    catch_TotalByFleet_current      = sum(CATCH_current, 2);	% total catch for each fleet in current MC model; (vertical vector: num_fleets X 1)
	landings_TotalByFleet_current	= sum(LANDINGS_current, 2);	% total landings for each fleet in current MC model; (vertical vector: num_fleets X 1)
	discards_TotalByFleet_current	= sum(DISCARDS_current, 2);	% total discards for each fleet in current MC model; (vertical vector: num_fleets X 1)

    catch_TotalByGroup_current      = sum(CATCH_current);       % total catch by group (summed across fleets); (horizontal vector: 1 X num_grps)
	landings_TotalByGroup_current	= sum(LANDINGS_current);	% total landings by group (summed across fleets); (horizontal vector: 1 X num_grps)
	discards_TotalByGroup_current	= sum(DISCARDS_current);	% total discards by group (summed across fleets); (horizontal vector: 1 X num_grps)

    % step 12d: revise EnergyBudget & ConsumptionBudget -------------------
    landings_fraction = landings_TotalByFleet_current ./ catch_TotalByFleet_current; % (vertical vector: num_fleets X 1)
    landings_fraction(isnan(landings_fraction)) = 0; % catch any div/0 error (only needed if a fleet has 0 catch)
    
    modify_ConsumptionBudget(7, looky_fleets) = landings_fraction'; % update em term (landings removed); (NOTE transpose)
    modify_ConsumptionBudget(1, looky_fleets) = (1 - landings_fraction'); % update feces term (discards); (NOTE transpose)
    
    discards_fraction = (1 - landings_fraction'); % (horizontal vector: 1 X num_fleets); (NOTE transpose)
    EnergyBudget_current(looky_ANYdetritus, looky_fleets) = modify_fate_feces(:, looky_fleets) .* repmat(discards_fraction, [num_ANYdetritus, 1]);
    % *********************************************************************

    
    
	% *********************************************************************
    % STEP 13: calculate production rates using modified, scenario matrix--
    Production_scenario(:, MonteCarlo_loop) = f_WebProductivity_03272019(current_TransferEfficiency, EnergyBudget_scenario(:, :, MonteCarlo_loop), InputProductionVector, PhysicalLossFraction); % production (actually CONSUMPTION) rates; (mmole N/m3/d); (2D matrix: num_grps X num_MC)
    BuildGroup(1:num_grps, MonteCarlo_loop) = (1:num_grps)';
    % *********************************************************************
    

    
	% *********************************************************************
    % step 14: build-up scenario results ----------------------------------
    EnergyBudget_scenario(:, :, MonteCarlo_loop)            = EnergyBudget_current; % (3D matrix: num_grps X num_grps X num_MC)
    ConsumptionBudget_scenario(:, :, MonteCarlo_loop)       = modify_ConsumptionBudget; % (3D matrix: 7 X num_grps X num_MC)
    
    LANDINGS_scenario(:, :, MonteCarlo_loop)                = LANDINGS_current';	% landings of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC); NOTE transpose
    DISCARDS_scenario(:, :, MonteCarlo_loop)                = DISCARDS_current';	% discards of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC); NOTE transpose
    CATCH_scenario(:, :, MonteCarlo_loop)                   = CATCH_current';       % catch of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC); NOTE transpose
    
    catch_TotalByFleet_scenario(1, :, MonteCarlo_loop)      = catch_TotalByFleet_current';      % total catch for each fleet in current MC model; (3D matrix: 1 X num_fleets X num_MC); NOTE transpose
    landings_TotalByFleet_scenario(1, :, MonteCarlo_loop)	= landings_TotalByFleet_current';	% total landings for each fleet in current MC model; (3D matrix: 1 X num_fleets X num_MC); NOTE transpose
    discards_TotalByFleet_scenario(1, :, MonteCarlo_loop)	= discards_TotalByFleet_current';	% total discards for each fleet in current MC model; (3D matrix: 1 X num_fleets X num_MC); NOTE transpose
    
    catch_TotalByGroup_scenario(:, 1, MonteCarlo_loop)      = catch_TotalByGroup_current';     % total catch by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC); NOTE transpose
    landings_TotalByGroup_scenario(:, 1, MonteCarlo_loop)	= landings_TotalByGroup_current';	% total landings by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC); NOTE transpose
    discards_TotalByGroup_scenario(:, 1, MonteCarlo_loop)	= discards_TotalByGroup_current';	% total discards by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC); NOTE transpose
    % *********************************************************************
    
end % (MonteCarlo_loop)



% *************************************************************************
% STEP 15: compile results-------------------------------------------------
StaticScenario_results.EnergyBudget_scenario            = EnergyBudget_scenario;            % (proportions); (3D matrix: num_grps X num_grps X num_MC)
StaticScenario_results.ConsumptionBudget_scenario       = ConsumptionBudget_scenario;       % (proportions); (3D matrix: 7 X num_grps X num_MC); NOTE: FFF ConsumptionBudget is not changed for scenarios, this is for the future

StaticScenario_results.Production_scenario              = Production_scenario;              % (biomass/km2/time); (2D matrix: num_grps X num_MC)

StaticScenario_results.fate_metabolism_scenario         = fate_metabolism_scenario;         % (proportions); (3D matrix: num_nutrients X num_grps X num_MC)
StaticScenario_results.fate_eggs_scenario               = fate_eggs_scenario;               % (proportions); (3D matrix: num_nutrients X num_grps X num_MC)
StaticScenario_results.fate_feces_scenario              = fate_feces_scenario;              % (proportions); (3D matrix: num_ANYdetritus X num_grps X num_MC)
StaticScenario_results.fate_predation_scenario          = fate_predation_scenario;          % (proportions); (3D matrix: num_livingANDfleets X num_grps X num_MC)
StaticScenario_results.fate_senescence_scenario         = fate_senescence_scenario;         % (proportions); (3D matrix: num_ANYdetritus X num_grps X num_MC)

StaticScenario_results.BuildGroup                       = BuildGroup;                       % group code numbers; (2D matrix: num_grps X num_MC)
StaticScenario_results.max_ScaleFactor                  = max_ScaleFactor;                  % maximum allowed scale factor for each group; (2D matrix: num_MC X num_grps)

StaticScenario_results.LANDINGS_scenario                = LANDINGS_scenario;                % corrected landings of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
StaticScenario_results.DISCARDS_scenario                = DISCARDS_scenario;                % corrected discards of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
StaticScenario_results.CATCH_scenario                   = CATCH_scenario;                   % catch of each group by each fleet; (mmole N/m3/d); (3D matrix: num_grps X num_fleets X num_MC)
StaticScenario_results.DiscardFraction_scenario         = DiscardFraction_scenario;         % discard fraction; (proportions); (3D matrix: num_grps X num_fleets X num_MC)
StaticScenario_results.catch_TotalByFleet_scenario      = catch_TotalByFleet_scenario;      % total catch for each fleet; (3D matrix: 1 X num_fleets X num_MC)
StaticScenario_results.landings_TotalByFleet_scenario	= landings_TotalByFleet_scenario;	% total landings for each fleet; (3D matrix: 1 X num_fleets X num_MC)
StaticScenario_results.discards_TotalByFleet_scenario	= discards_TotalByFleet_scenario;	% total discards for each fleet; (3D matrix: 1 X num_fleets X num_MC)
StaticScenario_results.catch_TotalByGroup_scenario      = catch_TotalByGroup_scenario;      % total catch by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)
StaticScenario_results.landings_TotalByGroup_scenario	= landings_TotalByGroup_scenario;	% total landings by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)
StaticScenario_results.discards_TotalByGroup_scenario	= discards_TotalByGroup_scenario;	% total discards by group (summed across fleets); (3D matrix: num_grps X 1 X num_MC)

StaticScenario_results.fname_ScenarioGenerator          = fname_ScenarioGenerator;          % name of f_ScenarioGenerator function
% *************************************************************************


% end m-file***************************************************************