function [REACH_vector, FOOTPRINT_trace, REACH_trace] = ...
         f_DietTrace_08032020(ECOTRAN, DietArray, TraceGroup)
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
% NOTE: The code uses the diet matrix to define the food web.
% NOTE: This function analyzes one (1) Monte Carlo static food web.
%
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%       ECOTRAN
%           num_grps
%           GroupType
%           GroupTypeDef_           (several)
%       DietArray                   (2D matrix: num_grps (producers) X num_grps (consumers))
%       TraceGroup               	(address of the model group to be traced)
%
% returns:
%       REACH_vector     - "REACH" of TraceGroup = producer (rows) to each consumer (columns)
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
%
% revision date: 8-3-2020


% *************************************************************************
% STEP 1: set operating parameters-----------------------------------------
iteration_tolerance         = 0.0001; % iterate until working_trace stabilizes
iteration_max               = 1000; % maximum number of iterations
% *************************************************************************



% *************************************************************************
% STEP 2: unpack ECOTRAN info----------------------------------------------
% step 2a: ----------------------------------------------------------------
num_grps                    = ECOTRAN.num_grps + 1; % NOTE: add 1 to account for import diet
GroupType                   = ECOTRAN.GroupType;
% -------------------------------------------------------------------------


% step 2b: find group addresses -------------------------------------------
looky_terminalBnthDetr      = find(GroupType == ECOTRAN.GroupTypeDef_terminalBnthDetr); % address of the single terminal benthic detritus group
looky_ANYconsumer           = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYConsumer); 
looky_fleets                = find(floor(GroupType) == ECOTRAN.GroupTypeDef_fleet);
looky_TraceGroup            = TraceGroup;
looky_import                = num_grps; % import is the last row and column
looky_ANYconsumer           = [looky_ANYconsumer; looky_fleets]; % append fleets to consumers; NOTE: this & original code do not include bacteria as consumers
looky_NonConsumer           = 1:num_grps;
looky_NonConsumer(looky_ANYconsumer) = []; % addresses of nutrients, eggs, detritus, bacteria, & import (= fleet landings)
% *************************************************************************



% *************************************************************************
% STEP 3: prep diet matrix-------------------------------------------------
working_Diet                        = DietArray; % (2D matrix: num_grps (producers) X num_grps (consumers))
working_Diet(:, looky_NonConsumer)	= 0; % set diet columns of non-predators to 0; (2D matrix: num_grps (producers) X num_grps (consumers))
% *************************************************************************



% *************************************************************************
% STEP 4: perform top-down diet trace--------------------------------------
% step 4a: initialize variables -------------------------------------------
working_trace                       = working_Diet; % QQQ initialize as zeros?; (2D matrix: num_grps (producers) X num_grps (consumers))
trace_fraction_previous             = ones(num_grps, 1); % initialize trace_fraction_previous; (vertical vector: num_grps X 1)
check_iteration                     = 1; % initialize check_iteration
iteration_count                     = 0; % initialize iteration_count    
% -------------------------------------------------------------------------


% step 4b: prepare TraceVector using the TraceGroup -----------------------
TraceVector                         = working_Diet(:, looky_TraceGroup)'; % TraceGroup-as-predator; amount of each producer in diet of TraceGroup (across multiple trophic steps); (horizontal vector: 1 X num_grps); NOTE transpose
TraceVector(1, looky_TraceGroup)	= 1; % contribution of TraceGroup to TraceGroup itself = 1; (horizontal vector: 1 X num_grps)
% -------------------------------------------------------------------------


% step 4c: multiply predation on TraceGroup through diet matrix, iteratively
%          NOTE: TraceGroup (and TraceVector) is now acting as a predator (consumer)
while (check_iteration > iteration_tolerance) && (iteration_count < iteration_max)

    % multiply TraceVector through the diet matrix row-wise
	for prey_loop = 1:num_grps
        working_trace(prey_loop, :) = TraceVector .* working_Diet(prey_loop, :); % (2D matrix: num_grps (producers) X num_grps (consumers))
	end % (prey_loop)
        
	trace_fraction_current          = sum(working_trace, 2); % fraction of trace group in diet of each predator; (vertical vector: num_grps X 1)
        
    % make new TraceVector
	TraceVector                     = trace_fraction_current'; % (horizontal vector: 1 X num_grps); NOTE transpose
    TraceVector(looky_TraceGroup)	= 1; % contribution of TraceGroup to TraceGroup itself = 1; (horizontal vector: 1 X num_grps)
    
    % check status of iterations
	check_iteration                 = max(abs(trace_fraction_previous - trace_fraction_current));
	iteration_count                 = iteration_count + 1;
	trace_fraction_previous         = trace_fraction_current;

end % while (check_iteration > iteration_tolerance) && (iteration_count < iteration_max)
% -------------------------------------------------------------------------


% step 4d: pack up results ------------------------------------------------
%          NOTE: TraceFraction_downward is identical (within several decimal places)
%                to REACH_vector (a.k.a, TraceFraction_upward) and is not used beyond this point.
FOOTPRINT_trace             = working_trace; % (2D matrix: num_grps (producers) X num_grps (consumers))
TraceFraction_downward      = trace_fraction_current; % (vertical vector: num_grps X 1)
TraceFraction_downward(looky_TraceGroup, 1) = 1; % set origin group to 1; contribution of TraceGroup to TraceGroup itself = 1; (vertical vector: num_grps X 1)
% *************************************************************************



% *************************************************************************
% STEP 5: perform bottom-up diet trace ------------------------------------
% step 5a: initialize variables -------------------------------------------
working_trace               = working_Diet; % QQQ initialize as zeros?; (2D matrix: num_grps (producers) X num_grps (consumers)
trace_fraction_previous     = ones(1, num_grps); % initialize; (horizontal vector: 1 X num_grps)
check_iteration             = 1; % initialize check_iteration
iteration_count             = 0; % initialize iteration_count
% -------------------------------------------------------------------------


% step 5b: prepare TraceVector using the TraceGroup -----------------------
%          NOTE: TraceGroup and TraceVector is now acting as prey (producer)
TraceVector                     = working_Diet(looky_TraceGroup, :)'; % TraceGroup-as-prey; amount of TraceGroup in diet of each predator (across multiple trophic steps); (vertical vector: num_grps X 1); NOTE transpose
TraceVector(looky_TraceGroup)	= 1; % contribution of TraceGroup to TraceGroup itself = 1; (vertical vector: num_grps X 1)
% -------------------------------------------------------------------------


% step 5c: multiply predation on TraceGroup through diet, iteratively -----
while (check_iteration > iteration_tolerance) && (iteration_count < iteration_max)
       
    % multiply TraceVector through the diet matrix column-wise
	for predator_loop = 1:num_grps
        working_trace(:, predator_loop) = TraceVector .* working_Diet(:, predator_loop); % (2D matrix: num_grps (producers) X num_grps (consumers)
    end
    
	trace_fraction_current          = sum(working_trace); % fraction of trace group in diet of each predator; (horizontal vector: 1 X num_grps)

	% make new TraceVector
	TraceVector                     = trace_fraction_current'; % (vertical vector: num_grps X 1); NOTE transpose
	TraceVector(looky_TraceGroup)	= 1; % contribution of TraceGroup to TraceGroup itself = 1; (vertical vector: num_grps X 1)

	% check status of iterations
    check_iteration         = max(abs(trace_fraction_previous - trace_fraction_current));
	iteration_count         = iteration_count + 1;
	trace_fraction_previous	= trace_fraction_current;

end % while (check_iteration > iteration_tolerance) && (iteration_count < iteration_max)
% -------------------------------------------------------------------------


% step 5d: pack up results ------------------------------------------------
REACH_trace                         = working_trace ./ working_Diet; % divide trace by diet to get true fraction of trace_group contribution in each linkage path; (2D matrix: num_grps (producers) X num_grps (consumers)
REACH_trace(isnan(REACH_trace))     = 0; % set div/0 errors to 0; (2D matrix: num_grps (producers) X num_grps (consumers)

REACH_vector                        = trace_fraction_current; % (horizontal vector: 1 X num_grps)
REACH_vector(1, looky_TraceGroup)	= 1; % set TraceGroup  to 1; contribution of TraceGroup to TraceGroup itself = 1; (horizontal vector: 1 X num_grps)
% *************************************************************************


% end m-file***************************************************************