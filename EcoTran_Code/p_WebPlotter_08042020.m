function p_WebPlotter_08042020(ECOTRAN, FOOTPRINT_array, FOOTPRINT_trace, REACH_array, REACH_trace, x_position, TL_position, GRPlabels, TraceGroup, grid_switch)
% plot a food web diagram with footprint and reach information
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%       ECOTRAN
%           biomass                                                                                             (vertical vector: num_grps X 1)
%           CONSUMPTION         consumption rate matrix, who eats how much of what                              (2D matrix: num_grps (producers) X num_grps (consumers))
%           consumption_total	domestic + import consumption rate of each group; use in place of production    (vertical vector: num_grps X 1)
%           consumption_import	import consumption rate of each group                                           (vertical vector: num_grps X 1)
%           num_grps
%           GroupType           type of functional group
%           GroupTypeDef_       several group type definition codes
%       FOOTPRINT_array         fraction of each PRODUCER's production flowing to TraceGroup = CONSUMER; use for web plotting footprint box colors                  (3D matrix: num_grps (consumers) X num_grps (producers) X 2)
%       FOOTPRINT_trace         fractional contribution of each linkage to production of CONSUMER = TraceGroup; use for web plotting footprint arrow colors         (3D matrix: num_grps+1 (producers, +import) X num_grps+1 (consumers, +import) X 2)
%       REACH_array             fraction of each CONSUMER production originating from PRODUCER = TraceGroup; use for web plotting reach box colors                  (3D matrix: num_grps+1 (producers; +import) X num_grps+1 (consumers; +import) X 2)
%       REACH_trace             fraction of energy within each linkage ultimately originating from PRODUCER = TraceGroup; use for web plotting reach arrow colors	(3D matrix: num_grps+1 (producers, +import) X num_grps+1 (consumers, +import) X 2)
%       x_position              x-position of each box
%       TL_position             y-position of each box
%       GRPlabels               labels for boxes
%       TraceGroup              the group (by position address in model) for which to plot footprint & reach metrics
%       grid_switch             'on' to turn on gridlines for better placing box x_position & TL_position
%
% returns:
%       1 horrendogram
%
% revision date: 8-4-2020


% *************************************************************************
% STEP 1: set color intensity; line width; minimum flow size to plot-------
BoxIntensity  = 0.4; % intensity of box trace color (smaller than 1 intensifies)
LineIntensity = 0.2; % intensity of flow lines (smaller than 1 intensifies)

minline       = 1;
maxline       = 15;
% minspot       = 7; % QQQ
% maxspot       = 20; % QQQ

minspot       = 10; % QQQ
maxspot       = 22; % QQQ

minflow       = 0.0000005; % log(0.01 + 1); % minimum flow size to plot (0.025)
grayness      = 0.9; % 1 is white, 0 is black

xmin          = 0;
xmax          = max(x_position);
ymin          = 0.5; % QQQ
% ymin          = 0.75;
ymax          = 5; % QQQ
% ymax          = 5.3;
% *************************************************************************



% *************************************************************************
% STEP 2: prepare variables------------------------------------------------
% step 2a: get variables from ECOTRAN structure ---------------------------
biomass             = ECOTRAN.biomass;              % (vertical vector: num_grps X 1)
CONSUMPTION         = ECOTRAN.CONSUMPTION;          % (2D matrix: num_grps (producers) X num_grps (consumers))
consumption_total	= ECOTRAN.consumption_total;    % domestic + import consumption rate of each group; use in place of production; (vertical vector: num_grps X 1)
consumption_import  = ECOTRAN.consumption_import;   % import consumption rate of each group; (vertical vector: num_grps X 1)
GroupType           = ECOTRAN.GroupType;
num_grps            = ECOTRAN.num_grps;
footprint_vector	= FOOTPRINT_array(TraceGroup, :, 1); % (horizontal vector: 1 X num_grps)
reach_vector        = REACH_array(TraceGroup, :, 1); % (horizontal vector: 1 X num_grps+1 (+import)); NOTE: import column can be deleted
% -------------------------------------------------------------------------


% step 2b: get groups addresses -------------------------------------------
looky_nutrients             = find(floor(GroupType)	== ECOTRAN.GroupTypeDef_ANYNitroNutr);	% row addresses of nutrients
looky_ANYPrimaryProducer	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYPrimaryProd);
looky_eggs               	= find(GroupType        == ECOTRAN.GroupTypeDef_eggs);
looky_ANYdetritus        	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYDetritus);
looky_fleets                = find(floor(GroupType) == ECOTRAN.GroupTypeDef_fleet);
looky_NONconsumers          = [looky_nutrients; looky_ANYPrimaryProducer; looky_eggs; looky_ANYdetritus];
looky_import                = num_grps + 1;
% *************************************************************************



% *************************************************************************
% STEP 3: append import prey info as needed--------------------------------
CONSUMPTION(looky_import, :)	= consumption_import'; % (2D matrix: num_grps+1 (producers; + import) X num_grps (consumers))
consumption_total(looky_import)	= sum(consumption_import); % (vertical vector: num_grps+1 (+import) X 1)
footprint_vector(looky_import)  = 0; % add footprint on import = 0; (horizontal vector: 1 X num_grps+1 (+null import))
% *************************************************************************



% *************************************************************************
% STEP 4: ZERO-out non-consumers in CONSUMPTION columns--------------------
CONSUMPTION(:, looky_NONconsumers)    = 0;
% *************************************************************************



% *************************************************************************
% STEP 5: PLOT WEB!--------------------------------------------------------
% step 5a: set up and initialize figure------------------------------------
num_producers   = num_grps + 1; % number of producers; +1 more for import diet)
num_consumers   = num_grps; % number of consumers
linerange       = log((max(max(CONSUMPTION)))+1);

fig1            = figure; hold on;
axis([xmin xmax ymin ymax])
% -------------------------------------------------------------------------


% step 5b: plot gray link scafolding --------------------------------------
for producer_loop = 1:num_producers
    for consumer_loop = 1:num_consumers
        current_line        = CONSUMPTION(producer_loop, consumer_loop);
        linesize            = minline + ((maxline - minline)/(linerange))*log(current_line + 1);
        current_prey_x      = x_position(producer_loop);
        current_prey_y      = TL_position(producer_loop);
        current_predator_x  = x_position(consumer_loop);
        current_predator_y  = TL_position(consumer_loop);
        flow_color          = [0.6 0.6 0.6];
        % plot lines
        if ~isnan(current_line) && current_line >= minflow
            plot([current_prey_x current_predator_x], [current_prey_y current_predator_y], 'color', flow_color, 'linewidth', linesize);
        end
    end
end
% -------------------------------------------------------------------------


% step 5c: plot FOOTPRINT -------------------------------------------------
%          takes:
%               footprint_vector
%               FOOTPRINT_trace

[footprint_intensity, clm_order] = sort(footprint_vector(1:(end-1))); % ignoring the appended null import value at end
if isempty(clm_order); clm_order = 1:num_consumers; end

% plot links in reverse order of intensity so that gray lines do not hide colored lines        
for consumer_loop = clm_order
	for producer_loop = 1:num_producers
        
        current_line       = CONSUMPTION(producer_loop, consumer_loop);
        linesize           = minline + ((maxline - minline)/(linerange))*log(current_line + 1);
        current_prey_x     = x_position(producer_loop);
        current_prey_y     = TL_position(producer_loop);
        current_predator_x = x_position(consumer_loop);
        current_predator_y = TL_position(consumer_loop);

        % set energy flow line color
        if ~isempty(FOOTPRINT_trace(:, :, 1))
        	current_trace = FOOTPRINT_trace(producer_loop, consumer_loop, 1).^LineIntensity; % FOOTPRINT_trace fraction with increased color intensity
            % current_trace = footprint_vector(producer_loop).^LineIntensity; % current_trace with increased color intensity
            
            if current_trace > 0
                flow_color = [(1-current_trace) 1 (1-current_trace)]; % color becomes more white as current_trace decreases
                flow_color = flow_color * grayness;
                % plot lines
                if ~isnan(current_line) && current_line >= minflow
                    plot([current_prey_x current_predator_x], [current_prey_y current_predator_y], 'color', flow_color, 'linewidth', linesize);
                end % (~isnan(current_line) && current_line >= minflow)
            end % (current_trace > 0)
        end % (~isempty(FOOTPRINT_trace(:, :, 1)))
    end % (producer_loop)
end % (consumer_loop)
% -------------------------------------------------------------------------


% step 5d: plot REACH -----------------------------------------------------
%          takes:
%               TraceFraction_upward
%               DietTrace_upward

[reach_intensity, row_order] = sort(reach_vector);
if isempty(row_order); row_order = 1:num_producers; end

% plot links in reverse order of intensity so that gray lines do not hide colored lines
for producer_loop = row_order
	for consumer_loop = 1:num_consumers

        current_line       = CONSUMPTION(producer_loop, consumer_loop);
        linesize           = minline + ((maxline - minline)/(linerange))*log(current_line + 1);
        current_prey_x     = x_position(producer_loop);
        current_prey_y     = TL_position(producer_loop);
        current_predator_x = x_position(consumer_loop);
        current_predator_y = TL_position(consumer_loop);

        % set energy flow line color
        if ~isempty(REACH_trace(:, :, 1))
        	current_trace = REACH_trace(producer_loop, consumer_loop, 1).^LineIntensity; % REACH_trace fraction with increased color intensity
            % current_trace = reach_vector(producer_loop).^LineIntensity; % current_trace with increased color intensity
            
            if current_trace > 0
                flow_color = [1 (1-current_trace) (1-current_trace)]; % color becomes more white as current_trace decreases
                flow_color = flow_color * grayness;
                % plot lines
                if ~isnan(current_line) && current_line >= minflow
                    plot([current_prey_x current_predator_x], [current_prey_y current_predator_y], 'color', flow_color, 'linewidth', linesize);
                end
            end
        end
     end
end
% -------------------------------------------------------------------------


% step 5e: plot group boxes -----------------------------------------------
% spotrange       = max(biomass);
spotrange       = max(consumption_total);
edgecolor       = [0.5 0.5 0.5];

if isempty(reach_vector)
	reach_vector    = zeros(1, num_consumers); % NOTE: final value is import "consumer" and can be deleted
end % (isempty(reach_vector))
    
if isempty(footprint_vector)
    footprint_vector         = zeros(1, num_producers);
end % (isempty(footprint_vector))

for producer_loop = 1:num_producers
    
	% current_spot        = biomass(producer_loop);
 	current_spot        = consumption_total(producer_loop);
	spotsize            = minspot + ((maxspot - minspot)/(spotrange))*abs(current_spot);
    
    % set box color (plot either reach or footprint color for a given box based on which is the greater value)
	current_reach       = reach_vector(producer_loop);
	current_footprint   = footprint_vector(producer_loop);
    
    if current_footprint > current_reach
        current_BoxTrace    = current_footprint.^BoxIntensity;                  % current_trace with increased color intensity
        spotcolor           = [(1-current_BoxTrace) 1 (1-current_BoxTrace)];    % color becomes more white as current_trace
    else
        current_BoxTrace    = current_reach.^BoxIntensity;                      % current_trace with increased color intensity
        spotcolor           = [1 (1-current_BoxTrace) (1-current_BoxTrace)];    % color becomes more white as current_trace
    end % (current_footprint > current_reach)

    if ~isnan(x_position(producer_loop)) && ~isnan(spotsize)
        % plot(producer_loop, dat.TL(producer_loop), 'o', 'markerfacecolor', spotcolor, 'markeredgecolor', edgecolor, 'markersize', spotsize)
        text(x_position(producer_loop), TL_position(producer_loop), GRPlabels{producer_loop}, 'horizontalalignment', 'center', 'fontsize', spotsize, 'BackgroundColor', spotcolor, 'margin', 2, 'edgecolor', edgecolor)
    end % (~isnan(x_position(producer_loop)) && ~isnan(spotsize))

end % (producer_loop)
% -------------------------------------------------------------------------


% step 5f: turn on for grid -----------------------------------------------
if strcmp(grid_switch, 'on')
    set(gca, 'xtick', [xmin:0.5:xmax], 'ytick', [ymin:0.5:ymax])
    grid on
    grid minor
% else
%     axis off
end
% -------------------------------------------------------------------------


% step 6g: set axis range and labels --------------------------------------
font_size = 12;
set(gca, 'ytick', [ymin 1:0.5:ymax])
set(gca, 'xticklabel', [], 'ticklength', [0 0], 'fontweight', 'bold', 'fontsize', font_size);
set(gca, 'yticklabel', [{''}, {'1.0'}, {'1.5'}, {'2.0'}, {'2.5'}, {'3.0'}, {'3.5'}, {'4.0'}, {'4.5'}, {'5.0'}, {'5.5'}]);
% set(gca, 'yticklabel', []);

ylabel('Trophic Level', 'fontweight', 'bold', 'fontsize', font_size)
% axis off
% -------------------------------------------------------------------------


% end m-file***************************************************************