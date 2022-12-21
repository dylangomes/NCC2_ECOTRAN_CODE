function [u_final, v_final] = f_read_NWPO3data_02012022(WindParams)
% this function calculates median winds over desired time period
% and using desired time resolution
% time resolution must be in days
% takes:
%         WindParams.dt
%         WindParams.t_grid
%         WindParams.datestart
%         WindParams.dateend
%         WindParams.file_winddata
% calls:
%         none
% revision date: 2-1-2022


% *************************************************************************
% STEP 1: unpack variables-------------------------------------------------
dt              = WindParams.dt;
t_grid          = WindParams.t_grid;
datestart       = WindParams.datestart;
dateend         = WindParams.dateend;
file_winddata	= WindParams.file_winddata;
% *************************************************************************



% *************************************************************************
% STEP 2:                 -------------------------------------------------
interval_length = 24 * dt; % how many days in 1 interval
u_final         = [];
v_final         = [];
start_time      = datevec(datestart);
end_time        = datevec(dateend);
year_list       = start_time:1:end_time;

for file_loop = year_list(1):year_list(end)
    raw_wind_file       = [file_winddata 'winds_' num2str(file_loop)];
    raw_wind_data       = xlsread(raw_wind_file, 1); % read first worksheet in excel file containing raw wind data

    % select data range
    if file_loop == year_list(end)
        time_index	= max(find((raw_wind_data(:, 2) == end_time(2)) & (raw_wind_data(:, 3) == end_time(3))));
    else
        time_index	= length(raw_wind_data(:, 1));
    end
    
    direction           = raw_wind_data(1:time_index, 5);
    speed               = raw_wind_data(1:time_index, 6);
    
    % break into u and v vectors (modified from Rob Suryan)
    direction           = 90 - direction; % change to Math convention
    looky               = find(direction <= -180);
    direction(looky)    = 360 + direction(looky);
    direction           = direction * (pi / 180); % change to radians
    [u, v]              = pol2cart(direction, speed);
    v                   = v * -1;
    u                   = u * -1; % change vector sign
    
    % find interval median winds
    clms_reshaped       = time_index / interval_length;
    u_temp              = reshape(u, interval_length, clms_reshaped);
    v_temp              = reshape(v, interval_length, clms_reshaped);    
    u_median            = nanmedian(u_temp);
    v_median            = nanmedian(v_temp);
    
    % append medians onto growing wind vector columns
    u_final             = [u_final; u_median'];
    v_final             = [v_final; v_median'];
    
end

% map winds onto t_grid
wind_time               = t_grid; % time grid that wind data lies upon
lookNaN                 = find(isnan(u_final));
wind_time(lookNaN)      = [];
u_final(lookNaN)        = [];
v_final(lookNaN)        = [];% remove NaNs
u_final                 = interp1(wind_time, u_final, t_grid, 'pchip'); % map u_final onto t_grid
v_final                 = interp1(wind_time, v_final, t_grid, 'pchip'); % map v_final onto t_grid

% catch extreme values created during interpolation
lookerror               = find(abs(u_final) >= 15.5 | abs(v_final) >= 15.5); % error addresses in lookNaN
u_final(lookerror)      = 0; % (m/s); (vertical vector)
v_final(lookerror)      = 0; % (m/s); (vertical vector)
% *************************************************************************


% end m-file***************************************************************