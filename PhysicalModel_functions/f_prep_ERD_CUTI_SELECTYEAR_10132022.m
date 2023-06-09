function ERD_CUTI = f_prep_ERD_CUTI_08052022(ERD_CUTI_input,CUTI_YEARS)
% prepare ERD_CUTI upwelling time-series
% by Jim Ruzicka
% takes: 
%       ERD_CUTI_input
%           readFile_ERD_CUTI	filename to load
%         	datestart            
%         	dateend              
%         	num_t                
%         	dt                  time step
%         	smoothing_window	moving average smoothing as time points for averaging before and after current time point
%           target_latitudes    choose 1 or more target latitude(s); [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
%
% returns:
%       ERD_CUTI                Coastal Upwelling Transport Index (vertical volume transport)'; (m2/s); 
%       ERD_latitude            latitudes; (2D matrix: num_lats X num_t)
%
% calls:
%   f_smooth
%
% revision: 8-5-2022


% *************************************************************************
% STEP 1: unpack operating parameters--------------------------------------
readFile_ERD_CUTI	= ERD_CUTI_input.readFile_ERD_CUTI;
datestart           = ERD_CUTI_input.datestart;
dateend             = ERD_CUTI_input.dateend;
num_t               = ERD_CUTI_input.num_t;
dt                  = ERD_CUTI_input.dt;
smoothing_window	= ERD_CUTI_input.smoothing_window; % set moving average smoothing as time points for averaging before and after current time point
target_latitudes	= ERD_CUTI_input.target_latitudes; % choose 1 or more target latitude(s); [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
% *************************************************************************


% *************************************************************************
% STEP 2: load ERD CUTI data-----------------------------------------------
ncid        = netcdf.open(readFile_ERD_CUTI, 'NC_NOWRITE'); % Open readFile
% ncid        = netcdf.open('CUTI/CUTI_daily - Copy.nc', 'NC_NOWRITE'); % Open readFile; For TESTING ONLY

varid       = netcdf.inqVarID(ncid, 'year');
CUTI_year	= netcdf.getVar(ncid, varid); % year; (vertical vector: 12600 X 1)

varid       = netcdf.inqVarID(ncid, 'month');
CUTI_month	= netcdf.getVar(ncid, varid); % month; (vertical vector: 12600 X 1)

varid       = netcdf.inqVarID(ncid, 'day');
CUTI_day	= netcdf.getVar(ncid, varid); % day; (degrees east); (vertical vector: 12600 X 1)

varid       = netcdf.inqVarID(ncid, 'time');
CUTI_time	= netcdf.getVar(ncid, varid); % time; (days since 1970-01-01); (vertical vector: 12600 X 1)

varid       = netcdf.inqVarID(ncid, 'latitude');
ERD_latitude	= netcdf.getVar(ncid, varid); % 'latitude'; (degrees east); (vertical vector: 12600 X 1)


if CUTI_YEARS == "AVG"
Avg = readtable('CUTI/CUTI_AVERAGE_daily.csv');
ERD_CUTI	= [table2array(Avg(:,5:21))]'; % 'Coastal Upwelling Transport Index (vertical volume transport)'; (m2/s); (2D matrix: 17 X 12600)

disp('  using an average CUTI time-series (1988-2021 avg, repeated)')

elseif CUTI_YEARS == "ALL"
varid       = netcdf.inqVarID(ncid, 'CUTI');
ERD_CUTI	= netcdf.getVar(ncid, varid); % 'Coastal Upwelling Transport Index (vertical volume transport)'; (m2/s); (2D matrix: 17 X 12600)
ERD_CUTI=ERD_CUTI(:,1:12419);

disp('using a real CUTI time-series (1988-2021, repeated if needed)')

else
Select = readtable('CUTI/CUTI_daily.csv');
DAYS=find(Select.year==str2num(CUTI_YEARS));

ERD_CUTI	= [table2array(Select(min(DAYS):max(DAYS),4:20))]'; % 'Coastal Upwelling Transport Index (vertical volume transport)'; (m2/s); (2D matrix: 17 X 12600)

ERD_CUTI=repmat(ERD_CUTI,35);
ERD_CUTI=ERD_CUTI(:,1:12419);

disp(strcat('using a repeated CUTI time-series (',CUTI_YEARS,")"))

end


% to trim 2022 off, because it is incomplete year
CUTI_year=CUTI_year(1:12419,:);
CUTI_month=CUTI_month(1:12419,:);
CUTI_day=CUTI_day(1:12419,:);
CUTI_time=CUTI_time(1:12419,:);


% Close the NetCDF file
netcdf.close(ncid)
% *************************************************************************



% *************************************************************************
% STEP 3: select ERD_CUTI time period--------------------------------------
%         extend time-series as necessary to reach required dateend

% step 3a: build matlab time vector ---------------------------------------
CUTI_date	= datenum(CUTI_year, CUTI_month, CUTI_day); % CUTI time in matlab format
num_t_CUTI	= length(CUTI_date);
% -------------------------------------------------------------------------


% step 3b: trim ERD_CUTI for time-series earlier than datestart -----------
if datestart < CUTI_date
    error('ERROR: chosen start date precedes available ERD_CUTI time-series')
else
    looky_early             = find(CUTI_date < datestart);
    CUTI_date(looky_early)	= [];
    num_t_CUTI              = length(CUTI_date);
    ERD_CUTI(:, looky_early)    = [];
end
% -------------------------------------------------------------------------


% step 3c: find desired end of CUTI time-series ---------------------------
if dateend > CUTI_date(end)
    
    disp('dateend is later than available ERD_CUTI time-series')
    disp('  replicating ERD_CUTI time-series to expand to dateend')
    
    num_replcates           = ceil(num_t ./ num_t_CUTI) + 1;
    
    ERD_CUTI              	= repmat(ERD_CUTI, [1, num_replcates]); % (2D matrix: num_lats X num_t_CUTI)
    
    [~, num_t_CUTI]      	= size(ERD_CUTI); % new time-series length
    CUTI_date               = CUTI_date(1):dt:(CUTI_date(1)+num_t_CUTI-1); % (horizontal vector: 1 X num_t_CUTI)
    CUTI_date               = CUTI_date'; % transpose; (vertical vector: num_CUTI_time X 1)
    
end
% -------------------------------------------------------------------------


% step 3d: find desired end of CUTI time-series ---------------------------
looky_late                  = find(CUTI_date > dateend);
ERD_CUTI(:, looky_late)     = [];
CUTI_date(looky_late)       = [];
% *************************************************************************



% *************************************************************************
% STEP 4: select target latitude(s)----------------------------------------
ERD_CUTI                    = ERD_CUTI'; % (2D matrix: num_t_CUTI X num_lats)
looky_latitude              = find(ismember(ERD_latitude, target_latitudes));
ERD_CUTI                	= ERD_CUTI(:, looky_latitude); % (m3/s per 1m of BoxWidth); % (2D matrix: num_t_CUTI X num_target_lats)
% *************************************************************************



% *************************************************************************
% STEP 5: optional moving average smoothing--------------------------------
if smoothing_window > 0
    disp(['NOTE: ERD_CUTI upwelling time-series is smoothed in window = ' num2str(smoothing_window)])
    ERD_CUTI	= f_smooth(ERD_CUTI, smoothing_window);
end

if length(ERD_CUTI) ~= num_t
    warning('WARNING: ERD_CUTI and t_grid have different time lengths');
end % (length(ERD_CUTI) ~= num_t)
% *************************************************************************


% end m-file***************************************************************