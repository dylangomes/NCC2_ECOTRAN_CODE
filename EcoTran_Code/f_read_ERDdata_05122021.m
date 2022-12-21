function [ERDdata] = f_read_ERDdata_05122021(ERDparams)
% calculate daily medians from ERD products
% ERD products from: http://http://www.pfeg.noaa.gov/products/PFELData/upwell/6_hourly/upwell45N125W
% queried: 10-16-2015
% NOTE: raw units from PFEG are metric tons/sec/100 m. of coastline and this is converted to m3 within this function
% takes:
%       ERDparams.file_ERDdata
%       ERDparams.datestart
%       ERDparams.dateend
%       ERDparams.t_grid
% returns:
%       ERDdata.ERD_date                                              (vertical vector)
%       ERDdata.ERD_ekman_offshore      (m3/s per 100m of coastline); (vertical vector)
%       ERDdata.latitude                                              (scaler)
% calls:
%       none
% revision date: 05-12-2021



% *************************************************************************
% STEP 1: unpack variables-------------------------------------------------
file_ERDdata = ERDparams.file_ERDdata;
datestart    = ERDparams.datestart;
dateend      = ERDparams.dateend;
t_grid       = ERDparams.t_grid;
rho_water    = 1026;           % seawater density; (kg/m3)
% *************************************************************************




% *************************************************************************
% STEP 2: initialize variables---------------------------------------------
startvec        = datevec(datestart);
    year_start  = startvec(1);
    month_start = startvec(2);
    day_start   = startvec(3);
endvec          = datevec(dateend);
    year_end    = endvec(1);
    month_end   = endvec(2);
    day_end     = endvec(3);
year_list       = year_start:1:year_end;
lat_list        = 45; % latitudes
lat_loop        = lat_list(1); % lat_loop is a place holder for time when I put in several latitudes to run through
   
ERD_date        = datenum(year_start, month_start, day_start, 0, 0, 0):1:datenum(year_end, month_end, day_end, 0, 0, 0);
ERD_date        = ERD_date'; % transpose into column vector
result_length   = length(ERD_date);
num_t           = length(t_grid);

if (result_length ~= num_t)
    warning('WARNING: there is a mis-match in size of the time vector')
end

ERD_data(1:num_t, 1)	= NaN; % initialize as NaN; (vertical vector: num_t X 1)
% *************************************************************************





% *************************************************************************
% STEP 3: load PFEL-ERD product and calculate daily medians----------------
% step 3a: load data ------------------------------------------------------
fid          = fopen(file_ERDdata);
textformat   = ['%f %f %f %f %f ', '%*[^\r]']; % clm space string string string string space parameters...
indata       = textscan(fid, textformat, 'delimiter', ',', 'CollectOutput', 1);
fclose(fid);
raw_ERD_data = cell2mat(indata);

raw_year     = raw_ERD_data(:, 1);
raw_month    = raw_ERD_data(:, 2);
raw_day      = raw_ERD_data(:, 3);
raw_hour     = raw_ERD_data(:, 4);
raw_data     = raw_ERD_data(:, 5);

looky_BadData           = find(raw_data == -9999);
raw_data(looky_BadData) = NaN;
% -------------------------------------------------------------------------


% step 3b: convert data from (t/s per 100 m) to (m3/2 per 100 m) ----------
raw_data = raw_data * (1/rho_water) * 1000;
% -------------------------------------------------------------------------


% step 3c: ----------------------------------------------------------------
for year_loop = year_list(1):year_list(end)
	looky_year    = find(raw_year == year_loop);
    
	this_year  = raw_year(looky_year);
	this_month = raw_month(looky_year);
	this_day   = raw_day(looky_year);
	this_data  = raw_data(looky_year);    

    for month_loop = 1:12
        looky_month   = find(this_month == month_loop);
        current_month = this_month(looky_month);
    	current_day   = this_day(looky_month);
    	current_data  = this_data(looky_month);
    	start_day     = min(current_day); 
        end_day       = max(current_day);
        
        for day_loop = start_day:1:end_day
            looky_day = find(current_day == day_loop);
            if ~isempty(looky_day) % skip if there is no data for this date
                immediate_day           = current_day(looky_day);
                immediate_data          = current_data(looky_day);
                immediate_date          = datenum(year_loop, month_loop, day_loop, 0, 0, 0);
                looky_date              = find(ERD_date == immediate_date);
                ERD_data(looky_date, 1) = nanmedian(immediate_data); % (vertical vector: num_t X 1)
                
            end
        end % end day_loop
    end % end month_loop
end % end year_loop
% *************************************************************************





% *************************************************************************
% STEP 4: map data onto t_grid and interpolate missing data---------------
temp_time          = t_grid; % time grid that ERD data lies upon; (vertical vector: num_t X 1)
ERD_temp           = ERD_data; % (vertical vector: num_t X 1)
lookNaN            = find(isnan(ERD_temp)); 
temp_time(lookNaN) = []; 
ERD_temp(lookNaN)  = []; % remove NaNs
ERD_final          = interp1(temp_time, ERD_temp, t_grid, 'PCHIP'); % map data onto t_grid; (vertical vector: num_t X 1)
ERD_data           = ERD_final; % (vertical vector: num_t X 1)
% *************************************************************************





% *************************************************************************
% STEP 5: pack results
ERDdata.ERD_date                 = ERD_date; % ERD_date is in matlab time format; (vertical vector: num_t X 1)
ERDdata.ERD_ekman_offshore       = ERD_data; % (m3/s per 100m of coastline); (vertical vector: num_t X 1)
ERDdata.latitude                 = lat_list; % latitude (scaler)
% *************************************************************************


% end m-file***************************************************************