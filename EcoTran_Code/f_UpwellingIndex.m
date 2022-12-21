function [UpwellingIndex] = f_UpwellingIndex(UpwellingParams)
% calculate the Bakun upwelling index from Rian Hoof's NWPO3 wind data
% formulae used here are from:
%   COASTAL UPWELLING INDICES
%   WEST COAST OF NORTH AMERICA
%   1946-95
%   NOAA-TM-NMFS-SWFSC-231 (1996)
%   Franklin B. Schwing
%   Michael OÌFarrell
%   John M. Steger
%   Kenneth Baltz
% takes:
%       UpwellingParams.northerlywinds  
%       UpwellingParams.easterlywinds
%       UpwellingParams.rho_air         air density; (kg/m3) (Schwing et al. 1996)
%       UpwellingParams.rho_water       seawater density; (kg/m3)
%       UpwellingParams.coriollis       coriollis factor; (1/s); (NH-Line NNPPZD model: Ruzicka et al. 2011)
% calls:
%       none
% revision date: 3-13-2014

% *************************************************************************
% STEP 1: unpack variables-------------------------------------------------
northerlywinds = UpwellingParams.northerlywinds;
easterlywinds  = UpwellingParams.easterlywinds;
rho_air        = UpwellingParams.rho_air; % air density; (kg/m3) (Schwing et al. 1996)
rho_water      = UpwellingParams.rho_water; % seawater density; (kg/m3)
coriollis      = UpwellingParams.coriollis; % coriollis factor (1/s); (NH-Line NNPPZD model: Ruzicka et al. 2011)





% *************************************************************************
% STEP 2: non-directional wind speed---------------------------------------
W   = sqrt(northerlywinds.^2 + easterlywinds.^2); % wind speed (m/s)





% *************************************************************************
% STEP 3: drag coefficient-------------------------------------------------
Cd(1:length(W), 1) = NaN; % initialize drag coefficient vector; (dimensionless)

for loop1 = 1:length(W);
    if W(loop1) <= 1
        Cd(loop1) = 0.00218;
    elseif (W(loop1) > 1) && (W(loop1) < 3)
        Cd(loop1) = (0.62 + 1.56 / W(loop1)) * 0.001;
    elseif W(loop1) >= 3 && W(loop1) < 10
        Cd(loop1) = 0.00114;
    elseif W(loop1) >= 10
        Cd(loop1) = (0.49 + 0.065 * W(loop1)) * 0.001;
    end
end

% Cd estimation method of PFEG
% The parameterization of the drag coefficient (CD) used in the wind stress calculation is taken as a function of wind speed (W) rather than a constant. :
%     * Historically, CD was taken as a constant value of 0.0013 to calculate indices from 6-hourly pressure fields and was increased to 0.0026 to calculate indices from a monthly pressure field.
%     * For the calculations on this page, a non-linear drag coefficient was used (based on Large and Pond (1981) modified for low wind speeds as in Trenberth et al. (1990)):
%       CD      = 0.00218 for W < or = 1m/s
%           = (0.62+1.56/W) x .001 for 1m/s < W < 3m/s
%           = 0.00114 for 3m/s < or = W < 10m/s
%           = (0.49 + 0.065W) x .001 for W > or = 10m/s
%           o References: Large, W., and S. Pond, 1981: Open Ocean momentum flux measurements in moderate to strong winds.  J. Phys. Oceanogr., 11, 324-336.
%           o Trenberth, K. E., W. G. Large and J. G. Olson, 1990: The mean annual cycle in global ocean wind stress. J. Phys. Oceanogr., 20, 1742-1760.
%
% Alternative Cd:
% Cd = 0.0013; % empirical drag coefficient (dimensionless) Bakun 1987
                % use drag = 0.0026 with monthly mean pressure data
                % (monthly mean geostrophic winds)

                
                
                
                
% *************************************************************************
% STEP 4: wind stress------------------------------------------------------
tau_y = rho_air * Cd .* northerlywinds .* W; % north-south wind stress; (kg/m/s2); (vertical vector)





% *************************************************************************
% STEP 5: upwelling index--------------------------------------------------
BUI = (tau_y / coriollis); % Ekman transport, negative offshore (kg/m/s)
BUI = BUI * (1/rho_water) * 100; % convert kg to mt and scale up to 100m of coastline; (mt/s per 100m coastline) or (m3/s per 100m coastline)
BUI = -BUI; % Bakun Upwelling Index (positive means upwelling)





% *************************************************************************
% STEP 6: pack results for export------------------------------------------
UpwellingIndex.BUI   = BUI; % (m3/s per 100m of coastline); (vertical vector)
UpwellingIndex.tau_y = tau_y; % north-south wind stress; (kg/m/s2); (vertical vector)

% end m-file***************************************************************