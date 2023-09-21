function filename = f_ReadROMSyear(var, baseFilename)
%   f_GetROMSfilename Determine the appropriate ROMS filename based on var and ROMStype.

    ROMStype = f_GetROMStype();
    
    
    % USCS yearly ROMS data are stored in a single file => just return baseFilename
    if strcmp(ROMStype, "UCSC")
        filename = baseFilename;
        
    elseif strcmp(ROMStype, "LiveOcean")
        
        % Determine the appropriate filename according to var
        if any(ismember(var, ["lat_rho", "lon_rho", "s_rho", "zeta", "phytoplankton", "hc", "ocean_time"]))
            suffix = "phytoplankton";
        elseif any(ismember(var, ["temp"]))
            suffix = "temp";
        elseif any(ismember(var, ["u", "lat_u", "lon_u"]))
            suffix = "u";
        elseif any(ismember(var, ["v", "lat_v", "lon_v"]))
            suffix = "v";
        elseif any(ismember(var, ["w", "s_w"]))
            suffix = "w";
        else
            error('Unrecognized variable in f_GetROMSfilename');
        end
        
        filename = strrep(baseFilename, '.nc', strcat('_', suffix, '.nc'));
        
    end
end