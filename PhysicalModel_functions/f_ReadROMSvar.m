function val = f_ReadROMSvar(var, fileType)
%   ReadROMSvar Read ROMS data for the specified var from the appropriate file.

    % Infer whether we're working with UCSC or LiveOcean files depending on which file paths are defined
    if ~isempty(f_GetFilePath("wc12_gr")) & ~isempty(f_GetFilePath("LiveOcean"))
        error('Cannot define both UCSC and LiveOcean ROMS files in f_GetFilePath. Comment out the unused file paths.');
    elseif ~isempty(f_GetFilePath("wc12_gr"))
        ROMStype = 'UCSC';
    elseif ~isempty(f_GetFilePath("LiveOcean"))
        ROMStype = 'LiveOcean';
    else
        error('Paths to ROMS data files are not properly specified in f_GetFilePath');
    end

    % Read from UCSC data files
    if strcmp(ROMStype, 'UCSC')
        readFile_DepthLevels = f_GetFilePath("depth_levels_trimmed");
        readFile_grid = f_GetFilePath("wc12_gr");
        readFile_ExampleYear = f_GetFilePath("wc12_avg_2005_trimmed");

        % Determine which file contains the variable
        if any(ismember(var, ["lat_rho", "lon_rho", "x_rho", "y_rho", "mask_rho", "h", ...
                "lat_v", "lon_v", "x_v", "y_v", "mask_v", ...
                "lat_u", "lon_u", "x_u", "y_u", "mask_u", ...
                "lat_psi", "lon_psi", "x_psi", "y_psi", "mask_psi"]))
            fileType = "grid";
        elseif any(ismember(var, ["z_rho", "z_w"]))
            fileType = "depth";
        elseif any(ismember(var, ["lat_rho_flow", "lon_rho_flow"]))
            var = strrep(var, "_flow", "");
            fileType = "exampleYear";
        end

        % Read variable
        switch fileType
            case "grid"
                ncid = netcdf.open(readFile_grid, 'NC_NOWRITE');

            case "depth"
                ncid = netcdf.open(readFile_DepthLevels, 'NC_NOWRITE');

            case "exampleYear"
                ncid = netcdf.open(readFile_ExampleYear, 'NC_NOWRITE');
        end
        
        varid = netcdf.inqVarID(ncid, var);
        val = netcdf.getVar(ncid, varid);
        netcdf.close(ncid);

    end

end