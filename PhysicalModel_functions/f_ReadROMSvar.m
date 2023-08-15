function val = f_ReadROMSvar(var, fileType)
%   f_ReadROMSvar Read ROMS data for the specified var from the appropriate file.

    ROMStype = f_GetROMStype();

    % Read from LiveOcean files
    if strcmp(ROMStype, "LiveOcean")
        gridFile = f_GetFilePath("LiveOceanGrid");
        metricsFile = f_GetFilePath("LiveOceanMetrics");
        exampleYearFile = f_GetFilePath("LiveOceanExampleYear");

        % Determine which file contains the variable
        if any(ismember(var, ["x_psi", "y_psi"]))
            fileType = "metrics";
        elseif any(ismember(var, ["lat_rho_flow", "lon_rho_flow"]))
            var = strrep(var, "_flow", "");
            fileType = "exampleYear";
        else
            fileType = "grid";
        end

        % Read variable
        switch fileType
            case "grid"
                ncid = netcdf.open(gridFile, 'NC_NOWRITE');
                varid = netcdf.inqVarID(ncid, var);
                val = netcdf.getVar(ncid, varid);
                netcdf.close(ncid);
            case "metrics"
                ROMSmetrics = load(metricsFile, "ROMSmetrics");
                val = ROMSmetrics.ROMSmetrics.(var);

            case "exampleYear"
                ncid = netcdf.open(exampleYearFile, 'NC_NOWRITE');
                varid = netcdf.inqVarID(ncid, var);
                val = netcdf.getVar(ncid, varid);
                netcdf.close(ncid);
        end

    % Read from UCSC data files
    elseif strcmp(ROMStype, "UCSC")
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
                varid = netcdf.inqVarID(ncid, var);
                val = netcdf.getVar(ncid, varid);
                netcdf.close(ncid);

            case "depth"
                ncid = netcdf.open(readFile_DepthLevels, 'NC_NOWRITE');
                varid = netcdf.inqVarID(ncid, var);
                val = netcdf.getVar(ncid, varid);
                netcdf.close(ncid);

            case "exampleYear"
                ncid = netcdf.open(readFile_ExampleYear, 'NC_NOWRITE');
                varid = netcdf.inqVarID(ncid, var);
                val = netcdf.getVar(ncid, varid);
                netcdf.close(ncid);
        end
    end

end