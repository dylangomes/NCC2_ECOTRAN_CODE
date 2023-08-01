function filePath = f_GetFilePath(key)

    projectPath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/ECOTRANprojects/ROMS_24jul23";
    codePath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/NCC2_ECOTRAN_CODE";
    ROMSdir = "/Volumes/QEDA_ARCHIVE/QEDA/NWFSC/ECOTRAN/fromJim";

    filePaths = containers.Map();
    filePaths("codePath") = codePath;
    filePaths("projectPath") = projectPath;
    filePaths("PhysicalModel_functions") = fullfile(codePath, "PhysicalModel_functions");
    filePaths("ODEsolver_functions") = fullfile(codePath, "ODEsolver_functions");
    filePaths("NutrientFile_directory") = fullfile(codePath, "PhysicalModel_functions");
    filePaths("ERD_CUTI_directory") = fullfile(projectPath, "CUTI", "CUTI_daily.nc");

    % ROMS data files
    filePaths("ROMSdir") = ROMSdir;
    filePaths("wc12_gr") = fullfile(ROMSdir, "wc12_grd.nc");
    filePaths("depth_levels_trimmed") = fullfile(ROMSdir, "depth_levels_trimmed.nc");
    filePaths("wc12_avg_2005_trimmed") = fullfile(ROMSdir, "wc12_avg_gfdl_2005_trimmed.nc");

    % Location of pre-processed ROMS data
    filePaths("preproDir") = fullfile(projectPath, "prepro");

    filePath = filePaths(key);
end