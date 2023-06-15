function filePath = f_GetFilePath(key)

    codePath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/NCC2_ECOTRAN_CODE";
    projectPath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/ECOTRANprojects/09jun23";

    filePaths = containers.Map();
    filePaths("codePath") = codePath;
    filePaths("projectPath") = projectPath;
    filePaths("PhysicalModel_functions") = fullfile(codePath, "PhysicalModel_functions");
    filePaths("ODEsolver_functions") = fullfile(codePath, "ODEsolver_functions");
    filePaths("NutrientFile_directory") = fullfile(codePath, "PhysicalModel_functions");
    filePaths("ERD_CUTI_directory") = fullfile(projectPath, "CUTI", "CUTI_daily.nc");

    % ROMS data files
    filePaths("wc12_gr") = "/Volumes/QEDA_ARCHIVE/QEDA/NWFSC/ECOTRAN/fromJim/wc12_grd.nc";
    filePaths("depth_levels_trimmed") = "/Volumes/QEDA_ARCHIVE/QEDA/NWFSC/ECOTRAN/fromJim/depth_levels_trimmed.nc";
    filePaths("ROMSdir") = "/Volumes/QEDA_ARCHIVE/QEDA/NWFSC/ECOTRAN/fromJim";
    filePaths("wc12_avg_2008_trimmed") = "/Volumes/QEDA_ARCHIVE/QEDA/NWFSC/ECOTRAN/fromJim/wc12_avg_gfdl_2005_trimmed.nc";


    filePath = filePaths(key);
end