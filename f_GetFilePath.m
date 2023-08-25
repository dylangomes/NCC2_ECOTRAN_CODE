function filePath = f_GetFilePath(key)

    projectPath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/ECOTRANprojects/ROMS_24jul23";
    codePath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/NCC2_ECOTRAN_CODE";

    filePaths = containers.Map();
    filePaths("codePath") = codePath;
    filePaths("projectPath") = projectPath;
    filePaths("PhysicalModel_functions") = fullfile(codePath, "PhysicalModel_functions");
    filePaths("ODEsolver_functions") = fullfile(codePath, "ODEsolver_functions");
    filePaths("NutrientFile_directory") = fullfile(codePath, "PhysicalModel_functions");
    filePaths("ERD_CUTI_directory") = fullfile(projectPath, "CUTI", "CUTI_daily.nc");

    % *************************************************************************
    % ROMS data files
    
    % Location of pre-processed ROMS data
    filePaths("preproDir") = fullfile(projectPath, "prepro");

    % NOTE: COMMENT OUT EITHER UCSC OR LIVEOCEAN PATHS DEPENDING ON WHICH FILE TYPES YOU'RE USING

    % USCS paths
    % ROMSdir = "/Volumes/QEDA_ARCHIVE/QEDA/NWFSC/ECOTRAN/fromJim";
    % filePaths("wc12_gr") = fullfile(ROMSdir, "wc12_grd.nc");
    % filePaths("depth_levels_trimmed") = fullfile(ROMSdir, "depth_levels_trimmed.nc");
    % filePaths("wc12_avg_2005_trimmed") = fullfile(ROMSdir, "wc12_avg_gfdl_2005_trimmed.nc");

    % LiveOcean paths
    ROMSdir = "/Volumes/QEDA_ARCHIVE/QEDA/NWFSC/ECOTRAN/LiveOcean";
    filePaths("LiveOceanGrid") = fullfile(ROMSdir, "LiveOcean_2017_grid.nc");
    filePaths("LiveOceanMetrics") = fullfile(ROMSdir, "LiveOcean_2017_ROMSmetrics.mat");
    filePaths("LiveOceanExampleYear") = fullfile(ROMSdir, "LiveOcean_2018.nc");
    
    % END OF SECTION TO COMMENT/UNCOMMENT

    filePaths("ROMSdir") = ROMSdir;

    % *************************************************************************
    % LiveOcean pre-processing paths
    filePaths("myromsDir") = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/matlab";

    if isKey(filePaths, key)
        filePath = filePaths(key);
    else
        filePath = [];
    end

end