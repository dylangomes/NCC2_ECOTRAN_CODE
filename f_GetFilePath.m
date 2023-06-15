function filePath = f_GetFilePath(key)

    codePath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/NCC2_ECOTRAN_CODE";
    projectPath = "/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/ECOTRANprojects/09jun23";

    filePaths = containers.Map();
    filePaths("NutrientFile_directory") = fullfile(codePath, "PhysicalModel_functions");
    filePaths("ERD_CUTI_directory") = fullfile(projectPath, "CUTI", "CUTI_daily.nc");

    filePath = filePaths(key);
end