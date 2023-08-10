% Calculate grid matrices (e.g., x_psi, y_psi) for LiveOcean ROMS files.
% Doug Jackson
% doug@QEDAconsulting.com

addpath("/Users/djackson/Documents/QEDA/NWFSC/ECOTRAN/programs/NCC2_ECOTRAN_CODE");

% Path to the myroms.org MATLAB scripts
myromsDir = f_GetFilePath("myromsDir");
addpath(myromsDir, fullfile(myromsDir, "grid"), fullfile(myromsDir, "netcdf"), ...
    fullfile(myromsDir, "utility"));

ROMSfiles = dir(fullfile(f_GetFilePath('ROMSdir'), "*.nc"));

for i = 1:numel(ROMSfiles)

    thisFile = ROMSfiles(i);
    thisPath = fullfile(thisFile.folder, thisFile.name);

    % Calculate using myroms.matlab.grid.roms_metrics
    thisROMSmetrics = roms_metrics(thisPath);

end