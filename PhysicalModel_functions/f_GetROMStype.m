function ROMStype = f_GetROMStype()
% f_GetROMStype Infer whether we're working with UCSC or LiveOcean files depending on which file paths are defined.

    if ~isempty(f_GetFilePath("wc12_gr")) & ~isempty(f_GetFilePath("LiveOceanGrid"))
        error('Cannot define both UCSC and LiveOcean ROMS files in f_GetFilePath. Comment out the unused file paths.');
    elseif ~isempty(f_GetFilePath("wc12_gr"))
        ROMStype = 'UCSC';
    elseif ~isempty(f_GetFilePath("LiveOceanGrid"))
        ROMStype = 'LiveOcean';
    else
        error('Paths to ROMS data files are not properly specified in f_GetFilePath');
    end

end