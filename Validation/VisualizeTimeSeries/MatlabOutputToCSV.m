% change to match path
PATH="Output/"
Files=dir(strcat(PATH,"*.mat"));


for i = 1:size(Files,1)
    if exist(strcat(PATH,regexprep(Files(i).name,".mat","_SR1.csv")))==0 % check if file exists in Step2 already, if not proceed:
FileData = load(strcat(PATH,Files(i).name)); % load .mat file

csvwrite(strcat(PATH,regexprep(Files(i).name,".mat","_SR1.csv")),... % rewrite part of file as csv
    FileData.re_Y(:,:,1)); % for 2D model, select subregion in 3rd position. options are 1:5
    end
end
