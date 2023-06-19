% ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023
% Use this to batch run ECOTRAN models under a defined range of conditions
% Output files and metadata will be saved in a new "Output/temp/" folder within the working directory you define below (setWD)
% User must add ECOTRAN_NCC2_package and all subdirectories to search path so functions can be located
%
% calls:
%       f_ECOTRANdynamic_NCC2_PARALLEL_05152023     functionalalized version of main ECOTRANdynamic script; use for batch runs; PARALLEL COMPUTING VERSION
% 
% revision date: 5/15/2023
%       5/15/2023  cleaned up indentation; changed f_ECOTRANdynamic function call f_ECOTRANdynamic_NCC2_PARALLEL_05152023 (revision of ECOTRAN_DynamicCode_Parallel_NCC2_02172023)

%% SETUP MODEL
tStart              = tic
setWD               = 'C:/Users/dgome/Documents/NCC2_ECOTRAN_CODE/';% Choose directory (project); you can create of copy of ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023 in a new project directory (setWD) to change input parameters without altering ECOTRANdynamic_NCC2_TEMPLATE_PARALLEL_05152023 in the shared GitHub repo.
cd(setWD)

Model_name          = 'NCC2_09032022.csv'; % Choose model csv 

ShowOutput          = 0;       % 0 = suppress some function notes output (saves computation time), 1 = show output; Can suppress more later, tedious
num_MC              = 1;     % SSS set this value to number of Monte Carlo models; test with value of 1; if switch_MonteCarlo = 'MonteCarlo_TypeModel' then this number doesn't matter
FileOffset          = 0;       % set this if adding to other files. e.g., if already ran 10 files of this type, and want the MC number to start from 11, set this to 10.
current_MC_offset   = 0; 

START               = '01-Jan-2000' % choose time start
END                 = '01-Jan-2001'  % choose time end; test with a smaller timeseries before full run
% END                 = '01-Jan-2075'  % choose time end

Region              = 'NA'; % Can choose a specific region to model (instead of the fully aggregated model): 'CR' = Columbia River, 'WA' = Washington coast, 'NOR' = Northern Oregon, 'SOR'=Southern Oregon, 'NCA'= Northern California
CUTI_YEARS          = "AVG"; % "ALL" "AVG" OR one year: "2008"
CUTI_LAT            = 45; % choose a latitude for CUTI data (may not matter for AVG_CUTI); If a specific Region is selected (above), can choose "Auto", which will result in auto-selecting latitudes -> WA=47,CR=46,NOR=45,SOR=43,NCA=41. Full list of target latitude options for upwelling: [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]

% set up scenarios
Treatment_1     = [8 10];               % vector; numbers represent Ecotran functional groups to change e.g., [9 10] is to change group 9 and then group 10 in separate scenarios (each will be run through the number of MC models and all Treatment_2 scenarios)
Treatment_2     = [1.2];                % vector; factor to scale above functional groups by, >1 will increase groups and <1 (but > 0) will decrease groups; this can be a vector of multiple values e.g., [0.25 0.5 0.75] which will run separate scenarios for each Treatment_1 specified
% Total scenarios run = ;ength of Treatment_1 x length of Treatment_2 (see num_Treatment_total below)
% Total simulations run = total scenarios run x num_MC


%% MODEL SWITCHES
switch_FoodWebScenario      = 'FoodWebScenario_ON'; % allow above treatments to be applied
% switch_FoodWebScenario       = 'FoodWebScenario_OFF';

switch_SubModel             = 'identical_SubModel';                 % OPTION 1: use the same food web across the shelf
% switch_SubModel             = 'independent_SubModel';               % OPTION 2: use independently defined webs for each shelf zone

% switch_INITIALproduction	= 'INITIALproduction_MichaelisMenton';	% METHOD 1: use for driving initial model conditions with primary production defined by Michaelis-Menton uptake
% switch_INITIALproduction	= 'INITIALproduction_pb';               % METHOD 2: use for driving initial model conditions with primary production defined by p = [(p/b) * b]
% switch_INITIALproduction	= 'INITIALproduction_SubModel';         % METHOD 3: use for driving initial model conditions with values loaded along with regional sub-model definitions 
% switch_INITIALproduction	= 'INITIALproduction_nutrients';      	% METHOD 4: use for driving initial model conditions with mean annual nutrient input rates
switch_INITIALproduction	= 'INITIALproduction_BioGeoChemicalModel';	% METHOD 5: use for driving initial model conditions with mean annual ROMS-BioGeoChemichal Primary Production rates

switch_MonteCarlo           = 'MonteCarlo_build';	% generate (and optionally save) a stack of MonteCarlo food webs
% switch_MonteCarlo           = 'MonteCarlo_load';	% load a saved stack of MonteCarlo food webs
% switch_MonteCarlo           = 'MonteCarlo_TypeModel';	% use NO MonteCarlo food webs


%% END MODEL SETUP AND SWITCHES (no editing required below)

% add paths to other functions within ECOTRAN_NCC2_package
folder = fileparts(mfilename('fullpath')); 
addpath(genpath(folder));


if ~exist(strcat("Output"), 'dir')	% checks for "Output" folder, and 
    mkdir(strcat("Output"));          % creates it if it doesn't exist
end
    if ~exist(strcat("Output/temp"), 'dir')	% checks for "temp" folder (to store outputs before they are organized by the user), and 
        mkdir(strcat("Output/temp"));          % creates it if it doesn't exist
    end
    
upwelling_driver = 'ERD_CUTI'; % this is the actual CUTI timeseries 

% STEP 2: batch process all variable combinations--------------------------
num_Treatment_1     = length(Treatment_1);
num_Treatment_2     = length(Treatment_2);
% num_Treatment_3     = length(Treatment_3);
num_Treatment_total	= num_Treatment_1 * num_Treatment_2; %* num_Treatment_3;

% counter = 0;
% z=[]; % initialize array to save metadata

% parallelize (parfor) loop with longest length (Treatment_1 vs Treatment_2).
if num_Treatment_1 >= num_Treatment_2
    parfor Treatment_1_loop = 1:num_Treatment_1

        current_Treatment_1     = Treatment_1(Treatment_1_loop);
    
        for Treatment_2_loop = 1:num_Treatment_2
            current_Treatment_2     = Treatment_2(Treatment_2_loop);
            run_Treatments          = [current_Treatment_1, current_Treatment_2]; % choose functional group (Treatment_1) to change by a factor of X (Treament 2); 
            rng('shuffle');

            f_ECOTRANdynamic_NCC2_PARALLEL_05152023(setWD, Model_name, START, END, Region, upwelling_driver, ...
                                                    CUTI_LAT, CUTI_YEARS, run_Treatments, switch_FoodWebScenario, ...
                                                    switch_SubModel, switch_INITIALproduction, switch_MonteCarlo, ...
                                                    num_MC, FileOffset, current_MC_offset, ShowOutput);

        end % (Treatment_2_loop)
    end % (Treatment_1_loop)
else % this else statement parallelizes the other loop, if the Treatment_2 vector is longer
    for Treatment_1_loop = 1:num_Treatment_1

        current_Treatment_1     = Treatment_1(Treatment_1_loop);
    
        parfor Treatment_2_loop = 1:num_Treatment_2
            current_Treatment_2     = Treatment_2(Treatment_2_loop);
            run_Treatments          = [current_Treatment_1, current_Treatment_2]; % choose functional group (Treatment_1) to change by a factor of X (Treament 2); 
            rng('shuffle');

            f_ECOTRANdynamic_NCC2_PARALLEL_05152023(setWD, Model_name, START, END, Region, upwelling_driver, ...
                                                    CUTI_LAT, CUTI_YEARS, run_Treatments, switch_FoodWebScenario, ...
                                                    switch_SubModel, switch_INITIALproduction, switch_MonteCarlo, ...
                                                    num_MC, FileOffset, current_MC_offset, ShowOutput);

        end % (Treatment_2_loop)
    end % (Treatment_1_loop)
end % (if num_Treatment_1 >= num_Treatment_2)

% writematrix(z, strcat("Output/","MetaData_",date,".csv"),'WriteMode','append') % to save metadata
toc(tStart)
% *************************************************************************


% end m-file***************************************************************