% ECOTRAN_DynamicScenarios_NCC2_2D_12082022
% Use this to batch run ECOTRAN models under a defined range of conditions
% Output files and metadata will be saved in a new "Output" folder within the working directory you define below (setWD)

%% SETUP MODEL
setWD = 'C:/Users/dgome/Documents/PROJECT_DIRECTORY/';% Choose directory (project)
cd(setWD)
Model_name = 'MODEL_FILE.csv'; % Choose model csv 

START = '01-Jan-1998' % choose time start
END = '01-Jan-1999'  % choose time end

Region = 'CR'
CUTI_LAT = "Auto"; % Can manually override automatic settings ("Auto" -> WA=47,CR=46,NOR=45,SOR=43,NCA=41) by choosing 1 target latitude for upwelling [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]

% set up scenarios
Treatment_1 = [9 10]; % vector; numbers represent Ecotran functional groups to change
Treatment_2 = [1.2];     % vector; factor to scale above functional groups by, >1 will increase groups and <1 (but > 0) will decrease groups
% Treatment_3 = ;

%% MODEL SWITCHES
switch_FoodWebScenario      = 'FoodWebScenario_ON';
% switch_FoodWebScenario       = 'FoodWebScenario_OFF';

% switch_SubModel             = 'identical_SubModel';                 % OPTION 1: use the same food web across the shelf
switch_SubModel             = 'independent_SubModel';               % OPTION 2: use independently defined webs for each shelf zone

% switch_INITIALproduction	= 'INITIALproduction_MichaelisMenton';	% METHOD 1: use for driving initial model conditions with primary production defined by Michaelis-Menton uptake
% switch_INITIALproduction	= 'INITIALproduction_pb';               % METHOD 2: use for driving initial model conditions with primary production defined by p = [(p/b) * b]
switch_INITIALproduction	= 'INITIALproduction_SubModel';         % METHOD 3: use for driving initial model conditions with values loaded along with regional sub-model definitions 
% switch_INITIALproduction	= 'INITIALproduction_nutrients';      	% METHOD 4: use for driving initial model conditions with mean annual nutrient input rates
%% END MODEL SETUP AND SWITCHES (no editing required below)



% create output folder & file_code_counter, based on previous model runs
MODEL = regexprep(Model_name,".csv","");
if ~exist(strcat(setWD,"Output"), 'dir');% checks for "Output" folder, and 
    mkdir(strcat(setWD,"Output"));% creates it if it doesn't exist
    
    file_code_counter = 0; % where to start file name number code
else
    S = dir(strcat(setWD,"Output/",MODEL,'*',date,'.mat')); % check if files exist already
    if size(S,1)>0; % if files exist...
        FILE = [];
        for i = 1:size(S,1);
    ADD = regexprep(S(i).name,strcat(MODEL,'_(\d{3})_',date,'.mat'),'$1'); % grab file code number
    FILE = [FILE; ADD];
        end
        file_code_counter = max(str2num(FILE)); % where to start file name number code if files already exist (take the maximum file number for a given model and date)
    else
        file_code_counter = 0; % where to start file name number code
    end
end

% STEP 2: batch process all variable combinations--------------------------
num_Treatment_1     = length(Treatment_1);
num_Treatment_2     = length(Treatment_2);
% num_Treatment_3     = length(Treatment_3);
num_Treatment_total	= num_Treatment_1 * num_Treatment_2 %* num_Treatment_3;

counter = 0;
z=[]; % initialize array to save metadata

for Treatment_1_loop = 1:num_Treatment_1
	current_Treatment_1 = Treatment_1(Treatment_1_loop);
    
    for Treatment_2_loop = 1:num_Treatment_2
        current_Treatment_2 = Treatment_2(Treatment_2_loop);
%         
%         for Treatment_3_loop = 1:num_Treatment_3
%             current_Treatment_3 = Treatment_3(Treatment_3_loop);
            
            run_Treatments      = [current_Treatment_1, current_Treatment_2]; % choose functional group (Treatment_1) to change by a factor of X (Treament 2); 
%           run_Treatments      = [current_Treatment_1, current_Treatment_2 current_Treatment_3]; % can add a 3rd treatment by un-commenting treatment 3 code and adding "current_Treatment_3" here
       
                
            file_code_counter	= file_code_counter+1; % simulation filename number
            counter = counter+1; % simulation number
            z = [z; file_code_counter	run_Treatments CUTI_LAT strcat(START,":",END) Model_name Region switch_FoodWebScenario]; % to save the metadata
            
            disp(' ')
            disp(' ')
            disp('-----------------------------')
            disp(['Now running simulation: ' num2str(counter) ' of ' num2str(num_Treatment_total)])
            disp('-----------------------------')
            disp(' ')
            disp(' ')
            

            ECOTRAN_DynamicCode_NCC2_12082022(setWD,Model_name,START,END,Region,CUTI_LAT,file_code_counter,run_Treatments,switch_FoodWebScenario,switch_SubModel,switch_INITIALproduction);

% 
%         end % (Treatment_3_loop)
    end % (Treatment_2_loop)
end % (Treatment_1_loop)
    
    
writematrix(z, strcat("Output/","MetaData_",date,".csv"),'WriteMode','append')

% end m-file***************************************************************