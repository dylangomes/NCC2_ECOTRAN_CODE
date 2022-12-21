% ECOTRAN_TuningMacro_NCC2_2D_09212022
% Use this macro to batch run ECOTRAN models under a defined range of conditions
% NOTE: the ECOTRANdynamic script/function needs to be set up to receive and interpret the variable operating conditions

setWD = 'C:/Users/dgome/Documents/PROJECT_DIRECTORY/';% Choose directory (project)
cd(setWD)
Model_name = 'MODEL_FILE.csv'; % Choose model csv 

START = '01-Jan-1998' % choose time start
END = '01-Jan-1999'  % choose time end

CUTI_LAT = 45; % choose 1 target latitude for upwelling [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]

switch_FoodWebScenario      = 'FoodWebScenario_ON';
% switch_FoodWebScenario       = 'FoodWebScenario_OFF';
Treatment_1 = [9 10]; % numbers represent Ecotran functional groups
Treatment_2 = 1.2;     % factor to scale above functional groups by. >1 will increase groups and <1 (but > 0) will decrease groups
% Treatment_3 = treatment_benthiccDetritusMetabolism;


% STEP 2: batch process all variable combinations--------------------------
num_Treatment_1     = length(Treatment_1);
num_Treatment_2     = length(Treatment_2);
% num_Treatment_3     = length(Treatment_3);
num_Treatment_total	= num_Treatment_1 * num_Treatment_2 %* num_Treatment_3;

file_code_counter = 0; % QQQ where to start file name number code
z=[]; % initialize array to save metadata

for Treatment_1_loop = 1:num_Treatment_1
	current_Treatment_1 = Treatment_1(Treatment_1_loop);
    
    for Treatment_2_loop = 1:num_Treatment_2
        current_Treatment_2 = Treatment_2(Treatment_2_loop);
%         
%         for Treatment_3_loop = 1:num_Treatment_3
%             current_Treatment_3 = Treatment_3(Treatment_3_loop);
            
            run_Treatments      = [current_Treatment_1, current_Treatment_2]; % choose functional group to change
            file_code_counter	= file_code_counter+1; % simulation filename number
            z = [z; file_code_counter	run_Treatments CUTI_LAT strcat(START,":",END) Model_name switch_FoodWebScenario]; % to save the metadata
            
            disp(' ')
            disp(' ')
            disp('-----------------------------')
            disp(['Now running simulation: ' num2str(file_code_counter) ' of ' num2str(num_Treatment_total)])
            disp('-----------------------------')
            disp(' ')
            disp(' ')
            

            ECOTRAN_DynamicCode_NCC2_12082022(setWD,Model_name,START,END,CUTI_LAT,file_code_counter,run_Treatments,switch_FoodWebScenario);

% 
%         end % (Treatment_3_loop)
    end % (Treatment_2_loop)
end % (Treatment_1_loop)

writematrix(z, strcat("Output\","MetaData.csv"))

% end m-file***************************************************************