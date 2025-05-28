%% ICC Values for Within Session Functional Connectivity and Brain Variability
%% ICC Calculation and Analysis (Figures 2, 3, 6)

% This file contains code to calculate within session ICC values (Part 1) and
% complete statistical analyses from figures 2, 3 and 6 (Part 2) from the paper. The
% code is meant to be run section by section IN ORDER, as there are
% dependencies between sections. Make sure to complete running Part 1
% before you move on to Part 2.

%% Add data paths

% Add paths to working directories.
data_dir = '/YOUR/PATH/HERE'; % Replace with the path to your CONN output .mat files are stored
output_dir = '/YOUR/PATH/HERE'; % Replace with where the output Excel files from this script will be stored.

% Add in the filenames of the CONN output files for all networks. 
results_whole = 'results_whole_ses01.mat'; % whole brain network
results_motor = 'results_motor_ses01.mat'; % motor network
results_limbic = 'results_limbic_ses01.mat'; % limbic network
results_associative = 'results_assoc_ses01.mat'; % associative network

% Add paths to functions.
addpath('YOUR/PATH/HERE/ICC.m'); % Add a path to where the ICC function is stored. (Arash Salarian (2025). Intraclass Correlation Coefficient (ICC) (https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc), MATLAB Central File Exchange.)
addpath('YOUR/PATH/HERE/daboxplot/'); % Add a path to where the daboxplot function is stored. (Povilas Karvelis (2025). daboxplot (https://github.com/frank-pk/DataViz/releases/tag/v3.2.3), GitHub.)

%% Initialize code

dataset_name = 'YOUR_DATASET_NAME_HERE'; % The name of the dataset which will be used to name the output Excel file.

num_subj = 16; % Input the number of subjects.

% Brain Target
stn_subj = [1 6 7 12 13 14 16]; % Input the STN patient CONN ID order
gpi_subj = [2 3 4 5 8 9 10 11 15]; %Input the GPi patient CONN IDs order

% Input MDS-UPDRS Scores for % Motor improvement, raw rigidity score off DBS, raw
% tremor score off DBS, and raw bradykinesia score off DBS

mdstotal=[39.02 28.30 55.88 48.28 57.14 48.48 38.98 33.33 35.85 23.53 39.34 31.25 26.19 37.04 25.49 42.86];
mdsrigidity_rawoff=[11 11 3	8 11 6 8 13	12 10 11 10 6 3 13 8];
mdstremor_rawoff=[4	4 5	6 6	6 19 14	5 10 13 4 8 3 6 15];
mdsbrady_rawoff=[22	33 22 14 17	17 27 31 29	25 31 29 25 18 27 28]; 

% Add in a file containing framewise displacement values.
FD = 'FD_testretest_motion_values.xlsx';

FD_OFF_test = readmatrix(FD, 'Range','INSERT_RANGE_HERE'); % FD values for all patients - OFF test
FD_OFF_retest = readmatrix(FD, 'Range','INSERT_RANGE_HERE'); % FD values for all patients - OFF retest

ICC_type = 'A-1'; % Specify the ICC type, alpha value and r0 values for the ICC function.
alpha = 0.05; 
r0 = 0; 

categories = {'Whole Brain', 'Associative', 'Limbic', 'Motor'}; % The networks used in the analysis.

%% -------------------------------------------------------------------------- %%

               %PART 1: WITHIN SESSION ICC CALCULATION

%% -------------------------------------------------------------------------- %%

%% Load connectomes

cd(data_dir)

% Whole brain
load(results_whole); %load functional connectivity data.
corrname_whole = corrname; %load ROI names.
allrois_idx_whole = allrois_idx; %load ROI index.
brain_var_whole = brain_var; %load brain variability data.
fcs_whole = fcs; 

for i = 1:size(fcs_whole,1) %loop over each subject
for l = 1:size(fcs_whole,2) %loop over each condition (DBS OFFR1, DBS OFFR2, DBS ONR1, DBS ONR2)
fcs_whole(i,l,:,:) = atanh(fcs_whole(i,l,:,:)); %fisher transform functional connectivity values for a normal distribution.
end 
end
clear('allrois_idx','fcs','corrname', 'brain_var');

% Motor
load(results_motor); %load functional connectivity data.
corrname_motor = corrname; %load ROI names.
allrois_idx_motor = allrois_idx; %load ROI index.
brain_var_motor = brain_var; %load brain variability data.
fcs_motor = fcs;

for i = 1:size(fcs_motor,1) %loop over each subject
for l = 1:size(fcs_motor,2) %loop over each condition (DBS OFFR1, DBS OFFR2, DBS ONR1, DBS ONR2)
fcs_motor(i,l,:,:)=atanh(fcs_motor(i,l,:,:)); %fisher transform functional connectivity values for a normal distribution.
end 
end
clear('allrois_idx','fcs','corrname','brain_var');

% Limbic
load(results_limbic); %load functional connectivity data.
corrname_limbic = corrname; %load ROI names.
allrois_idx_limbic = allrois_idx; %load ROI index.
brain_var_limbic = brain_var;
fcs_limbic = fcs;

for i = 1:size(fcs_limbic,1)
for l = 1:size(fcs_limbic,2)
fcs_limbic(i,l,:,:) = atanh(fcs_limbic(i,l,:,:)); %fisher transform
end 
end
clear('allrois_idx','fcs','corrname','brain_var');

% Associative
load(results_associative); %load functional connectivity data.
corrname_associative = corrname; %load ROI names.
allrois_idx_associative = allrois_idx; %load ROI index.
brain_var_associative = brain_var;
fcs_associative=fcs;

for i = 1:size(fcs_associative,1)
for l = 1:size(fcs_associative,2)
fcs_associative(i,l,:,:) = atanh(fcs_associative(i,l,:,:)); %fisher transform
end 
end

clear('allrois_idx','fcs','corrname');


%% Prepare functional connectivity matrices as vectors for correlation plots

% Extract lower triangle from the matricies and remove center line of r=1

fcs_whole_tril = zeros(size(fcs_whole));
fcs_motor_tril = zeros(size(fcs_motor));
fcs_limbic_tril = zeros(size(fcs_limbic));
fcs_associative_tril = zeros(size(fcs_associative));

for i = 1:size(fcs_whole,1)
    for j = 1:size(fcs_whole,2)
fcs_whole_tril(i,j,:,:) = tril(squeeze(fcs_whole(i,j,:,:)),-1);
fcs_motor_tril(i,j,:,:) = tril(squeeze(fcs_motor(i,j,:,:)),-1);
fcs_limbic_tril(i,j,:,:) = tril(squeeze(fcs_limbic(i,j,:,:)),-1);
fcs_associative_tril(i,j,:,:) = tril(squeeze(fcs_associative(i,j,:,:)),-1);
    end
end

%% QC: Verify tril worked. If the plots show that only the lower triangular portion (excluding the diagonal) has non-zero colors and the upper triangular part is all zeros (a single uniform color representing zero), then the tril function has worked correctly.
figure;
% Plot the 'motor' plot
subplot(2,2,1);
imagesc(squeeze(fcs_motor_tril(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Motor Network'); 

% Plot the 'whole' plot
subplot(2,2,2);
imagesc(squeeze(fcs_whole_tril(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Whole Brain');

% Plot the 'associative' plot
subplot(2,2,3);
imagesc(squeeze(fcs_associative_tril(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Associative Network'); 

% Plot the 'limbic' plot
subplot(2,2,4);
imagesc(squeeze(fcs_limbic_tril(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Limbic Network');


%% calculate ICC values for functional connectivity and brain variability.

% Initialize matrices to store ICC results for each connectome
% FC = functional connectivity
% BV = brain variability
% 7 = 7 ICC output variables
FC_ICC_DBSOFF_whole = zeros(num_subj, 7); % FC ICC results for 'whole' connectome, DBS OFF
FC_ICC_DBSON_whole = zeros(num_subj, 7);  % FC ICC results for 'whole' connectome, DBS ON
BV_ICC_DBSOFF_whole = zeros(num_subj, 7); % BV ICC results for 'whole' connectome, DBS OFF
BV_ICC_DBSON_whole = zeros(num_subj, 7);  % BV ICC results for 'whole' connectome, DBS ON


FC_ICC_DBSOFF_motor = zeros(num_subj, 7); % FC ICC results for 'motor' connectome, DBS OFF
FC_ICC_DBSON_motor = zeros(num_subj, 7);  % FC ICC results for 'motor' connectome, DBS ON
BV_ICC_DBSOFF_motor = zeros(num_subj, 7); % BV ICC results for 'motor' connectome, DBS OFF
BV_ICC_DBSON_motor = zeros(num_subj, 7);  % BV ICC results for 'motor' connectome, DBS ON

FC_ICC_DBSOFF_limbic = zeros(num_subj, 7); % FC ICC results for 'limbic' connectome, DBS OFF
FC_ICC_DBSON_limbic = zeros(num_subj, 7);  % FC ICC results for 'limbic' connectome, DBS ON
BV_ICC_DBSOFF_limbic = zeros(num_subj, 7); % BV ICC results for 'limbic' connectome, DBS OFF
BV_ICC_DBSON_limbic = zeros(num_subj, 7);  % BV ICC results for 'limbic' connectome, DBS ON

FC_ICC_DBSOFF_associative = zeros(num_subj, 7); % FC ICC results for 'associative' connectome, DBS OFF
FC_ICC_DBSON_associative = zeros(num_subj, 7);  % FC ICC results for 'associative' connectome, DBS ON
BV_ICC_DBSOFF_associative = zeros(num_subj, 7); % BV ICC results for 'associative' connectome, DBS OFF
BV_ICC_DBSON_associative = zeros(num_subj, 7);  % BV ICC results for 'associative' connectome, DBS ON


% Arrange rows (measures) and columns (rater or test-retest)  &
% perform ICC
for i = 1:num_subj
    % Functional connectivity - whole brain
    M_OFF = [nonzeros(fcs_whole_tril(i, 1, :, :)) nonzeros(fcs_whole_tril(i, 2, :, :))];
    disp(M_OFF)
    M_ON = [nonzeros(fcs_whole_tril(i, 3, :, :)) nonzeros(fcs_whole_tril(i, 4, :, :))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0); %calc ICC
    FC_ICC_DBSOFF_whole(i, :) = [r_off, LB, UB, F, df1, df2, p]; %assign to results vector
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    FC_ICC_DBSON_whole(i, :) = [r_on, LB, UB, F, df1, df2, p];
   
    % Brain variability - whole brain
    M_OFF = [nonzeros(brain_var_whole(i, 1, :)) nonzeros(brain_var_whole(i, 2, :))];
    disp(M_OFF)
    M_ON = [nonzeros(brain_var_whole(i, 3, :)) nonzeros(brain_var_whole(i, 4, :))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    BV_ICC_DBSOFF_whole(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    BV_ICC_DBSON_whole(i, :) = [r_on, LB, UB, F, df1, df2, p];
    
    % Functional connectivity - Motor areas
    M_OFF = [nonzeros(fcs_motor_tril(i, 1, :, :)) nonzeros(fcs_motor_tril(i, 2, :, :))];
    M_ON = [nonzeros(fcs_motor_tril(i, 3, :, :)) nonzeros(fcs_motor_tril(i, 4, :, :))];
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    FC_ICC_DBSOFF_motor(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    FC_ICC_DBSON_motor(i, :) = [r_on, LB, UB, F, df1, df2, p];
    
    % Brain variability - Motor areas
    M_OFF = [nonzeros(brain_var_motor(i, 1, :, :)) nonzeros(brain_var_motor(i, 2, :, :))];
    disp(M_OFF)
    M_ON = [nonzeros(brain_var_motor(i, 3, :, :)) nonzeros(brain_var_motor(i, 4, :, :))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    BV_ICC_DBSOFF_motor(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    BV_ICC_DBSON_motor(i, :) = [r_on, LB, UB, F, df1, df2, p];

    % Functional connectivity - Limbic areas
    M_OFF = [nonzeros(fcs_limbic_tril(i, 1, :, :)) nonzeros(fcs_limbic_tril(i, 2, :, :))];
    M_ON = [nonzeros(fcs_limbic_tril(i, 3, :, :)) nonzeros(fcs_limbic_tril(i, 4, :, :))];
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    FC_ICC_DBSOFF_limbic(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    FC_ICC_DBSON_limbic(i, :) = [r_on, LB, UB, F, df1, df2, p];
    
    % Brain variability - Limbic areas
    M_OFF = [nonzeros(brain_var_limbic(i, 1, :, :)) nonzeros(brain_var_limbic(i, 2, :, :))];
    disp(M_OFF)
    M_ON = [nonzeros(brain_var_limbic(i, 3, :, :)) nonzeros(brain_var_limbic(i, 4, :, :))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    BV_ICC_DBSOFF_limbic(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    BV_ICC_DBSON_limbic(i, :) = [r_on, LB, UB, F, df1, df2, p];

    % Functional connectivity - Associative areas
    M_OFF = [nonzeros(fcs_associative_tril(i, 1, :, :)) nonzeros(fcs_associative_tril(i, 2, :, :))];
    M_ON = [nonzeros(fcs_associative_tril(i, 3, :, :)) nonzeros(fcs_associative_tril(i, 4, :, :))];
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    FC_ICC_DBSOFF_associative(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    FC_ICC_DBSON_associative(i, :) = [r_on, LB, UB, F, df1, df2, p];

    % Brain variability - Associative areas
    M_OFF = [nonzeros(brain_var_associative(i, 1, :, :)) nonzeros(brain_var_associative(i, 2, :, :))];
    disp(M_OFF)
    M_ON = [nonzeros(brain_var_associative(i, 3, :, :)) nonzeros(brain_var_associative(i, 4, :, :))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    BV_ICC_DBSOFF_associative(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    BV_ICC_DBSON_associative(i, :) = [r_on, LB, UB, F, df1, df2, p];

end

%% Save ICC Values to Excel

cd(output_dir)

% Functional Connectivity ICC Values
% Convert ICC matrices to tables

% Define header for ICC table
fc_icc_header = {'FC ICC', 'Lower Bound', 'Upper Bound', 'F', 'df1', 'df2', 'p-value'};

% Create tables for each connectome type and DBS condition
table_fc_whole_dbs_off = array2table(FC_ICC_DBSOFF_whole, 'VariableNames', fc_icc_header);
table_fc_whole_dbs_on = array2table(FC_ICC_DBSON_whole, 'VariableNames', fc_icc_header);

table_fc_motor_dbs_off = array2table(FC_ICC_DBSOFF_motor, 'VariableNames', fc_icc_header);
table_fc_motor_dbs_on = array2table(FC_ICC_DBSON_motor, 'VariableNames', fc_icc_header);

table_fc_limbic_dbs_off = array2table(FC_ICC_DBSOFF_limbic, 'VariableNames', fc_icc_header);
table_fc_limbic_dbs_on = array2table(FC_ICC_DBSON_limbic, 'VariableNames', fc_icc_header);

table_fc_associative_dbs_off = array2table(FC_ICC_DBSOFF_associative, 'VariableNames', fc_icc_header);
table_fc_associative_dbs_on = array2table(FC_ICC_DBSON_associative, 'VariableNames', fc_icc_header);

% Write tables to Excel file

% Define the filename for the Excel file
filename = [dataset_name '_FC_ICC_Results.xlsx'];

% Write each table to a separate sheet in the Excel file
writetable(table_fc_whole_dbs_off, filename, 'Sheet', 'Whole_DBS_OFF');
writetable(table_fc_whole_dbs_on, filename, 'Sheet', 'Whole_DBS_ON');

writetable(table_fc_motor_dbs_off, filename, 'Sheet', 'Motor_DBS_OFF');
writetable(table_fc_motor_dbs_on, filename, 'Sheet', 'Motor_DBS_ON');

writetable(table_fc_limbic_dbs_off, filename, 'Sheet', 'Limbic_DBS_OFF');
writetable(table_fc_limbic_dbs_on, filename, 'Sheet', 'Limbic_DBS_ON');

writetable(table_fc_associative_dbs_off, filename, 'Sheet', 'Associative_DBS_OFF');
writetable(table_fc_associative_dbs_on, filename, 'Sheet', 'Associative_DBS_ON');

%Brain Variability
% Convert ICC matrices to tables

% Define header for ICC table
bv_icc_header = {'BV ICC', 'Lower Bound', 'Upper Bound', 'F', 'df1', 'df2', 'p-value'};

% Create tables for each connectome type and DBS condition
table_bv_whole_dbs_off = array2table(BV_ICC_DBSOFF_whole, 'VariableNames', bv_icc_header);
table_bv_whole_dbs_on = array2table(BV_ICC_DBSON_whole, 'VariableNames', bv_icc_header);

table_bv_motor_dbs_off = array2table(BV_ICC_DBSOFF_motor, 'VariableNames', bv_icc_header);
table_bv_motor_dbs_on = array2table(BV_ICC_DBSON_motor, 'VariableNames', bv_icc_header);

table_bv_limbic_dbs_off = array2table(BV_ICC_DBSOFF_limbic, 'VariableNames', bv_icc_header);
table_bv_limbic_dbs_on = array2table(BV_ICC_DBSON_limbic, 'VariableNames', bv_icc_header);

table_bv_associative_dbs_off = array2table(BV_ICC_DBSOFF_associative, 'VariableNames', bv_icc_header);
table_bv_associative_dbs_on = array2table(BV_ICC_DBSON_associative, 'VariableNames', bv_icc_header);

% Write tables to Excel file

% Define the filename for the Excel file
filename = [dataset_name '_BV_ICC_Results.xlsx'];

% Write each table to a separate sheet in the Excel file
writetable(table_bv_whole_dbs_off, filename, 'Sheet', 'Whole_DBS_OFF');
writetable(table_bv_whole_dbs_on, filename, 'Sheet', 'Whole_DBS_ON');

writetable(table_bv_motor_dbs_off, filename, 'Sheet', 'Motor_DBS_OFF');
writetable(table_bv_motor_dbs_on, filename, 'Sheet', 'Motor_DBS_ON');

writetable(table_bv_limbic_dbs_off, filename, 'Sheet', 'Limbic_DBS_OFF');
writetable(table_bv_limbic_dbs_on, filename, 'Sheet', 'Limbic_DBS_ON');

writetable(table_bv_associative_dbs_off, filename, 'Sheet', 'Associative_DBS_OFF');
writetable(table_bv_associative_dbs_on, filename, 'Sheet', 'Associative_DBS_ON');


%% -------------------------------------------------------------------------- %%

               %PART 2: ANALYSIS (FIGURES 2, 3 and 6)

%% -------------------------------------------------------------------------- %%

%% Extract ICC values for plotting
FC_ICC_DBSON_r_whole = FC_ICC_DBSON_whole(:,1);
FC_ICC_DBSOFF_r_whole = FC_ICC_DBSOFF_whole(:,1);
BV_ICC_DBSON_r_whole = BV_ICC_DBSON_whole(:,1);
BV_ICC_DBSOFF_r_whole = BV_ICC_DBSOFF_whole(:,1);

FC_ICC_DBSON_r_motor = FC_ICC_DBSON_motor(:,1);
FC_ICC_DBSOFF_r_motor = FC_ICC_DBSOFF_motor(:,1);
BV_ICC_DBSON_r_motor = BV_ICC_DBSON_motor(:,1);
BV_ICC_DBSOFF_r_motor = BV_ICC_DBSOFF_motor(:,1);

FC_ICC_DBSON_r_limbic = FC_ICC_DBSON_limbic(:,1);
FC_ICC_DBSOFF_r_limbic = FC_ICC_DBSOFF_limbic(:,1);
BV_ICC_DBSON_r_limbic = BV_ICC_DBSON_limbic(:,1);
BV_ICC_DBSOFF_r_limbic = BV_ICC_DBSOFF_limbic(:,1);

FC_ICC_DBSON_r_associative = FC_ICC_DBSON_associative(:,1);
FC_ICC_DBSOFF_r_associative = FC_ICC_DBSOFF_associative(:,1);
BV_ICC_DBSON_r_associative = BV_ICC_DBSON_associative(:,1);
BV_ICC_DBSOFF_r_associative = BV_ICC_DBSOFF_associative(:,1);

% Create a string that labels each subject as STN or GPi.
grp = strings(1,num_subj);
grp(stn_subj)='STN';
grp(gpi_subj)='GPi';

% Define categories and corresponding ICC data
FC_ICC_DBSOFF = {FC_ICC_DBSOFF_r_whole, FC_ICC_DBSOFF_r_associative, FC_ICC_DBSOFF_r_limbic, FC_ICC_DBSOFF_r_motor};
FC_ICC_DBSON = {FC_ICC_DBSON_r_whole, FC_ICC_DBSON_r_associative, FC_ICC_DBSON_r_limbic, FC_ICC_DBSON_r_motor};
BV_ICC_DBSOFF = {BV_ICC_DBSOFF_r_whole, BV_ICC_DBSOFF_r_associative, BV_ICC_DBSOFF_r_limbic, BV_ICC_DBSOFF_r_motor};
BV_ICC_DBSON = {BV_ICC_DBSON_r_whole, BV_ICC_DBSON_r_associative, BV_ICC_DBSON_r_limbic, BV_ICC_DBSON_r_motor};

for k = 1:length(categories)
    % Extract data for current category.
    DBS_OFF_data = FC_ICC_DBSOFF{k}; 
    DBS_ON_data = FC_ICC_DBSON{k};    

    % Combine DBS OFF and DBS ON data for the current category.
    FC_combined_ICC = [DBS_OFF_data, DBS_ON_data];  % Combine as a numeric array

end


for k = 1:length(categories)
    % Extract data for current category
    DBS_OFF_data = BV_ICC_DBSOFF{k};  
    DBS_ON_data = BV_ICC_DBSON{k};    

    % Combine DBS OFF and DBS ON data for the current category
    BV_combined_ICC = [DBS_OFF_data, DBS_ON_data];  % Combine as a numeric array

end

%% Reshape data for box plots of ICC Values for all patients and networks (grouped)

M_FC_ICC_dbsoff=zeros(num_subj,4,1); % create empty reshaped matrices for functional connectivity ICC values.
M_FC_ICC_dbson=zeros(num_subj,4,1);
M_BV_ICC_dbsoff=zeros(num_subj,4,1); % create empty reshaped matrices for brain variability ICC values.
M_BV_ICC_dbson=zeros(num_subj,4,1);

for r = 1:4 %categories = {'Whole Brain', 'Associative', 'Limbic', 'Motor'};
    M_FC_ICC_dbsoff(:, r, 1) = FC_ICC_DBSOFF{r};
    M_FC_ICC_dbson(:, r, 1) = FC_ICC_DBSON{r};
    M_BV_ICC_dbsoff(:, r, 1) = BV_ICC_DBSOFF{r};
    M_BV_ICC_dbson(:, r, 1) = BV_ICC_DBSON{r};
end

% stack DBS ON/DBS OFF for grouping

FC_ICC_stacked = [M_FC_ICC_dbson ; M_FC_ICC_dbsoff];
BV_ICC_stacked = [M_BV_ICC_dbson ; M_BV_ICC_dbsoff];

group_names = {'ON', 'OFF'};
group_inx = [ones(1,num_subj), 2.*ones(1,num_subj)]; % This creates an index vector indicating the group for each subject. 1 refers to DBS ON and 2 refers to DBS OFF.

%% Figure 2A - Functional connectivity ICC box plot

figure('Color', 'w');
group_colors = [[0.5, 0.8, 1]; [0, 0.45, 0.7]]; % Set colors for the DBS conditions
daboxplot(FC_ICC_stacked,'groups',group_inx,'color', group_colors,'whiskers' , 0, 'scatter', 1, 'scattersize', 80, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories,'fill',1,'legend',group_names,'boxspacing',0.6, 'boxwidth', 1.5);
legend('Location', 'southeast');
ylabel('ICC', 'FontSize',25);
ylim([0.3,0.9]);
yticks(0.3:0.2:0.9);
xlim([0.5, 4.5]);
set(gca, 'XTickLabel', categories, 'FontSize', 20, 'XTickLabelRotation', 0);
for j=1:4 %categories = {'whole', 'associative', 'limbic', 'motor'};
[p, h] = signrank(squeeze(M_FC_ICC_dbsoff(:,j)), squeeze(M_FC_ICC_dbson(:,j)));
    str=sprintf('p= %1.2f',p);
    x_pos = j;
    T = text(x_pos, max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'center');
effect = meanEffectSize(squeeze(M_FC_ICC_dbsoff(:,j)), squeeze(M_FC_ICC_dbson(:,j)),Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(effect);
end

%% Figure 2B - Functional Connectivity ICC Scatter Plot
% Define the number of rows and columns and background color for the FC subplots
rows = 2;
cols = 2;
figure('Color', 'w'); 

% Initialize an empty array to hold the scatter plot handles for the legend
hLegend = [];

% Plot FC ICC values
for k = 1:length(categories)

    % Create a new subplot in a 2x2 grid
    subplot(rows, cols, k); % Specify the grid location for the current plot

    % Separate data by group
    isSTN = strcmp(grp, 'STN'); % Logical index for STN group
    isGPi = strcmp(grp, 'GPi'); % Logical index for GPi group

    % Plot STN group
    hSTN = scatter(FC_ICC_DBSOFF{k}(isSTN), FC_ICC_DBSON{k}(isSTN), 200, 'o', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
    hold on;

    % Plot GPi group
    hGPi = scatter(FC_ICC_DBSOFF{k}(isGPi), FC_ICC_DBSON{k}(isGPi), 200, 'v', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);

    % Customize plot appearance
    set(gca, 'FontSize', 20);
    xlabel('ICC: OFF');
    ylabel('ICC: ON');
    title(categories{k}, 'FontWeight', 'bold', 'fontsize', 20);

    % Set the same axis limits for all subplots
    xlim([0.3 0.9]);
    xticks(0.3:0.2:0.9);
    ylim([0.3 0.9]);
    yticks(0.3:0.2:0.9);
    ax = gca; % Get current axes handle
    ax.Box = 'off'; % Turn off the box

    % Calculate and display the correlation coefficient.
    [r,p,rl,ru] = corrcoef(FC_ICC_DBSOFF{k}, FC_ICC_DBSON{k});
    str_r = sprintf('r = %1.2f', r(1,2)); % Format the correlation coefficient as a string
    r_print = text(0.33, 0.89, str_r); % Display the correlation on the plot
    set(r_print, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    str_p = sprintf('p = %1.2f', p(1,2));
    p_print = text(0.33, 0.83, str_p); 
    set(p_print, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    disp(['Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
    
    % Add a linear fit line to the plot.
    P = polyfit(FC_ICC_DBSOFF{k}, FC_ICC_DBSON{k}, 1); % Perform linear fit
    Bfit = polyval(P, FC_ICC_DBSOFF{k}); % Evaluate the fit
    h_line_best_fit = plot(FC_ICC_DBSOFF{k}, Bfit, '-k', 'LineWidth', 1.5); % Plot the fit line in black

    % Ensure that the plot is held for multiple annotations
    hold off;

    % Store handles for legend
    hLegend = [hLegend; hSTN; hGPi];
end

% Create a common legend for all subplots
legend([hSTN hGPi], {'STN', 'GPi'}, 'Location', 'southwest', 'FontSize', 15);

%% Figure 2C: Motor Network Connectomes of Patients with the Highest and Lowest Motor ICC Values

%Identify patients with lowest and highest ICC for Motor Network: DBS ON

cd(output_dir)

motor_on_ICC = readtable([dataset_name "_FC_ICC_Results.xlsx"], 'Sheet', 'Motor_DBS_ON');

% Find minimum and maximum values
[min_ICC, minIdx] = min(motor_on_ICC{:,"FCICC"});
[max_ICC, maxIdx] = max(motor_on_ICC{:,"FCICC"});

% Reorder ROIs so that they always appear L then R

new_order = [1:8, 10,9,12,11,14,13,16,15,30,29,32,31,34,33,18,17,20,19,22,21,24,23,26,25,28,27,35:60];
corrname_motor_reordered = corrname_motor(new_order);
disp(corrname_motor_reordered);

% Remove L and R from the ROI names so the labels fit in figure

match=" " + wildcardPattern + "r";
match_idx = contains(corrname_motor_reordered, match);
corrname_motor_reordered(match_idx)='';
match=" " + wildcardPattern + "l";
 corrname_motor_reordered = erase(corrname_motor_reordered, match);
display(corrname_motor_reordered); 

%plot pt with minimum ICC value
for i=3:4 % 3, 4 = on test, on retest
    figure('Position',[0 0 900 750]);
    set(gcf,'color','w');
    fc=squeeze(fcs_motor(minIdx,:,:,:)); %choose pt with lowest ICC
    imagesc(squeeze(fc(i,new_order,new_order))); %plots only one run at a time
    set(gca,'XTickLabel',corrname_motor_reordered,'XTickLabelRotation', 90);
    set(gca,'YTickLabel',corrname_motor_reordered);
    set(gca, 'XTick', (1:length(corrname_motor_reordered)), 'FontSize', 15);
    set(gca, 'YTick', (1:length(corrname_motor_reordered)), 'FontSize', 15);
    grid off
    box off
    colormap("parula");
    if i==1
    title('OFF (Test)' ,'FontSize', 25);
    elseif i==2
    title('OFF (Retest)' ,'FontSize', 25);
    elseif i==3
    title('ON (Test)' ,'FontSize', 25);
    elseif i==4
    title('ON (Retest)' ,'FontSize', 25);
    end
    a = colorbar; clim([-1, 1]);
    ylabel(a,'Functional Connectivity','FontSize',25,'Rotation',90);

end

%plot difference between test and retest for minimum ICC value
fc=squeeze(fcs_motor(minIdx,:,:,:)); %choose pt with lowest ICC
difference = abs(squeeze(fc(4,new_order,new_order))-squeeze(fc(3, new_order, new_order))); 
difference(isnan(difference) | isinf(difference)) = 0;
figure('Position',[0 0 900 750]);
set(gcf,'color','w');
imagesc(difference); 
set(gca,'XTickLabel',corrname_motor_reordered,'XTickLabelRotation', 90);
set(gca,'YTickLabel',corrname_motor_reordered);
set(gca, 'XTick', (1:length(corrname_motor_reordered)), 'FontSize', 15);
set(gca, 'YTick', (1:length(corrname_motor_reordered)), 'FontSize', 15);
grid off
box off
colormap("parula");
title('ON (Test-Retest)' ,'FontSize', 25);

a = colorbar; clim([-1, 1]);
ylabel(a,'Functional Connectivity','FontSize',25,'Rotation',90);


%plot pt with max ICC value
for i=3:4
    figure('Position',[0 0 900 750]);
    set(gcf,'color','w');
    fc=squeeze(fcs_motor(maxIdx,:,:,:)); %choose pt with lowest ICC
    imagesc(squeeze(fc(i,new_order,new_order))); %plots only one run at a time
    set(gca,'XTickLabel',corrname_motor_reordered,'XTickLabelRotation', 90);
    set(gca,'YTickLabel',corrname_motor_reordered);
    set(gca, 'XTick', 1:length(corrname_motor_reordered), 'FontSize', 15);
    set(gca, 'YTick', 1:length(corrname_motor_reordered), 'FontSize', 15);
    grid off
    box off
    colormap("parula");
    if i==1
    title('OFF (Test)' ,'FontSize', 25);
    elseif i==2
    title('OFF (Retest)' ,'FontSize', 25);
    elseif i==3
    title('ON (Test)' ,'FontSize', 25);
    elseif i==4
    title('ON (Retest)' ,'FontSize', 25);
    end
    a = colorbar; caxis([-1, 1]);
    ylabel(a,'Functional Connectivity','FontSize',25,'Rotation',90);

end

%plot difference between test and retest for max ICC value
fc=squeeze(fcs_motor(maxIdx,:,:,:)); %choose pt with highest ICC
difference = abs(squeeze(fc(4,new_order,new_order))-squeeze(fc(3, new_order, new_order))); 
difference(isnan(difference) | isinf(difference)) = 0;
figure('Position',[0 0 900 750]);
set(gcf,'color','w');
imagesc(difference); 
set(gca,'XTickLabel',corrname_motor_reordered,'XTickLabelRotation', 90);
set(gca,'YTickLabel',corrname_motor_reordered);
set(gca, 'XTick', (1:length(corrname_motor_reordered)), 'FontSize', 15);
set(gca, 'YTick', (1:length(corrname_motor_reordered)), 'FontSize', 15);
grid off
box off
colormap("parula");
title('ON (Test-Retest)' ,'FontSize', 25);

a = colorbar; clim([-1, 1]);
ylabel(a,'Functional Connectivity','FontSize',25,'Rotation',90);

% Make scatterplot of each patient's FC values test vs. retest.
fc_max_on_test=squeeze(fcs_motor(maxIdx,3,:,:)); 
fc_max_on_test=fc_max_on_test(:); % flatten to a 1D vector for plotting.
fc_max_on_retest=squeeze(fcs_motor(maxIdx,4,:,:)); 
fc_max_on_retest=fc_max_on_retest(:); % flatten to a 1D vector for plotting.

fc_min_on_test=squeeze(fcs_motor(minIdx,3,:,:)); 
fc_min_on_test=fc_min_on_test(:); % flatten to a 1D vector for plotting.
fc_min_on_retest=squeeze(fcs_motor(minIdx,4,:,:)); 
fc_min_on_retest=fc_min_on_retest(:); % flatten to a 1D vector for plotting.

% Remove NaN and Inf values
% Find finite, real values within the range [-1, 1] in both datasets
valid_idx_max = isfinite(fc_max_on_test) & isfinite(fc_max_on_retest) & ...
            ~isnan(fc_max_on_test) & ~isnan(fc_max_on_retest) & ...
            fc_max_on_test >= -1 & fc_max_on_test <= 1 & ...
            fc_max_on_retest >= -1 & fc_max_on_retest <= 1;

valid_idx_min = isfinite(fc_min_on_test) & isfinite(fc_min_on_retest) & ...
            ~isnan(fc_min_on_test) & ~isnan(fc_min_on_retest) & ...
            fc_min_on_test >= -1 & fc_min_on_test <= 1 & ...
            fc_min_on_retest >= -1 & fc_min_on_retest <= 1;

% Filter out the corresponding NaN/Inf values
fc_max_on_test_clean = fc_max_on_test(valid_idx_max);
fc_max_on_retest_clean = fc_max_on_retest(valid_idx_max);

fc_min_on_test_clean = fc_min_on_test(valid_idx_min);
fc_min_on_retest_clean = fc_min_on_retest(valid_idx_min);

% max ICC scatterplot
figure('Position', [100, 100, 400, 400]);
scatter(fc_max_on_test_clean, fc_max_on_retest_clean, 40, 'o', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1, 'LineWidth', 1.5);
xlabel('FC: ON (Test)', 'FontSize', 20);
ylabel('FC: ON (Retest)', 'FontSize', 20);
set(gca, 'FontSize', 20);

title('Max ICC', 'FontSize', 25);
[r,p] = corr(fc_max_on_test_clean, fc_max_on_retest_clean);
str_r = sprintf('r = %1.2f', r); % Format the correlation coefficient as a string
r_print = text(-0.9, 0.8, str_r, "FontSize", 20); % Display the correlation on the plot
str_p = sprintf('p < 0.001'); % Format the correlation coefficient as a string
p_print = text(-0.9, 0.6, str_p, "FontSize", 20); % Display the correlation on the plot

hold on

P = polyfit(fc_max_on_test_clean, fc_max_on_retest_clean, 1); % Perform linear fit
Bfit = polyval(P, fc_max_on_test_clean); % Evaluate the fit
h_line_best_fit = plot(fc_max_on_test_clean, Bfit, '-k', 'LineWidth', 1.5); % Plot the fit line in black

 % Set axes limits
xlim([-1, 1]);
ylim([-1, 1]);

xticks([-1, 0, 1]);
yticks([-1,0, 1]);

% Min ICC scatterplot
figure('Position', [100, 100, 400, 400]);
scatter(fc_min_on_test_clean, fc_min_on_retest_clean,40, 'o', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1, 'LineWidth', 1.5);
xlabel('FC: ON (Test)', 'FontSize', 20);
ylabel('FC: ON (Retest)', 'FontSize', 20);
set(gca, 'FontSize', 20);

title('Min ICC', 'FontSize', 25);
[r,p] = corr(fc_min_on_test_clean, fc_min_on_retest_clean);
str_r = sprintf('r = %1.2f', r); % Format the correlation coefficient as a string
r_print = text(-0.9, 0.8, str_r,"FontSize", 20); % Display the correlation on the plot
str_p = sprintf('p < 0.001'); % Format the correlation coefficient as a string
p_print = text(-0.9, 0.6, str_p,"FontSize", 20); % Display the correlation on the plot

hold on

P = polyfit(fc_min_on_test_clean, fc_min_on_retest_clean, 1); % Perform linear fit
Bfit = polyval(P, fc_min_on_test_clean); % Evaluate the fit
h_line_best_fit = plot(fc_min_on_test_clean, Bfit, '-k', 'LineWidth', 1.5); % Plot the fit line in black

 % Set axes limits
xlim([-1, 1]);
ylim([-1, 1]);

xticks([-1, 0, 1]);
yticks([-1, 0, 1]);

%% Figure 3A - Brain variability ICC Box plot

figure('Color', 'w'); 
group_colors = [[0.5, 0.8, 1]; [0, 0.45, 0.7]]; % Set colors for the DBS conditions
daboxplot(BV_ICC_stacked,'groups',group_inx,'color', group_colors,'whiskers' , 0, 'scatter', 1, 'scattersize', 80, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories,'fill',1,'legend',group_names,'boxspacing',0.6, 'boxwidth', 1.5);
legend('Location', 'southeast');
set(gca,'FontSize',15);
ylabel('ICC');
ylim([0.6,1]);
yticks(0.6:0.2:1);
xlim([0.5, 4.5]);
set(gca, 'XTickLabel', categories, 'FontSize', 20, 'XTickLabelRotation', 0);
for j=1:4 %categories = {'whole', 'associative', 'limbic', 'motor'};
[p, h] = signrank(squeeze(M_BV_ICC_dbsoff(:,j)), squeeze(M_BV_ICC_dbson(:,j)));
    str=sprintf('p= %1.2f',p);
    x_pos = j;
    T = text(x_pos, max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 15, 'verticalalignment', 'top', 'horizontalalignment', 'center');
effect = meanEffectSize(squeeze(M_BV_ICC_dbsoff(:,j)), squeeze(M_BV_ICC_dbson(:,j)),Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(effect);

end

%% Figure 3B - Brain Variability ICC Scatter Plot
% Define the number of rows and columns and background color for the FC subplots
rows = 2;
cols = 2;
figure('Color', 'w'); 

% Initialize an empty array to hold the scatter plot handles for the legend
hLegend = [];

% Plot FC ICC values
for k = 1:length(categories)

    % Create a new subplot in a 2x2 grid
    subplot(rows, cols, k); % Specify the grid location for the current plot

    % Separate data by group
    isSTN = strcmp(grp, 'STN'); % Logical index for STN group
    isGPi = strcmp(grp, 'GPi'); % Logical index for GPi group

    % Plot STN group
    hSTN = scatter(BV_ICC_DBSOFF{k}(isSTN), BV_ICC_DBSON{k}(isSTN), 200, 'o', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
    hold on;

    % Plot GPi group
    hGPi = scatter(BV_ICC_DBSOFF{k}(isGPi), BV_ICC_DBSON{k}(isGPi), 200, 'v', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);

    % Customize plot appearance
    set(gca, 'FontSize', 20);
    xlabel('ICC: OFF');
    ylabel('ICC: ON');
    title(categories{k}, 'FontWeight', 'bold', 'fontsize', 20);

    % Set the same axis limits for all subplots
    xlim([0.6 1]);
    xticks(0.6:0.2:1);
    ylim([0.6 1]);
    yticks(0.6:0.2:1);
    ax = gca; % Get current axes handle
    ax.Box = 'off'; % Turn off the box

    % Calculate and display the correlation coefficient.
    [r,p] = corr(BV_ICC_DBSOFF{k}, BV_ICC_DBSON{k});
    str_r = sprintf('r = %1.2f', r); % Format the correlation coefficient as a string
    r_print = text(0.61, 0.99, str_r); % Display the correlation on the plot
    set(r_print, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    if p < 0.01
        str_p = sprintf('p = %.2e', p);
    else
        str_p = sprintf('p = %1.2f', p);
    end
    p_print = text(0.61, 0.94, str_p); 
    set(p_print, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    [r,p,rl, ru] = corrcoef(BV_ICC_DBSOFF{k}, BV_ICC_DBSON{k});

    disp(['Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);

    % Add a linear fit line to the plot.
    P = polyfit(BV_ICC_DBSOFF{k}, BV_ICC_DBSON{k}, 1); % Perform linear fit
    Bfit = polyval(P, BV_ICC_DBSOFF{k}); % Evaluate the fit
    h_line_best_fit = plot(BV_ICC_DBSOFF{k}, Bfit, '-k', 'LineWidth', 1.5); % Plot the fit line in black
   
    hold off;

    % Store handles for legend
    hLegend = [hLegend; hSTN; hGPi];
end

% Create a common legend for all subplots
legend([hSTN hGPi], {'STN', 'GPi'}, 'Location', 'southwest', 'FontSize', 15);


%% Figure 3C: Heatmaps of whole brain brain variability values

% Define the original order for ROIs 109-142
roi_range_to_swap = 109:142;

% Create a new order for 109-142 by swapping every other ROI
roi_swapped = roi_range_to_swap;
roi_swapped(1:2:end) = roi_range_to_swap(2:2:end); % Swap odd indices with even
roi_swapped(2:2:end) = roi_range_to_swap(1:2:end); % Swap even indices with odd

% Extract the swapped subranges from roi_swapped
roi_109_118 = roi_swapped(1:10);
roi_119_120 = roi_swapped(11:12); % ROIs 119 and 120 (caudate)
roi_137_142 = roi_swapped(29:34); % ROIs 137 to 142 (GPi/GPe/STN)
roi_121_136 = roi_swapped(13:28); % ROIs 121 to 136 (TL/amyg/hippo)

% Define the final reordered ROIs - move dentate AFTER vermis, Gpi/STN before TL
reordered_rois = [1:108, roi_119_120, roi_109_118, roi_137_142, roi_121_136]; 

% Apply the new order to brain_var_whole
brain_var_whole_reordered = brain_var_whole(:, :, reordered_rois);

% Reorder corrname_whole to match the new ROI order
corrname_whole_reordered = corrname_whole(reordered_rois);

% Extract the conditions
offtest_bv_whole = squeeze(brain_var_whole_reordered(:, 1, :)); % Off Test
offretest_bv_whole = squeeze(brain_var_whole_reordered(:, 2, :)); % Off Retest
ontest_bv_whole = squeeze(brain_var_whole_reordered(:, 3, :)); % On Test
onretest_bv_whole = squeeze(brain_var_whole_reordered(:, 4, :)); % On Retest

% Create a figure
figure;

% Adjust margins and spacing
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot the heatmaps
% OFF (Test)
nexttile;
imagesc(offtest_bv_whole);
clim([0 4]);
set(gca, 'XTick', [],'FontSize', 18); % Remove x-axis tick marks for this plot

% OFF (Retest)
nexttile;
imagesc(offretest_bv_whole);
clim([0 4]);
set(gca, 'XTick', [],'FontSize', 18); % Remove x-axis tick marks for this plot

% ON (Test)
nexttile;
imagesc(ontest_bv_whole);
clim([0 4]);
set(gca, 'XTick', [],'FontSize', 18); % Remove x-axis tick marks for this plot

% ON (Retest)
nexttile;
imagesc(onretest_bv_whole);
clim([0 4]);
set(gca, 'XTick', corrname_whole_reordered,'FontSize', 18);

% Set colormap for all
colormap('parula'); 

% Add a single colorbar
cb = colorbar;
cb.Layout.Tile = 'east'; % Place it on the right side of all subplots

% Add a single x-axis label
xlabel('Whole Brain ROIs', 'FontSize', 20, 'FontWeight', 'bold');

% Add a shared y-axis label
% Create a dummy axis spanning the whole figure
ax = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
ylabel(ax, 'Patients', 'FontSize', 20, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');

%% Figure 6A - Scatter Plots of Sources of between-subject variance in DBS OFF ICC 
target_stn=ones(1,16); 
target_stn(gpi_subj)=0; %binary; 1=stn, 2=gpi

% Initialize arrays to store averages
FD_OFF_test_95 = zeros(1,16); % 95th percentile
FD_OFF_retest_95 = zeros(1,16); % 95th percentile


% Loop through each patient (column)
for i = 1:16
    % Extract the data for the current patient
    patient_test = FD_OFF_test(:, i);
    patient_retest = FD_OFF_retest(:,i);
    
    % Calculate the 75th and 95th percentiles
    FD_OFF_test_95(i) = prctile(patient_test, 95);
    FD_OFF_retest_95(i) = prctile(patient_retest, 95);

end

FD_dif_95 = abs(FD_OFF_test_95-FD_OFF_retest_95);

% Set up a 2x5 grid of subplots
figure('Color', 'w', 'Position', [0 0 900 500]);

% Plot 1: FC ICC OFF vs FD difference (95%)
subplot(2, 3, 1);
scatter(FD_dif_95, FC_ICC_DBSOFF_r_whole', 300, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
ylabel('FC ICC: OFF', 'FontSize', 25);
xlim([0 0.4]); xticks(0:0.2:0.4);
ylim([0.3 0.9]); yticks(0.3:0.2:0.9);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(FD_dif_95', FC_ICC_DBSOFF_r_whole);
text(0.05, 0.9, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(0.05, 0.82, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(FD_dif_95', FC_ICC_DBSOFF_r_whole);
disp(['6A Plot 1 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(FD_dif_95, FC_ICC_DBSOFF_r_whole, 1);
plot(FD_dif_95, polyval(P, FD_dif_95), '-k', 'LineWidth', 1.5);
hold off;

% Plot 2: BV ICC OFF vs FD difference (95%)
subplot(2, 3, 4);
scatter(FD_dif_95, BV_ICC_DBSOFF_r_whole', 300, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
xlabel('FD test – FD retest (95%)', 'FontSize', 25);
ylabel('BV ICC: OFF', 'FontSize', 25);
xlim([0 0.4]); xticks(0:0.2:0.4);
ylim([0.6 1]); yticks(0.6:0.2:1);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(FD_dif_95', BV_ICC_DBSOFF_r_whole);
text(0.25, 0.75, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(0.25, 0.7, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(FD_dif_95', BV_ICC_DBSOFF_r_whole);
disp(['6A Plot 2 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(FD_dif_95, BV_ICC_DBSOFF_r_whole, 1);
plot(FD_dif_95, polyval(P, FD_dif_95), '-k', 'LineWidth', 1.5);
hold off;

% Plot 3: FC ICC OFF vs Tremor raw off
subplot(2, 3, 2);
scatter(mdstremor_rawoff', FC_ICC_DBSOFF_r_whole, 300, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
%xlim([-1 3]); xticks(-1:2:3);
ylim([0.3 0.9]); yticks(0.3:0.2:0.9);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(mdstremor_rawoff', FC_ICC_DBSOFF_r_whole);
text(2, 0.9, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(2, 0.82, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(mdstremor_rawoff', FC_ICC_DBSOFF_r_whole);
disp(['6A Plot 3 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(mdstremor_rawoff, FC_ICC_DBSOFF_r_whole, 1);
plot(mdstremor_rawoff, polyval(P, mdstremor_rawoff), '-k', 'LineWidth', 1.5);
hold off;

% Plot 4: BV ICC OFF vs Tremor
subplot(2, 3, 5);
scatter(mdstremor_rawoff', BV_ICC_DBSOFF_r_whole, 250, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
xlabel('Raw Tremor Score', 'FontSize', 25);
%xlim([-1 3]); xticks(-1:2:3);
ylim([0.6 1]); yticks(0.6:0.2:1);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(mdstremor_rawoff', BV_ICC_DBSOFF_r_whole);
text(12, 0.75, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(12, 0.7, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(mdstremor_rawoff', BV_ICC_DBSOFF_r_whole);
disp(['6A Plot 4 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(mdstremor_rawoff, BV_ICC_DBSOFF_r_whole, 1);
plot(mdstremor_rawoff, polyval(P, mdstremor_rawoff), '-k', 'LineWidth', 1.5);
hold off;

% Plot 5: FC ICC OFF vs Rigidity + Bradykinesia
subplot(2, 3, 3);
scatter(mdsrigidity_rawoff' + mdsbrady_rawoff', FC_ICC_DBSOFF_r_whole, 250, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
%xlim([-3 3]); xticks(-3:2:3);
ylim([0.3 0.9]); yticks(0.3:0.2:0.9);
ax = gca; ax.Box = 'off'; 
ax.LineWidth = 3;
[r, p] = corr(mdsrigidity_rawoff' + mdsbrady_rawoff', FC_ICC_DBSOFF_r_whole);
text(22, 0.9, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(22, 0.82, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(mdsrigidity_rawoff' + mdsbrady_rawoff', FC_ICC_DBSOFF_r_whole);
disp(['6A Plot 5 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(mdsrigidity_rawoff' + mdsbrady_rawoff', FC_ICC_DBSOFF_r_whole, 1);
plot(mdsrigidity_rawoff' + mdsbrady_rawoff', polyval(P, mdsrigidity_rawoff' + mdsbrady_rawoff'), '-k', 'LineWidth', 1.5);
hold off;

% Plot 6: BV ICC OFF vs Rigidity + Bradykinesia
subplot(2, 3, 6);
scatter(mdsrigidity_rawoff' + mdsbrady_rawoff', BV_ICC_DBSOFF_r_whole, 250, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
xlabel('Raw Rigidity + Brady Score', 'FontSize', 25);
%xlim([-3 3]); xticks(-3:2:3);
ylim([0.6 1]); yticks(0.6:0.2:1);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(mdsrigidity_rawoff' + mdsbrady_rawoff', BV_ICC_DBSOFF_r_whole);
text(36, 0.75, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(36, 0.7, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(mdsrigidity_rawoff' + mdsbrady_rawoff', BV_ICC_DBSOFF_r_whole);
disp(['6A Plot 6 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(mdsrigidity_rawoff' + mdsbrady_rawoff', BV_ICC_DBSOFF_r_whole, 1);
plot(mdsrigidity_rawoff' + mdsbrady_rawoff', polyval(P, mdsrigidity_rawoff' + mdsbrady_rawoff'), '-k', 'LineWidth', 1.5);
hold off;


%% Figure 6B - Scatter Plots of Sources of between-subject variance in DBS ON ICC 

% Plot 1: FC ICC ON vs % Total Motor Improvement
figure('Color', 'w', 'Position', [0 0 600 500]);
subplot(2, 2, 1);
scatter(mdstotal, FC_ICC_DBSON_r_whole, 300, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
ylabel('FC ICC: ON', 'FontSize', 25);
ylim([0.3 0.9]); yticks(0.3:0.2:0.9);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(mdstotal', FC_ICC_DBSON_r_whole);
text(22, 0.9, ['r = ', num2str(r, '%.2e')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(22, 0.82, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(mdstotal', FC_ICC_DBSON_r_whole);
disp(['6B Plot 1 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(mdstotal', FC_ICC_DBSON_r_whole, 1);
plot(mdstotal, polyval(P, mdstotal), '-k', 'LineWidth', 1.5);
hold off;

% Plot 2: BV ICC ON vs % Total Motor Improvement
subplot(2, 2, 3);
scatter(mdstotal, BV_ICC_DBSON_r_whole, 300, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
xlabel('% Total Motor Improvement', 'FontSize', 25);
ylabel('BV ICC: ON', 'FontSize', 25);
ylim([0.6 1]); yticks(0.6:0.2:1);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(mdstotal', BV_ICC_DBSON_r_whole);
text(45, 0.75, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(45, 0.7, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(mdstotal', BV_ICC_DBSON_r_whole);
disp(['6B Plot 2 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(mdstotal', BV_ICC_DBSON_r_whole, 1);
plot(mdstotal, polyval(P, mdstotal), '-k', 'LineWidth', 1.5);
hold off;

% Plot 3: FC ICC ON vs Total Brain Response

fcs_whole_tril_reshape = zeros([size(fcs_whole,1) size(fcs_whole,2) nnz(fcs_whole_tril(i,j,:,:))]);

for i=1:size(fcs_whole,1)
    for j=1:size(fcs_whole,2)
temp = squeeze(fcs_whole_tril(i,j,:,:));
fcs_whole_tril_reshape(i,j,:) = nonzeros(temp(:));
    end
end

mean_ontestretest = (squeeze(fcs_whole_tril_reshape(:,3,:))+squeeze(fcs_whole_tril_reshape(:,4,:)))/2; % create a new matrix containing mean FC values for whole brain for all 16 patients.
mag_response=range(mean_ontestretest,2); % Find the range of FC values for each patient.

subplot(2, 2, 2);
scatter(abs(mag_response) .* mdstotal', FC_ICC_DBSON_r_whole, 300, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
ylim([0.3 0.9]); yticks(0.3:0.2:0.9);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(abs(mag_response) .* mdstotal', FC_ICC_DBSON_r_whole);
text(55, 0.9, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(55, 0.82, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(abs(mag_response) .* mdstotal', FC_ICC_DBSON_r_whole);
disp(['6B Plot 3 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(abs(mag_response) .* mdstotal', FC_ICC_DBSON_r_whole, 1);
plot(abs(mag_response) .* mdstotal', polyval(P, abs(mag_response) .* mdstotal'), '-k', 'LineWidth', 1.5);
hold off;

% Plot 4: BV ICC ON vs Total Brain Response

subplot(2, 2, 4);
scatter(abs(mag_response) .* mdstotal', BV_ICC_DBSON_r_whole, 300, 'o', ...
    'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1.5);
set(gca, 'FontSize', 25);
xlabel('%Total Motor x Range FC', 'FontSize', 25);
ylim([0.6 1]); yticks(0.6:0.2:1);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
[r, p] = corr(abs(mag_response) .* mdstotal', BV_ICC_DBSON_r_whole);
text(95, 0.75, ['r = ', num2str(r, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
text(95, 0.7, ['p = ', num2str(p, '%.2f')], 'FontSize', 25, ...
    'verticalalignment', 'top', 'horizontalalignment', 'left');
hold on;
[r,p,rl,ru] = corrcoef(abs(mag_response) .* mdstotal', BV_ICC_DBSON_r_whole);
disp(['6B Plot 4 Statistical Values:''r=' r(1,2) 'p=' p(1,2) '95% CI LB=' rl(1,2) '95% CI UB=' ru(1,2)]);
P = polyfit(abs(mag_response) .* mdstotal', BV_ICC_DBSON_r_whole, 1);
plot(abs(mag_response) .* mdstotal', polyval(P, abs(mag_response) .* mdstotal'), '-k', 'LineWidth', 1.5);
hold off;

%% Figure 6C - Brain target(STN vs. GPi) vs. ICC box plot

% Reshape data for box plots of ICC Values for all patients (STN vs. GPi) and networks (grouped)

M_FC_ICC_dbson_stn = M_FC_ICC_dbson(stn_subj,:,:);
M_FC_ICC_dbson_gpi = M_FC_ICC_dbson(gpi_subj,:,:);

M_FC_ICC_dbsoff_stn = M_FC_ICC_dbsoff(stn_subj,:,:);
M_FC_ICC_dbsoff_gpi = M_FC_ICC_dbsoff(gpi_subj,:,:);

M_BV_ICC_dbson_stn = M_BV_ICC_dbson(stn_subj,:,:);
M_BV_ICC_dbson_gpi = M_BV_ICC_dbson(gpi_subj,:,:);

M_BV_ICC_dbsoff_stn = M_BV_ICC_dbsoff(stn_subj,:,:);
M_BV_ICC_dbsoff_gpi = M_BV_ICC_dbsoff(gpi_subj,:,:);

% Pad STN group with NaN to match size of GPi group for plotting 
padding = NaN(2,4);
M_FC_ICC_dbson_stn_padded = cat(1, M_FC_ICC_dbson_stn, padding);
M_FC_ICC_dbsoff_stn_padded = cat(1, M_FC_ICC_dbsoff_stn, padding);

M_BV_ICC_dbson_stn_padded = cat(1, M_BV_ICC_dbson_stn, padding);
M_BV_ICC_dbsoff_stn_padded = cat(1, M_BV_ICC_dbsoff_stn, padding);

% stack DBS ON/DBS OFF for grouping

FC_ICC_ON_stacked = [M_FC_ICC_dbson_stn_padded ; M_FC_ICC_dbson_gpi];
FC_ICC_OFF_stacked = [M_FC_ICC_dbsoff_stn_padded ; M_FC_ICC_dbsoff_gpi];
BV_ICC_ON_stacked = [M_BV_ICC_dbson_stn_padded ; M_BV_ICC_dbson_gpi];
BV_ICC_OFF_stacked = [M_BV_ICC_dbsoff_stn_padded ; M_BV_ICC_dbsoff_gpi];

group_inx_braintarget = [ones(1,length(gpi_subj)), 2.*ones(1,length(gpi_subj))]; 

group_names = {'STN', 'GPi'};

%FC ICC ON plot
figure('Color', 'w', 'Position',[0 0 1000 700]); 

subplot(2,2,1)
group_colors = [[.80, .345, .369]; [124/244, 177/255, 58/255]]; % Set colors for the DBS conditions
daboxplot(FC_ICC_ON_stacked,'groups',group_inx_braintarget,'color', group_colors,'whiskers' , 0, 'scatter', 1, 'scattersize', 80, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories,'fill',1,'legend',group_names,'boxspacing',0.6, 'boxwidth', 1.5);
legend('Location', 'southeast');
set(gca,'FontSize',15);
ylabel('FC ICC: ON');
ylim([0.3,0.9]);yticks(0.3:0.2:0.9);
xlim([0.5, 4.5]);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
set(gca, 'XTickLabel', categories, 'FontSize', 25, 'XTickLabelRotation', 30);
for j=1:4 %categories = {'whole', 'associative', 'limbic', 'motor'};
[p, h] = ranksum(squeeze(M_FC_ICC_dbson_stn(:,j)), squeeze(M_FC_ICC_dbson_gpi(:,j)));
if p<0.01
    str=sprintf('p= %1.2e',p);
else
    str=sprintf('p= %1.2f',p);
end
    effect = meanEffectSize(squeeze(M_FC_ICC_dbson_stn(:,j)), squeeze(M_FC_ICC_dbson_gpi(:,j)), Effect="cliff", paired=false, ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
    disp(effect)
    x_pos = j;
    T = text(x_pos, max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'center');
end

%FC ICC OFF plot
subplot(2,2,2)
daboxplot(FC_ICC_OFF_stacked,'groups',group_inx_braintarget,'color', group_colors,'whiskers' , 0, 'scatter', 1, 'scattersize', 80, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories,'fill',1,'legend',group_names,'boxspacing',0.6, 'boxwidth', 1.5);
legend('Location', 'southeast');
set(gca,'FontSize',15);
ylabel('FC ICC: OFF');
ylim([0.3,0.9]);yticks(0.3:0.2:0.9);
xlim([0.5, 4.5]);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
set(gca, 'XTickLabel', categories, 'FontSize', 25, 'XTickLabelRotation', 30);
for j=1:4 %categories = {'whole', 'associative', 'limbic', 'motor'};
    [p, h] = ranksum(squeeze(M_FC_ICC_dbsoff_stn(:,j)), squeeze(M_FC_ICC_dbsoff_gpi(:,j)));

    if p<0.01
        str=sprintf('p= %1.2e',p);
    else
        str=sprintf('p= %1.2f',p);
    end

    effect = meanEffectSize(squeeze(M_FC_ICC_dbsoff_stn(:,j)), squeeze(M_FC_ICC_dbsoff_gpi(:,j)), Effect="cliff", paired=false, ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
    disp(effect)

    x_pos = j;
    T = text(x_pos, max(get(gca, 'ylim')), str);
    set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'center');
end

%BV ICC ON plot
subplot(2,2,3)
group_colors = [[.80, .345, .369]; [124/244, 177/255, 58/255]]; % Set colors for the DBS conditions
daboxplot(BV_ICC_ON_stacked,'groups',group_inx_braintarget,'color', group_colors,'whiskers' , 0, 'scatter', 1, 'scattersize', 80, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories,'fill',1,'legend',group_names,'boxspacing',0.6, 'boxwidth', 1.5);
legend('Location', 'southeast');
set(gca,'FontSize',15);
ylabel('BV ICC: ON');
ylim([0.6,1.05]);yticks(0.6:0.2:1);
xlim([0.5, 4.5]);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
set(gca, 'XTickLabel', categories, 'FontSize', 25, 'XTickLabelRotation', 30);
for j=1:4 %categories = {'whole', 'associative', 'limbic', 'motor'};
[p, h] = ranksum(squeeze(M_BV_ICC_dbson_stn(:,j)), squeeze(M_BV_ICC_dbson_gpi(:,j)));

if p<0.01
    str=sprintf('p= %1.2e',p);
else
    str=sprintf('p= %1.2f',p);
end

    effect = meanEffectSize(squeeze(M_BV_ICC_dbson_stn(:,j)), squeeze(M_BV_ICC_dbson_gpi(:,j)), Effect="cliff", paired=false, ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
    disp(effect)
    x_pos = j;
    T = text(x_pos, max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'center');
end

%BV ICC OFF plot
subplot(2,2,4)
daboxplot(BV_ICC_OFF_stacked,'groups',group_inx_braintarget,'color', group_colors,'whiskers' , 0, 'scatter', 1, 'scattersize', 80, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories,'fill',1,'legend',group_names,'boxspacing',0.6, 'boxwidth', 1.5);
legend('Location', 'southeast');
set(gca,'FontSize',15);
ylabel('BV ICC: OFF');
ylim([0.6,1.05]);yticks(0.6:0.2:1);
xlim([0.5, 4.5]);
ax = gca; ax.Box = 'off';
ax.LineWidth = 3;
set(gca, 'XTickLabel', categories, 'FontSize', 25, 'XTickLabelRotation', 30);
for j=1:4 %categories = {'whole', 'associative', 'limbic', 'motor'};
[p, h] = ranksum(squeeze(M_BV_ICC_dbsoff_stn(:,j)), squeeze(M_BV_ICC_dbsoff_gpi(:,j)));
 
if p<0.01
    str=sprintf('p= %1.2e',p);
else
    str=sprintf('p= %1.2f',p);
end
    effect = meanEffectSize(squeeze(M_BV_ICC_dbsoff_stn(:,j)), squeeze(M_BV_ICC_dbsoff_gpi(:,j)), Effect="cliff", paired=false, ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
    disp(effect)
    x_pos = j;
    T = text(x_pos, max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 20, 'verticalalignment', 'top', 'horizontalalignment', 'center');
end