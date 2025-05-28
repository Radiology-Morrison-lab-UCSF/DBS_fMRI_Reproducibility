%% ICC Values for Functional Connectivity and Brain Variability Between Sessions
%% And Analysis of Longitudinal ICCs (Figure 4)

% This script calculates ICC values across longitudinal fMRI scans. % This script 
% is meant to be run section by section in order and used AFTER running the
% testretest_ICC.m script.


%% Add data paths

data_dir = '/YOUR/PATH/HERE'; % Replace with the path to your CONN output .mat files are stored

% Add paths to functions.
addpath('YOUR/PATH/HERE/ICC.m'); % Add a path to where the ICC function is stored. (Arash Salarian (2025). Intraclass Correlation Coefficient (ICC) (https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc), MATLAB Central File Exchange.)
addpath('YOUR/PATH/HERE/daboxplot/'); % Add a path to where the daboxplot function is stored. (Povilas Karvelis (2025). daboxplot (https://github.com/frank-pk/DataViz/releases/tag/v3.2.3), GitHub.)

% Add file names of the CONN output files for all networks for session 1 and session 2. 
results_whole_1 = 'results_whole_ses1.mat'; % whole brain network
results_motor_1 = 'results_motor_ses1.mat'; % motor network
results_limbic_1 = 'results_limbic_ses1.mat'; % limbic network
results_assoc_1 = 'results_assoc_ses1.mat'; % associative network
results_whole_2 = 'results_whole_ses1.mat'; % whole brain network
results_motor_2 = 'results_motor_ses1.mat'; % motor network
results_limbic_2 = 'results_limbic_ses1.mat'; % limbic network
results_assoc_2 = 'results_assoc_ses1.mat'; % associative network

%% Initialize code

dataset_1_name = 'V1'; % The name of the dataset for session 1 will be used to name the output Excel file.
dataset_2_name = 'V2'; % The name of dataset for session 2 which will be used to name the output Excel file.

brain_networks = {'whole','motor','limbic','associative'};

num_subj = 3; % Input the number of subjects.

ses1_order = [1,3,9]; % Specify the patient IDs in each session. 
ses2_order = [1,2,3]; 

ICC_type = 'A-1'; % Specify the ICC type, alpha value and r0 values.
alpha = 0.05; 
r0 = 0; 

%% -------------------------------------------------------------------------- %%

               %PART 1: BETWEEN SESSION ICC CALCULATION

%% -------------------------------------------------------------------------- %%

cd(data_dir)

%% Load connectomes

% whole brain v1
load(results_whole_1); % load functional connectivity data
brain_var_whole_1 = brain_var; % load brain variability data
fcs_whole_1 = fcs; 

for i = 1:size(fcs_whole_1,1) % loop over each subject
for l = 1:size(fcs_whole_1,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_whole_1(i,l,:,:) = atanh(fcs_whole_1(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

% whole brain v2
load(results_whole_2); % load functional connectivity data
brain_var_whole_2 = brain_var; % load brain variability data
fcs_whole_2 = fcs; 

for i = 1:size(fcs_whole_2,1) % loop over each subject
for l = 1:size(fcs_whole_2,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_whole_2(i,l,:,:) = atanh(fcs_whole_2(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

% motor v1
load(results_motor_1); % load functional connectivity data
brain_var_motor_1 = brain_var; % load brain variability data
fcs_motor_1 = fcs;

for i = 1:size(fcs_motor_1,1) % loop over each subject
for l = 1:size(fcs_motor_1,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_motor_1(i,l,:,:)=atanh(fcs_motor_1(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

% motor v2
load(results_motor_2); % load functional connectivity data
brain_var_motor_2 = brain_var; % load brain variability data
fcs_motor_2 = fcs;

for i = 1:size(fcs_motor_2,1) % loop over each subject
for l = 1:size(fcs_motor_2,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_motor_2(i,l,:,:)=atanh(fcs_motor_2(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

% limbic v1
load(results_limbic_1); % load functional connectivity data
brain_var_limbic_1 = brain_var; % load brain variability data
fcs_limbic_1 = fcs;

for i = 1:size(fcs_limbic_1,1) % loop over each subject
for l = 1:size(fcs_limbic_1,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_limbic_1(i,l,:,:) = atanh(fcs_limbic_1(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

% limbic v2
load(results_limbic_2); % load functional connectivity data
brain_var_limbic_2 = brain_var; % load brain variability data
fcs_limbic_2 = fcs;

for i = 1:size(fcs_limbic_2,1) % loop over each subject
for l = 1:size(fcs_limbic_2,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_limbic_2(i,l,:,:) = atanh(fcs_limbic_2(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

% associative v1
load(results_assoc_1); % load functional connectivity data
brain_var_associative_1 = brain_var; % load brain variability data
fcs_assoc_1=fcs;

for i = 1:size(fcs_assoc_1,1) % loop over each subject
for l = 1:size(fcs_assoc_1,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_assoc_1(i,l,:,:) = atanh(fcs_assoc_1(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

% associative v2
load(results_assoc_2); % load functional connectivity data
brain_var_associative_2 = brain_var; % load brain variability data
fcs_assoc_2=fcs;

for i = 1:size(fcs_assoc_2,1) % loop over each subject
for l = 1:size(fcs_assoc_2,2) % loop over each condition (OFFR1, OFFR2, ONR1, ONR2)
fcs_assoc_2(i,l,:,:) = atanh(fcs_assoc_2(i,l,:,:)); % fisher transform functional connectivity values for a normal distribution
end 
end
clear('fcs','brain_var');

%% Prepare functional connectivity matrices as vectors for correlation plots

% extract lower triangle from the matricies and remove center line of r=1

% v1
fcs_whole_tril_1 = zeros(size(fcs_whole_1));
fcs_motor_tril_1 = zeros(size(fcs_motor_1));
fcs_limbic_tril_1 = zeros(size(fcs_limbic_1));
fcs_associative_tril_1 = zeros(size(fcs_assoc_1));

% v2
fcs_whole_tril_2 = zeros(size(fcs_whole_2));
fcs_motor_tril_2 = zeros(size(fcs_motor_2));
fcs_limbic_tril_2 = zeros(size(fcs_limbic_2));
fcs_associative_tril_2 = zeros(size(fcs_assoc_2));

% v1
for i = 1:size(fcs_whole_1,1)
    for j = 1:size(fcs_whole_1,2)
fcs_whole_tril_1(i,j,:,:) = tril(squeeze(fcs_whole_1(i,j,:,:)),-1);
fcs_motor_tril_1(i,j,:,:) = tril(squeeze(fcs_motor_1(i,j,:,:)),-1);
fcs_limbic_tril_1(i,j,:,:) = tril(squeeze(fcs_limbic_1(i,j,:,:)),-1);
fcs_associative_tril_1(i,j,:,:) = tril(squeeze(fcs_assoc_1(i,j,:,:)),-1);
    end
end

% v2
for i = 1:size(fcs_whole_2,1)
    for j = 1:size(fcs_whole_2,2)
fcs_whole_tril_2(i,j,:,:) = tril(squeeze(fcs_whole_2(i,j,:,:)),-1);
fcs_motor_tril_2(i,j,:,:) = tril(squeeze(fcs_motor_2(i,j,:,:)),-1);
fcs_limbic_tril_2(i,j,:,:) = tril(squeeze(fcs_limbic_2(i,j,:,:)),-1);
fcs_associative_tril_2(i,j,:,:) = tril(squeeze(fcs_assoc_2(i,j,:,:)),-1);
    end
end

% plot 'motor' v1
figure;
imagesc(squeeze(fcs_motor_tril_1(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Motor Network V1'); 

% plot 'motor' v2
figure;
imagesc(squeeze(fcs_motor_tril_2(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Motor Network V2'); 

% plot 'whole brain' v1
figure;
imagesc(squeeze(fcs_whole_tril_1(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Wholebrain Network V1');

% plot 'whole brain' v2
figure;
imagesc(squeeze(fcs_whole_tril_2(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Wholebrain Network V2');

% plot 'associative' v1
figure; 
imagesc(squeeze(fcs_associative_tril_1(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Associative Network V1'); 

% plot 'associative' v2
figure; 
imagesc(squeeze(fcs_associative_tril_2(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Associative Network V2'); 

% plot 'limbic' v1
figure; 
imagesc(squeeze(fcs_limbic_tril_1(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Limbic Network V1');

% plot 'limbic' v2
figure; 
imagesc(squeeze(fcs_limbic_tril_2(i, j, :, :))); 
colormap default;  
colorbar; 
clim([-1, 1]);
title('Limbic Network V2');

% QA: Verify tril worked. If the plots show that only the lower triangular portion (excluding the diagonal) has non-zero colors and the upper triangular part is all zeros (a single uniform color representing zero), then the tril function has worked correctly.

% reshape as vectors instead of matrices excluding zeros v1
i=size(fcs_whole_1,1);
fcs_whole_tril_v1_reshape = zeros([size(fcs_whole_1,1) size(fcs_whole_1,2) nnz(fcs_whole_tril_1(i,j,:,:))]);
fcs_motor_tril_v1_reshape = zeros([size(fcs_motor_1,1) size(fcs_motor_1,2) nnz(fcs_motor_tril_1(i,j,:,:))]);
fcs_limbic_tril_v1_reshape = zeros([size(fcs_limbic_1,1) size(fcs_limbic_1,2) nnz(fcs_limbic_tril_1(i,j,:,:))]);
fcs_associative_tril_v1_reshape = zeros([size(fcs_assoc_1,1) size(fcs_assoc_1,2) nnz(fcs_associative_tril_1(i,j,:,:))]);

for i=1:size(fcs_whole_1,1)
    for j=1:size(fcs_whole_1,2)
temp = squeeze(fcs_whole_tril_1(i,j,:,:));
fcs_whole_tril_v1_reshape(i,j,:) = nonzeros(temp(:));
temp = squeeze(fcs_motor_tril_1(i,j,:,:));
fcs_motor_tril_v1_reshape(i,j,:) = nonzeros(temp(:));
temp = squeeze(fcs_limbic_tril_1(i,j,:,:));
fcs_limbic_tril_v1_reshape(i,j,:) = nonzeros(temp(:));
temp = squeeze(fcs_associative_tril_1(i,j,:,:));
fcs_associative_tril_v1_reshape(i,j,:) = nonzeros(temp(:));
    end
end

% reshape as vectors instead of matrices excluding zeros v2
i=size(fcs_whole_2,1);
fcs_whole_tril_v2_reshape = zeros([size(fcs_whole_2,1) size(fcs_whole_2,2) nnz(fcs_whole_tril_2(i,j,:,:))]);
fcs_motor_tril_v2_reshape = zeros([size(fcs_motor_2,1) size(fcs_motor_2,2) nnz(fcs_motor_tril_2(i,j,:,:))]);
fcs_limbic_tril_v2_reshape = zeros([size(fcs_limbic_2,1) size(fcs_limbic_2,2) nnz(fcs_limbic_tril_2(i,j,:,:))]);
fcs_associative_tril_v2_reshape = zeros([size(fcs_assoc_2,1) size(fcs_assoc_2,2) nnz(fcs_associative_tril_2(i,j,:,:))]);

for i=1:size(fcs_whole_2,1)
    for j=1:size(fcs_whole_2,2)
temp = squeeze(fcs_whole_tril_2(i,j,:,:));
fcs_whole_tril_v2_reshape(i,j,:) = nonzeros(temp(:));
temp = squeeze(fcs_motor_tril_2(i,j,:,:));
fcs_motor_tril_v2_reshape(i,j,:) = nonzeros(temp(:));
temp = squeeze(fcs_limbic_tril_2(i,j,:,:));
fcs_limbic_tril_v2_reshape(i,j,:) = nonzeros(temp(:));
temp = squeeze(fcs_associative_tril_2(i,j,:,:));
fcs_associative_tril_v2_reshape(i,j,:) = nonzeros(temp(:));
    end
end

%% Calculate and plot ICC values for functional connectivity and brain variability

% step 1: initialize matrices to store ICC results for each connectome
% FC = functional connectivity
% BV = brain variability
% 7 = 7 ICC output variables

% whole brain
FC_ICC_DBSOFF_whole_v1tv2rt = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
FC_ICC_DBSON_whole_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON
BV_ICC_DBSOFF_whole_v1tv2rt = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
BV_ICC_DBSON_whole_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON

FC_ICC_DBSOFF_whole_v1tv2t = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
FC_ICC_DBSON_whole_v1tv2t = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON
BV_ICC_DBSOFF_whole_v1tv2t = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
BV_ICC_DBSON_whole_v1tv2t = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON

FC_ICC_DBSOFF_whole_v1rtv2t = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
FC_ICC_DBSON_whole_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON
BV_ICC_DBSOFF_whole_v1rtv2t = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
BV_ICC_DBSON_whole_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON

FC_ICC_DBSOFF_whole_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
FC_ICC_DBSON_whole_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON
BV_ICC_DBSOFF_whole_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
BV_ICC_DBSON_whole_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON

% motor
FC_ICC_DBSOFF_motor_v1tv2rt = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
FC_ICC_DBSON_motor_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON
BV_ICC_DBSOFF_motor_v1tv2rt = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
BV_ICC_DBSON_motor_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON

FC_ICC_DBSOFF_motor_v1tv2t = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
FC_ICC_DBSON_motor_v1tv2t = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON
BV_ICC_DBSOFF_motor_v1tv2t = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
BV_ICC_DBSON_motor_v1tv2t = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON

FC_ICC_DBSOFF_motor_v1rtv2t = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
FC_ICC_DBSON_motor_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON
BV_ICC_DBSOFF_motor_v1rtv2t = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
BV_ICC_DBSON_motor_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON

FC_ICC_DBSOFF_motor_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
FC_ICC_DBSON_motor_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON
BV_ICC_DBSOFF_motor_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'motor' connectome, DBS OFF
BV_ICC_DBSON_motor_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'motor' connectome, DBS ON

% limbic
FC_ICC_DBSOFF_limbic_v1tv2rt = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
FC_ICC_DBSON_limbic_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON
BV_ICC_DBSOFF_limbic_v1tv2rt = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
BV_ICC_DBSON_limbic_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON

FC_ICC_DBSOFF_limbic_v1tv2t = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
FC_ICC_DBSON_limbic_v1tv2t = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON
BV_ICC_DBSOFF_limbic_v1tv2t = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
BV_ICC_DBSON_limbic_v1tv2t = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON

FC_ICC_DBSOFF_limbic_v1rtv2t = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
FC_ICC_DBSON_limbic_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON
BV_ICC_DBSOFF_limbic_v1rtv2t = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
BV_ICC_DBSON_limbic_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON

FC_ICC_DBSOFF_limbic_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
FC_ICC_DBSON_limbic_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON
BV_ICC_DBSOFF_limbic_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'limbic' connectome, DBS OFF
BV_ICC_DBSON_limbic_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'limbic' connectome, DBS ON

% associative
FC_ICC_DBSOFF_associative_v1tv2rt = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
FC_ICC_DBSON_associative_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON
BV_ICC_DBSOFF_associative_v1tv2rt = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
BV_ICC_DBSON_associative_v1tv2rt = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON

FC_ICC_DBSOFF_associative_v1tv2t = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
FC_ICC_DBSON_associative_v1tv2t = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON
BV_ICC_DBSOFF_associative_v1tv2t = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
BV_ICC_DBSON_associative_v1tv2t = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON

FC_ICC_DBSOFF_associative_v1rtv2t = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
FC_ICC_DBSON_associative_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON
BV_ICC_DBSOFF_associative_v1rtv2t = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
BV_ICC_DBSON_associative_v1rtv2t = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON

FC_ICC_DBSOFF_associative_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
FC_ICC_DBSON_associative_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON
BV_ICC_DBSOFF_associative_v1rtv2rt = zeros(num_subj, 7); % ICC results for 'associative' connectome, DBS OFF
BV_ICC_DBSON_associative_v1rtv2rt = zeros(num_subj, 7);  % ICC results for 'associative' connectome, DBS ON

% step 2: arrange rows (measures) and columns (rater or test-retest) & perform ICC

for k = 1:num_subj
    for net = 1:length(brain_networks)
        current_network = brain_networks{net}; % Get the current network name
        % functional Connectivity
        % construct variable names based on the current network
        M_OFF_v1tv2rt = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 1, :, :)'])) ...
                         nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 2, :, :)']))];
        disp(M_OFF_v1tv2rt)
        
        M_OFF_v1tv2t = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 1, :, :)'])) ...
                        nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 1, :, :)']))];
        disp(M_OFF_v1tv2t)
        
        M_OFF_v1rtv2t = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 2, :, :)'])) ...
                         nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 1, :, :)']))];
        disp(M_OFF_v1rtv2t)
        
        M_OFF_v1rtv2rt = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 2, :, :)'])) ...
                          nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 2, :, :)']))];
        disp(M_OFF_v1rtv2rt)

        % repeat the same for the ON condition
        M_ON_v1tv2rt = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 3, :, :)'])) ...
                        nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 4, :, :)']))];
        disp(M_ON_v1tv2rt)
        
        M_ON_v1tv2t = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 3, :, :)'])) ...
                       nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 3, :, :)']))];
        disp(M_ON_v1tv2t)
        
        M_ON_v1rtv2t = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 4, :, :)'])) ...
                        nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 3, :, :)']))];
        disp(M_ON_v1rtv2t)
        
        M_ON_v1rtv2rt = [nonzeros(eval(['fcs_' current_network '_tril_1(V1_order(k), 4, :, :)'])) ...
                         nonzeros(eval(['fcs_' current_network '_tril_2(V2_order(k), 4, :, :)']))];
        disp(M_ON_v1rtv2rt)

        % calculate ICC for OFF condition
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1tv2rt, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSOFF_' current_network '_v1tv2rt(k, :) = [r_off, LB, UB, F, df1, df2, p];']);
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1tv2t, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSOFF_' current_network '_v1tv2t(k, :) = [r_off, LB, UB, F, df1, df2, p];']);
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1rtv2t, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSOFF_' current_network '_v1rtv2t(k, :) = [r_off, LB, UB, F, df1, df2, p];']);
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1rtv2rt, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSOFF_' current_network '_v1rtv2rt(k, :) = [r_off, LB, UB, F, df1, df2, p];']);

        % calculate ICC for ON condition
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1tv2rt, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSON_' current_network '_v1tv2rt(k, :) = [r_on, LB, UB, F, df1, df2, p];']);
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1tv2t, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSON_' current_network '_v1tv2t(k, :) = [r_on, LB, UB, F, df1, df2, p];']);
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1rtv2t, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSON_' current_network '_v1rtv2t(k, :) = [r_on, LB, UB, F, df1, df2, p];']);
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1rtv2rt, ICC_type, alpha, r0);
        eval(['FC_ICC_DBSON_' current_network '_v1rtv2rt(k, :) = [r_on, LB, UB, F, df1, df2, p];']);
    
        % brain variability
        % construct variable names based on the current network
        M_OFF_v1tv2rt = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 1, :, :)'])) ...
                         nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 2, :, :)']))];
        disp(M_OFF_v1tv2rt)
        
        M_OFF_v1tv2t = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 1, :, :)'])) ...
                        nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 1, :, :)']))];
        disp(M_OFF_v1tv2t)
        
        M_OFF_v1rtv2t = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 2, :, :)'])) ...
                         nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 1, :, :)']))];
        disp(M_OFF_v1rtv2t)
        
        M_OFF_v1rtv2rt = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 2, :, :)'])) ...
                          nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 2, :, :)']))];
        disp(M_OFF_v1rtv2rt)

        % repeat the same for the ON condition
        M_ON_v1tv2rt = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 3, :, :)'])) ...
                        nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 4, :, :)']))];
        disp(M_ON_v1tv2rt)
        
        M_ON_v1tv2t = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 3, :, :)'])) ...
                       nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 3, :, :)']))];
        disp(M_ON_v1tv2t)
        
        M_ON_v1rtv2t = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 4, :, :)'])) ...
                        nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 3, :, :)']))];
        disp(M_ON_v1rtv2t)
        
        M_ON_v1rtv2rt = [nonzeros(eval(['brain_var_' current_network '_1(V1_order(k), 4, :, :)'])) ...
                         nonzeros(eval(['brain_var_' current_network '_2(V2_order(k), 4, :, :)']))];
        disp(M_ON_v1rtv2rt)

        % calculate ICC for OFF condition
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1tv2rt, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSOFF_' current_network '_v1tv2rt(k, :) = [r_off, LB, UB, F, df1, df2, p];']);
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1tv2t, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSOFF_' current_network '_v1tv2t(k, :) = [r_off, LB, UB, F, df1, df2, p];']);
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1rtv2t, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSOFF_' current_network '_v1rtv2t(k, :) = [r_off, LB, UB, F, df1, df2, p];']);
        [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF_v1rtv2rt, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSOFF_' current_network '_v1rtv2rt(k, :) = [r_off, LB, UB, F, df1, df2, p];']);

        % calculate ICC for ON condition
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1tv2rt, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSON_' current_network '_v1tv2rt(k, :) = [r_on, LB, UB, F, df1, df2, p];']);
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1tv2t, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSON_' current_network '_v1tv2t(k, :) = [r_on, LB, UB, F, df1, df2, p];']);
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1rtv2t, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSON_' current_network '_v1rtv2t(k, :) = [r_on, LB, UB, F, df1, df2, p];']);
        [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON_v1rtv2rt, ICC_type, alpha, r0);
        eval(['BV_ICC_DBSON_' current_network '_v1rtv2rt(k, :) = [r_on, LB, UB, F, df1, df2, p];']);

    end
end

%% Save functional connectivity ICC values to excel file

% functional connectivity ICC values
% step 1: convert ICC matrices to tables

% define header for ICC table
fc_icc_header = {'FC ICC', 'Lower Bound', 'Upper Bound', 'F', 'df1', 'df2', 'p-value'};

% create tables for each connectome type and DBS condition

% whole
table_FC_ICC_DBSOFF_whole_v1rtv2rt = array2table(FC_ICC_DBSOFF_whole_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_whole_v1tv2t = array2table(FC_ICC_DBSOFF_whole_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_whole_v1rtv2t = array2table(FC_ICC_DBSOFF_whole_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_whole_v1tv2rt = array2table(FC_ICC_DBSOFF_whole_v1rtv2rt, 'VariableNames', fc_icc_header);

table_FC_ICC_DBSON_whole_v1rtv2rt = array2table(FC_ICC_DBSON_whole_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_whole_v1tv2t = array2table(FC_ICC_DBSON_whole_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_whole_v1rtv2t = array2table(FC_ICC_DBSON_whole_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_whole_v1tv2rt = array2table(FC_ICC_DBSON_whole_v1rtv2rt, 'VariableNames', fc_icc_header);

% motor
table_FC_ICC_DBSOFF_motor_v1rtv2rt = array2table(FC_ICC_DBSOFF_motor_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_motor_v1tv2t = array2table(FC_ICC_DBSOFF_motor_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_motor_v1rtv2t = array2table(FC_ICC_DBSOFF_motor_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_motor_v1tv2rt = array2table(FC_ICC_DBSOFF_motor_v1rtv2rt, 'VariableNames', fc_icc_header);

table_FC_ICC_DBSON_motor_v1rtv2rt = array2table(FC_ICC_DBSON_motor_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_motor_v1tv2t = array2table(FC_ICC_DBSON_motor_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_motor_v1rtv2t = array2table(FC_ICC_DBSON_motor_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_motor_v1tv2rt = array2table(FC_ICC_DBSON_motor_v1rtv2rt, 'VariableNames', fc_icc_header);

% limbic
table_FC_ICC_DBSOFF_limbic_v1rtv2rt = array2table(FC_ICC_DBSOFF_limbic_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_limbic_v1tv2t = array2table(FC_ICC_DBSOFF_limbic_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_limbic_v1rtv2t = array2table(FC_ICC_DBSOFF_limbic_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_limbic_v1tv2rt = array2table(FC_ICC_DBSOFF_limbic_v1rtv2rt, 'VariableNames', fc_icc_header);

table_FC_ICC_DBSON_limbic_v1rtv2rt = array2table(FC_ICC_DBSON_limbic_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_limbic_v1tv2t = array2table(FC_ICC_DBSON_limbic_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_limbic_v1rtv2t = array2table(FC_ICC_DBSON_limbic_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_limbic_v1tv2rt = array2table(FC_ICC_DBSON_limbic_v1rtv2rt, 'VariableNames', fc_icc_header);

% associative
table_FC_ICC_DBSOFF_associative_v1rtv2rt = array2table(FC_ICC_DBSOFF_associative_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_associative_v1tv2t = array2table(FC_ICC_DBSOFF_associative_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_associative_v1rtv2t = array2table(FC_ICC_DBSOFF_associative_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSOFF_associative_v1tv2rt = array2table(FC_ICC_DBSOFF_associative_v1rtv2rt, 'VariableNames', fc_icc_header);

table_FC_ICC_DBSON_associative_v1rtv2rt = array2table(FC_ICC_DBSON_associative_v1tv2rt, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_associative_v1tv2t = array2table(FC_ICC_DBSON_associative_v1tv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_associative_v1rtv2t = array2table(FC_ICC_DBSON_associative_v1rtv2t, 'VariableNames', fc_icc_header);
table_FC_ICC_DBSON_associative_v1tv2rt = array2table(FC_ICC_DBSON_associative_v1rtv2rt, 'VariableNames', fc_icc_header);

% step 2: write tables to Excel file

% define the filename for the Excel file
filename = [dataset_1_name '_' dataset_2_name '_FC_ICC_Between_Results.xlsx'];

% write each table to a separate sheet in the Excel file

% whole
writetable(table_FC_ICC_DBSOFF_whole_v1rtv2rt, filename, 'Sheet', 'whole_DBSOFF_v1rtv2rt');
writetable(table_FC_ICC_DBSOFF_whole_v1tv2t, filename, 'Sheet', 'whole_DBSOFF_v1tv2t');
writetable(table_FC_ICC_DBSOFF_whole_v1rtv2t, filename, 'Sheet', 'whole_DBSOFF_v1rtv2t');
writetable(table_FC_ICC_DBSOFF_whole_v1tv2rt, filename, 'Sheet', 'whole_DBSOFF_v1tv2rt');

writetable(table_FC_ICC_DBSON_whole_v1rtv2rt, filename, 'Sheet', 'whole_DBSON_v1rtv2rt');
writetable(table_FC_ICC_DBSON_whole_v1tv2t, filename, 'Sheet', 'whole_DBSON_v1tv2t');
writetable(table_FC_ICC_DBSON_whole_v1rtv2t, filename, 'Sheet', 'whole_DBSON_v1rtv2t');
writetable(table_FC_ICC_DBSON_whole_v1tv2rt, filename, 'Sheet', 'whole_DBSON_v1tv2rt');

% motor
writetable(table_FC_ICC_DBSOFF_motor_v1rtv2rt, filename, 'Sheet', 'motor_DBSOFF_v1rtv2rt');
writetable(table_FC_ICC_DBSOFF_motor_v1tv2t, filename, 'Sheet', 'motor_DBSOFF_v1tv2t');
writetable(table_FC_ICC_DBSOFF_motor_v1rtv2t, filename, 'Sheet', 'motor_DBSOFF_v1rtv2t');
writetable(table_FC_ICC_DBSOFF_motor_v1tv2rt, filename, 'Sheet', 'motor_DBSOFF_v1tv2rt');

writetable(table_FC_ICC_DBSON_motor_v1rtv2rt, filename, 'Sheet', 'motor_DBSON_v1rtv2rt');
writetable(table_FC_ICC_DBSON_motor_v1tv2t, filename, 'Sheet', 'motor_DBSON_v1tv2t');
writetable(table_FC_ICC_DBSON_motor_v1rtv2t, filename, 'Sheet', 'motor_DBSON_v1rtv2t');
writetable(table_FC_ICC_DBSON_motor_v1tv2rt, filename, 'Sheet', 'motor_DBSON_v1tv2rt');

% limbic
writetable(table_FC_ICC_DBSOFF_limbic_v1rtv2rt, filename, 'Sheet', 'limbic_DBSOFF_v1rtv2rt');
writetable(table_FC_ICC_DBSOFF_limbic_v1tv2t, filename, 'Sheet', 'limbic_DBSOFF_v1tv2t');
writetable(table_FC_ICC_DBSOFF_limbic_v1rtv2t, filename, 'Sheet', 'limbic_DBSOFF_v1rtv2t');
writetable(table_FC_ICC_DBSOFF_limbic_v1tv2rt, filename, 'Sheet', 'limbic_DBSOFF_v1tv2rt');

writetable(table_FC_ICC_DBSON_limbic_v1rtv2rt, filename, 'Sheet', 'limbic_DBSON_v1rtv2rt');
writetable(table_FC_ICC_DBSON_limbic_v1tv2t, filename, 'Sheet', 'limbic_DBSON_v1tv2t');
writetable(table_FC_ICC_DBSON_limbic_v1rtv2t, filename, 'Sheet', 'limbic_DBSON_v1rtv2t');
writetable(table_FC_ICC_DBSON_limbic_v1tv2rt, filename, 'Sheet', 'limbic_DBSON_v1tv2rt');

% associative
writetable(table_FC_ICC_DBSOFF_associative_v1rtv2rt, filename, 'Sheet', 'associative_DBSOFF_v1rtv2rt');
writetable(table_FC_ICC_DBSOFF_associative_v1tv2t, filename, 'Sheet', 'associative_DBSOFF_v1tv2t');
writetable(table_FC_ICC_DBSOFF_associative_v1rtv2t, filename, 'Sheet', 'associative_DBSOFF_v1rtv2t');
writetable(table_FC_ICC_DBSOFF_associative_v1tv2rt, filename, 'Sheet', 'associative_DBSOFF_v1tv2rt');

writetable(table_FC_ICC_DBSON_associative_v1rtv2rt, filename, 'Sheet', 'associative_DBSON_v1rtv2rt');
writetable(table_FC_ICC_DBSON_associative_v1tv2t, filename, 'Sheet', 'associative_DBSON_v1tv2t');
writetable(table_FC_ICC_DBSON_associative_v1rtv2t, filename, 'Sheet', 'associative_DBSON_v1rtv2t');
writetable(table_FC_ICC_DBSON_associative_v1tv2rt, filename, 'Sheet', 'associative_DBSON_v1tv2rt');

%% Save brain variability ICC values to excel file

% brain variability ICC values
% step 1: convert ICC matrices to tables

% define header for ICC table
bv_icc_header = {'BV ICC', 'Lower Bound', 'Upper Bound', 'F', 'df1', 'df2', 'p-value'};

% create tables for each connectome type and DBS condition

% whole
table_BV_ICC_DBSOFF_whole_v1rtv2rt = array2table(BV_ICC_DBSOFF_whole_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_whole_v1tv2t = array2table(BV_ICC_DBSOFF_whole_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_whole_v1rtv2t = array2table(BV_ICC_DBSOFF_whole_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_whole_v1tv2rt = array2table(BV_ICC_DBSOFF_whole_v1rtv2rt, 'VariableNames', bv_icc_header);

table_BV_ICC_DBSON_whole_v1rtv2rt = array2table(BV_ICC_DBSON_whole_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_whole_v1tv2t = array2table(BV_ICC_DBSON_whole_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_whole_v1rtv2t = array2table(BV_ICC_DBSON_whole_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_whole_v1tv2rt = array2table(BV_ICC_DBSON_whole_v1rtv2rt, 'VariableNames', bv_icc_header);

% motor
table_BV_ICC_DBSOFF_motor_v1rtv2rt = array2table(BV_ICC_DBSOFF_motor_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_motor_v1tv2t = array2table(BV_ICC_DBSOFF_motor_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_motor_v1rtv2t = array2table(BV_ICC_DBSOFF_motor_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_motor_v1tv2rt = array2table(BV_ICC_DBSOFF_motor_v1rtv2rt, 'VariableNames', bv_icc_header);

table_BV_ICC_DBSON_motor_v1rtv2rt = array2table(BV_ICC_DBSON_motor_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_motor_v1tv2t = array2table(BV_ICC_DBSON_motor_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_motor_v1rtv2t = array2table(BV_ICC_DBSON_motor_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_motor_v1tv2rt = array2table(BV_ICC_DBSON_motor_v1rtv2rt, 'VariableNames', bv_icc_header);

% limbic
table_BV_ICC_DBSOFF_limbic_v1rtv2rt = array2table(BV_ICC_DBSOFF_limbic_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_limbic_v1tv2t = array2table(BV_ICC_DBSOFF_limbic_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_limbic_v1rtv2t = array2table(BV_ICC_DBSOFF_limbic_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_limbic_v1tv2rt = array2table(BV_ICC_DBSOFF_limbic_v1rtv2rt, 'VariableNames', bv_icc_header);

table_BV_ICC_DBSON_limbic_v1rtv2rt = array2table(BV_ICC_DBSON_limbic_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_limbic_v1tv2t = array2table(BV_ICC_DBSON_limbic_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_limbic_v1rtv2t = array2table(BV_ICC_DBSON_limbic_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_limbic_v1tv2rt = array2table(BV_ICC_DBSON_limbic_v1rtv2rt, 'VariableNames', bv_icc_header);

% associative
table_BV_ICC_DBSOFF_associative_v1rtv2rt = array2table(BV_ICC_DBSOFF_associative_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_associative_v1tv2t = array2table(BV_ICC_DBSOFF_associative_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_associative_v1rtv2t = array2table(BV_ICC_DBSOFF_associative_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSOFF_associative_v1tv2rt = array2table(BV_ICC_DBSOFF_associative_v1rtv2rt, 'VariableNames', bv_icc_header);

table_BV_ICC_DBSON_associative_v1rtv2rt = array2table(BV_ICC_DBSON_associative_v1tv2rt, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_associative_v1tv2t = array2table(BV_ICC_DBSON_associative_v1tv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_associative_v1rtv2t = array2table(BV_ICC_DBSON_associative_v1rtv2t, 'VariableNames', bv_icc_header);
table_BV_ICC_DBSON_associative_v1tv2rt = array2table(BV_ICC_DBSON_associative_v1rtv2rt, 'VariableNames', bv_icc_header);

% step 2: write tables to Excel file

% define the filename for the Excel file
filename = [dataset_1_name '_' dataset_2_name '_BV_ICC_Between_Results.xlsx'];

% write each table to a separate sheet in the Excel file

% whole
writetable(table_BV_ICC_DBSOFF_whole_v1rtv2rt, filename, 'Sheet', 'whole_DBSOFF_v1rtv2rt');
writetable(table_BV_ICC_DBSOFF_whole_v1tv2t, filename, 'Sheet', 'whole_DBSOFF_v1tv2t');
writetable(table_BV_ICC_DBSOFF_whole_v1rtv2t, filename, 'Sheet', 'whole_DBSOFF_v1rtv2t');
writetable(table_BV_ICC_DBSOFF_whole_v1tv2rt, filename, 'Sheet', 'whole_DBSOFF_v1tv2rt');

writetable(table_BV_ICC_DBSON_whole_v1rtv2rt, filename, 'Sheet', 'whole_DBSON_v1rtv2rt');
writetable(table_BV_ICC_DBSON_whole_v1tv2t, filename, 'Sheet', 'whole_DBSON_v1tv2t');
writetable(table_BV_ICC_DBSON_whole_v1rtv2t, filename, 'Sheet', 'whole_DBSON_v1rtv2t');
writetable(table_BV_ICC_DBSON_whole_v1tv2rt, filename, 'Sheet', 'whole_DBSON_v1tv2rt');

% motor
writetable(table_BV_ICC_DBSOFF_motor_v1rtv2rt, filename, 'Sheet', 'motor_DBSOFF_v1rtv2rt');
writetable(table_BV_ICC_DBSOFF_motor_v1tv2t, filename, 'Sheet', 'motor_DBSOFF_v1tv2t');
writetable(table_BV_ICC_DBSOFF_motor_v1rtv2t, filename, 'Sheet', 'motor_DBSOFF_v1rtv2t');
writetable(table_BV_ICC_DBSOFF_motor_v1tv2rt, filename, 'Sheet', 'motor_DBSOFF_v1tv2rt');

writetable(table_BV_ICC_DBSON_motor_v1rtv2rt, filename, 'Sheet', 'motor_DBSON_v1rtv2rt');
writetable(table_BV_ICC_DBSON_motor_v1tv2t, filename, 'Sheet', 'motor_DBSON_v1tv2t');
writetable(table_BV_ICC_DBSON_motor_v1rtv2t, filename, 'Sheet', 'motor_DBSON_v1rtv2t');
writetable(table_BV_ICC_DBSON_motor_v1tv2rt, filename, 'Sheet', 'motor_DBSON_v1tv2rt');

% limbic
writetable(table_BV_ICC_DBSOFF_limbic_v1rtv2rt, filename, 'Sheet', 'limbic_DBSOFF_v1rtv2rt');
writetable(table_BV_ICC_DBSOFF_limbic_v1tv2t, filename, 'Sheet', 'limbic_DBSOFF_v1tv2t');
writetable(table_BV_ICC_DBSOFF_limbic_v1rtv2t, filename, 'Sheet', 'limbic_DBSOFF_v1rtv2t');
writetable(table_BV_ICC_DBSOFF_limbic_v1tv2rt, filename, 'Sheet', 'limbic_DBSOFF_v1tv2rt');

writetable(table_BV_ICC_DBSON_limbic_v1rtv2rt, filename, 'Sheet', 'limbic_DBSON_v1rtv2rt');
writetable(table_BV_ICC_DBSON_limbic_v1tv2t, filename, 'Sheet', 'limbic_DBSON_v1tv2t');
writetable(table_BV_ICC_DBSON_limbic_v1rtv2t, filename, 'Sheet', 'limbic_DBSON_v1rtv2t');
writetable(table_BV_ICC_DBSON_limbic_v1tv2rt, filename, 'Sheet', 'limbic_DBSON_v1tv2rt');

% associative
writetable(table_BV_ICC_DBSOFF_associative_v1rtv2rt, filename, 'Sheet', 'associative_DBSOFF_v1rtv2rt');
writetable(table_BV_ICC_DBSOFF_associative_v1tv2t, filename, 'Sheet', 'associative_DBSOFF_v1tv2t');
writetable(table_BV_ICC_DBSOFF_associative_v1rtv2t, filename, 'Sheet', 'associative_DBSOFF_v1rtv2t');
writetable(table_BV_ICC_DBSOFF_associative_v1tv2rt, filename, 'Sheet', 'associative_DBSOFF_v1tv2rt');

writetable(table_BV_ICC_DBSON_associative_v1rtv2rt, filename, 'Sheet', 'associative_DBSON_v1rtv2rt');
writetable(table_BV_ICC_DBSON_associative_v1tv2t, filename, 'Sheet', 'associative_DBSON_v1tv2t');
writetable(table_BV_ICC_DBSON_associative_v1rtv2t, filename, 'Sheet', 'associative_DBSON_v1rtv2t');
writetable(table_BV_ICC_DBSON_associative_v1tv2rt, filename, 'Sheet', 'associative_DBSON_v1tv2rt');


%% -------------------------------------------------------------------------- %%

               %PART 2: Analysis - Figure 4

%% -------------------------------------------------------------------------- %%

%% Figure 4B - Plot a scatter plot of ICC ses01 against ICC ses02 ON OFF functional connectivity/brain variability within values

% define the Excel file name
fileName_1 = 'V1_FC_ICC_Results.xlsx'; % change to BV/FC to plot other variable
fileName_2 = 'V2_FC_ICC_Results.xlsx'; % change to BV/FC to plot other variable
% list of sheet names
sheetNames = {'Whole_DBS_OFF', 'Whole_DBS_ON', 'Motor_DBS_OFF', 'Motor_DBS_ON',...
    'Limbic_DBS_OFF', 'Limbic_DBS_ON', 'Associative_DBS_OFF', 'Associative_DBS_ON'};
labels = {'Whole', 'Whole', 'Motor', 'Motor',...
    'Limbic', 'Limbic', 'Associative', 'Associative'};

% preallocate cell arrays to store data
dataV1 = cell(1, length(sheetNames));
dataV2 = cell(1, length(sheetNames));

% define patient indices to match between V1 and V2
patientsV1_order = [1, 3, 9]; 
patientsV2_order = [1, 2, 3]; 

% loop through each sheet and extract relevant patients' data
for i = 1:length(sheetNames)
    % read the first column from the current sheet of both files
    v1_data = readmatrix(fileName_1, 'Sheet', sheetNames{i}, 'Range', 'A2:A17');  % start from A2 due to header
    v2_data = readmatrix(fileName_2, 'Sheet', sheetNames{i}, 'Range', 'A2:A17');  % start from A2 due to header

    % extract ICC values based on patient order
    dataV1{i} = v1_data(patientsV1_order); 
    dataV2{i} = v2_data(patientsV2_order);
end

% define the color for each network
colors_OFF = [0 0.45 0.7]; % darker blue for OFF state
colors_ON = [0.5 0.8 1];  % light blue for ON state

% define marker styles for different networks
markers = {'o', 'o', 'square', 'square', '^', '^', 'diamond', 'diamond'};  % markers for whole, motor, limbic, and associative networks

% create scatter plot
figure;
hold on;

% initialize arrays for regression analysis
allV1_ON = [];
allV2_ON = [];
allV1_OFF = [];
allV2_OFF = [];

% loop through each sheet to plot data
for i = 1:length(sheetNames)
    % check if the sheet corresponds to 'ON' or 'OFF'
    if contains(sheetNames{i}, 'ON')
        % plot using black corresponding markers for 'ON'
        markerType = markers{mod(i-1, length(markers)) + 1};  % select marker based on network type
        scatter(dataV2{i}, dataV1{i}, 400, 'filled', markerType, 'DisplayName', labels{i}, 'MarkerFaceColor', colors_ON);
        % collect data for regression
        allV1_ON = [allV1_ON; dataV1{i}];
        allV2_ON = [allV2_ON; dataV2{i}];

    elseif contains(sheetNames{i}, 'OFF')
        % plot using red corresponding markers for 'OFF'
        markerType = markers{mod(i-1, length(markers)) + 1};  % select marker based on network type
        scatter(dataV2{i}, dataV1{i}, 400, 'filled', markerType, 'DisplayName', labels{i}, 'MarkerFaceColor', colors_OFF);
        % collect data for regression
        allV1_OFF = [allV1_OFF; dataV1{i}];
        allV2_OFF = [allV2_OFF; dataV2{i}];
    end
end

% perform linear fit for ON and OFF combined
P_ON = polyfit(allV2_ON, allV1_ON, 1); % linear fit for ON
yfit_ON = polyval(P_ON, allV2_ON); % evaluate the fit

P_OFF = polyfit(allV2_OFF, allV1_OFF, 1); % linear fit for OFF
yfit_OFF = polyval(P_OFF, allV2_OFF); % evaluate the fit

% plot regression lines
plot(allV2_OFF, yfit_OFF, 'Color', [0 0.45 0.7], 'LineWidth', 1.5, 'DisplayName', 'Fit OFF'); % fit line for OFF
plot(allV2_ON, yfit_ON, 'Color', [0.5 0.8 1], 'LineWidth', 1.5, 'DisplayName', 'Fit ON'); % fit line for ON

% calculate and display correlation coefficient (r) and p-value
[r_OFF, p_OFF] = corr(allV2_OFF, allV1_OFF);
str_off = sprintf('r_{OFF} = %.2f, p_{OFF} = %.1e', r_OFF, p_OFF);
T_off = text(0.7, 1, str_off, 'FontSize', 19);

[r_ON, p_ON] = corr(allV2_ON ,allV1_ON);
str_on = sprintf('r_{ON} = %.2f, p_{ON} = %.1e', r_ON, p_ON);
T_on = text(0.7, 0.95, str_on, 'FontSize', 19);

% add labels and legend
xlabel('ICC: V2', 'FontSize', 25);
xticks([0.7 0.75 0.8 0.85 0.9 0.95 1]); % edit depending on axis size
xlim([0.7 1]); % edit depending on axis size
ylabel('ICC: V1', 'FontSize', 25);
yticks([0.7 0.75 0.8 0.85 0.9 0.95 1]); % edit depending on axis size
ylim([0.7 1]);

% set other properties
set(gca, 'FontSize', 25); % increase axis font size
set(gcf, 'Color', 'w'); % optional: set figure background color


[r_ON, p_ON, rl_ON, ru_ON] = corrcoef(allV2_ON, allV1_ON);
    disp([r_ON(1,2) p_ON(1,2) rl_ON(1,2) ru_ON(1,2)]);

[r_OFF,p_OFF,rl_OFF, ru_OFF] = corrcoef(allV2_OFF, allV1_OFF);
    disp([r_OFF(1,2) p_OFF(1,2) rl_OFF(1,2) ru_OFF(1,2)]);


%% Legend for scatter plot above

% create a figure
figure;

% hold the plot to add legend entries
hold on;

% create dummy plot objects to populate the legend
h1 = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 15, 'LineWidth', 2); % for whole brain
h2 = plot(NaN, NaN, 'square', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 15, 'LineWidth', 2);  % for motor
h3 = plot(NaN, NaN, '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none','MarkerSize', 15, 'LineWidth', 2); % for limbic
h4 = plot(NaN, NaN, 'diamond', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none','MarkerSize', 15, 'LineWidth', 2); % for associative

% create the legend
legend([h1, h2, h3, h4], {'Whole Brain', 'Motor', 'Limbic', 'Associative'}, ...
    'Location', 'northeastoutside', 'FontSize', 20);

% add labels and axes properties
xlabel('X-axis Label', 'FontSize', 20);
ylabel('Y-axis Label', 'FontSize', 20);
set(gca, 'FontSize', 12); % set axis font size

%% Plot box plot for functional connectivity and/or brain variability of within and between ICC values of ON OFF DBS 
% x = [w, b (on)] [w, b (off)] 
% y = ICC

% define the Excel file name
filename_1 = 'V1_BV_ICC_Results.xlsx'; % within FC/BV ICC values visit 1
filename_2 = 'V2_BV_ICC_Results.xlsx'; % within FC/BV ICC values visit 2
filename_3 = 'V1_V2_BV_ICC_Between_Results.xlsx'; % between FC/BV ICC values visit 1 & 2

% display sheet names for each file:
sheets_1 = sheetnames('V1_BV_ICC_Results.xlsx'); % for file 1
disp('Sheet names:');
disp(sheets_1);
sheets_2 = sheetnames('V2_BV_ICC_Results.xlsx'); % for file 2
disp('Sheet names:');
disp(sheets_2);
sheets_3 = sheetnames('V1_V2_BV_ICC_Between_Results.xlsx'); % for file 3
disp('Sheet names:');
disp(sheets_3);

% list of sheet names
sheetnames_1_2 = {'Whole_DBS_OFF', 'Whole_DBS_ON', 'Motor_DBS_OFF', 'Motor_DBS_ON',...
    'Limbic_DBS_OFF', 'Limbic_DBS_ON', 'Associative_DBS_OFF', 'Associative_DBS_ON'}; % filename_1 & filemame_2
sheetnames_3 = {'whole_DBSOFF_v1rtv2rt', 'whole_DBSOFF_v1tv2t',	'whole_DBSOFF_v1rtv2t',...
    'whole_DBSOFF_v1tv2rt',	'whole_DBSON_v1rtv2rt',	'whole_DBSON_v1tv2t','whole_DBSON_v1rtv2t',...
    'whole_DBSON_v1tv2rt',	'motor_DBSOFF_v1rtv2rt', 'motor_DBSOFF_v1tv2t', 'motor_DBSOFF_v1rtv2t',...
    'motor_DBSOFF_v1tv2rt',	'motor_DBSON_v1rtv2rt',	'motor_DBSON_v1tv2t', 'motor_DBSON_v1rtv2t',...
    'motor_DBSON_v1tv2rt',	'limbic_DBSOFF_v1rtv2rt', 'limbic_DBSOFF_v1tv2t', 'limbic_DBSOFF_v1rtv2t',...
    'limbic_DBSOFF_v1tv2rt', 'limbic_DBSON_v1rtv2rt', 'limbic_DBSON_v1tv2t', 'limbic_DBSON_v1rtv2t',...
    'limbic_DBSON_v1tv2rt',	'associative_DBSOFF_v1rtv2rt', 'associative_DBSOFF_v1tv2t', 'associative_DBSOFF_v1rtv2t',...
    'associative_DBSOFF_v1tv2rt', 'associative_DBSON_v1rtv2rt',	'associative_DBSON_v1tv2t', 'associative_DBSON_v1rtv2t',...
    'associative_DBSON_v1tv2rt'}; % filename_3

% labels for sheet names
labels_1_2 = {'Whole DBS OFF', 'Whole DBS ON', 'Motor DBS OFF', 'Motor DBS ON',...
    'Limbic DBS OFF', 'Limbic DBS ON', 'Associative DBS OFF', 'Associative DBS ON'}; % filename_1 & filemame_2
labels_3 = {'Whole DBSOFF V1RTV2RT', 'Whole DBSOFF V1TV2T', 'Whole DBSOFF V1RTV2T', 'Whole DBSOFF V1TV2RT',... 
'Whole DBSON V1RTV2RT', 'Whole DBSON V1TV2T', 'Whole DBSON V1RTV2T', 'Whole DBSON V1TV2RT',... 
'Motor DBSOFF V1RTV2RT', 'Motor DBSOFF V1TV2T', 'Motor DBSOFF V1RTV2T', 'Motor DBSOFF V1TV2RT', ...
'Motor DBSON V1RTV2RT', 'Motor DBSON V1TV2T', 'Motor DBSON V1RTV2T', 'Motor DBSON V1TV2RT',... 
'Limbic DBSOFF V1RTV2RT', 'Limbic DBSOFF V1TV2T', 'Limbic DBSOFF V1RTV2T', 'Limbic DBSOFF V1TV2RT', ...
'Limbic DBSON V1RTV2RT', 'Limbic DBSON V1TV2T', 'Limbic DBSON V1RTV2T', 'Limbic DBSON V1TV2RT',... 
'Associative DBSOFF V1RTV2RT', 'Associative DBSOFF V1TV2T', 'Associative DBSOFF V1RTV2T',... 
'Associative DBSOFF V1TV2RT', 'Associative DBSON V1RTV2RT', 'Associative DBSON V1TV2T', ...
'Associative DBSON V1RTV2T', 'Associative DBSON V1TV2RT'}; % filename_3

% initialize arrays to hold the data
within_ON = [];
between_ON = [];
within_OFF = [];
between_OFF = [];

% read within ON values from both file filename_1 and filename_2 and combine
for i = 1:length(sheetnames_1_2)
    % check if the sheet corresponds to 'ON'
    if contains(sheetnames_1_2{i}, 'ON')
        % read data from filename_1
        data_ON_filename_1 = readmatrix(filename_1, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % read data from filename_2
        data_ON_filename_2 = readmatrix(filename_2, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % combine values into within_ON array
        within_ON = [within_ON; data_ON_filename_1; data_ON_filename_2]; % makes within_ON 24x1 double
    end
end

% read within OFF values from both file filename_1 and filename_2 and combine
for i = 1:length(sheetnames_1_2)
    % check if the sheet corresponds to 'OFF'
    if contains(sheetnames_1_2{i}, 'OFF')
        % read data from filename_1
        data_OFF_filename_1 = readmatrix(filename_1, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % read data from filename_2
        data_OFF_filename_2 = readmatrix(filename_2, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % combine values into within_OFF array
        within_OFF = [within_OFF; data_OFF_filename_1; data_OFF_filename_2]; % makes within_OFF 24x1 double
    end
end

% read between ON values
for i = 1:length(sheetnames_3)
    % check for 'ON' in sheet names
    if contains(sheetnames_3{i}, 'ON')
        data_BETWEEN_ON = readmatrix(filename_3, 'Sheet', sheetnames_3{i}, 'Range', 'A2:A4');
        between_ON = [between_ON; data_BETWEEN_ON]; % makes between_ON 48x1 double
    end
end

% read between OFF values
for i = 1:length(sheetnames_3)
    % check for 'OFF' in sheet names
    if contains(sheetnames_3{i}, 'OFF')
        data_BETWEEN_OFF = readmatrix(filename_3, 'Sheet', sheetnames_3{i}, 'Range', 'A2:A4');
        between_OFF = [between_OFF; data_BETWEEN_OFF];  % makes between_OFF 48x1 double
    end
end

% determine maximum length and pad arrays with NaNs
max_length = max([length(within_ON), length(between_ON), length(within_OFF), length(between_OFF)]);
within_ON = [within_ON; nan(max_length - length(within_ON), 1)];
between_ON = [between_ON; nan(max_length - length(between_ON), 1)];
within_OFF = [within_OFF; nan(max_length - length(within_OFF), 1)];
between_OFF = [between_OFF; nan(max_length - length(between_OFF), 1)];

% combine data into a single cell array for boxplot
iccdata_mat = {within_ON, between_ON, within_OFF, between_OFF};
labels = {'W', 'B', 'W', 'B'};
group_index = [1 1 2 2];

% create the boxplot
figure('Position', [1, 160, 400, 600]);
daboxplot(iccdata_mat,'groups', group_index, 'colors', [0.5, 0.8, 1; 0.5, 0.8, 1; 0, 0.45, 0.7; 0, 0.45, 0.7],'whiskers', 0, 'scatter', 1, 'scattersize', 200, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5, 'outsymbol', 'kx', 'boxwidth', 2.5, 'boxspacing', 0.1);
set(gca, 'XTickLabel', labels, 'FontSize', 50); % set x-tick labels to your specified labels
ylabel('ICC','FontSize', 60, 'Color', 'k');
% yticks(0.3:0.1:0.8); % edit depending on axis size
% ylim([0.3 0.8]);

% run Wilxcoxon rank sum test with uneven samples
[p_ON,h_ON] = ranksum(within_ON, between_ON); 
fprintf('P-value: %.7e\n', p_ON); 
fprintf('H-value: %.0f\n', h_ON); 
[p_OFF,h_OFF] = ranksum(within_OFF, between_OFF);
fprintf('P-value: %.7e\n', p_OFF); 
fprintf('H-value: %.0f\n', h_OFF); 

% get the current axis limits for positioning the p-value text
ylimits = ylim;

% set the y position for the p-values slightly below the top of the plot
y_offset = 0.02 * (ylimits(2) - ylimits(1));  % offset for positioning

% display p-values on the plot just below ON and OFF labels
text(1.5, ylimits(2) - y_offset, sprintf('p = %.1e', p_ON), 'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', 'k');
text(3.5, ylimits(2) - y_offset, sprintf('p = %.1e', p_OFF), 'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', 'k');

% set other properties
set(gca, 'FontSize', 30); % increase axis font size
set(gcf, 'Color', 'w'); % optional: set figure background color

effect_ON = meanEffectSize(within_ON, between_ON, Effect="cliff", paired=false, ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(p_ON);
disp(effect_ON);

effect_OFF = meanEffectSize(within_OFF, between_OFF, Effect="cliff", paired=false, ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(p_OFF);
disp(effect_OFF);

%% Plot box plot for functional connectivity and/or brain variability of within and between ICC values of ON OFF DBS 
% x = [on, off (within)] [on, off (between)] 
% y = ICC

% define the Excel file name
filename_1 = 'V1_BV_ICC_Results.xlsx'; % within FC/BV ICC values visit 1
filename_2 = 'V2_BV_ICC_Results.xlsx'; % within FC/BV ICC values visit 2
filename_3 = 'V1_V2_BV_ICC_Between_Results.xlsx'; % between FC/BV ICC values visit 1 & 2

% display sheet names for each file:
sheets_1 = sheetnames('V1_BV_ICC_Results.xlsx'); % for file 1
disp('Sheet names:');
disp(sheets_1);
sheets_2 = sheetnames('V2_BV_ICC_Results.xlsx'); % for file 2
disp('Sheet names:');
disp(sheets_2);
sheets_3 = sheetnames('V1_V2_BV_ICC_Between_Results.xlsx'); % for file 3
disp('Sheet names:');
disp(sheets_3);

% list of sheet names
sheetnames_1_2 = {'Whole_DBS_OFF', 'Whole_DBS_ON', 'Motor_DBS_OFF', 'Motor_DBS_ON',...
    'Limbic_DBS_OFF', 'Limbic_DBS_ON', 'Associative_DBS_OFF', 'Associative_DBS_ON'}; % filename_1 & filemame_2
sheetnames_3 = {'whole_DBSOFF_v1rtv2rt', 'whole_DBSOFF_v1tv2t',	'whole_DBSOFF_v1rtv2t',...
    'whole_DBSOFF_v1tv2rt',	'whole_DBSON_v1rtv2rt',	'whole_DBSON_v1tv2t','whole_DBSON_v1rtv2t',...
    'whole_DBSON_v1tv2rt',	'motor_DBSOFF_v1rtv2rt', 'motor_DBSOFF_v1tv2t', 'motor_DBSOFF_v1rtv2t',...
    'motor_DBSOFF_v1tv2rt',	'motor_DBSON_v1rtv2rt',	'motor_DBSON_v1tv2t', 'motor_DBSON_v1rtv2t',...
    'motor_DBSON_v1tv2rt',	'limbic_DBSOFF_v1rtv2rt', 'limbic_DBSOFF_v1tv2t', 'limbic_DBSOFF_v1rtv2t',...
    'limbic_DBSOFF_v1tv2rt', 'limbic_DBSON_v1rtv2rt', 'limbic_DBSON_v1tv2t', 'limbic_DBSON_v1rtv2t',...
    'limbic_DBSON_v1tv2rt',	'associative_DBSOFF_v1rtv2rt', 'associative_DBSOFF_v1tv2t', 'associative_DBSOFF_v1rtv2t',...
    'associative_DBSOFF_v1tv2rt', 'associative_DBSON_v1rtv2rt',	'associative_DBSON_v1tv2t', 'associative_DBSON_v1rtv2t',...
    'associative_DBSON_v1tv2rt'}; % filename_3

% labels for sheet names
labels_1_2 = {'Whole DBS OFF', 'Whole DBS ON', 'Motor DBS OFF', 'Motor DBS ON',...
    'Limbic DBS OFF', 'Limbic DBS ON', 'Associative DBS OFF', 'Associative DBS ON'}; % filename_1 & filemame_2
labels_3 = {'Whole DBSOFF V1RTV2RT', 'Whole DBSOFF V1TV2T', 'Whole DBSOFF V1RTV2T', 'Whole DBSOFF V1TV2RT',... 
'Whole DBSON V1RTV2RT', 'Whole DBSON V1TV2T', 'Whole DBSON V1RTV2T', 'Whole DBSON V1TV2RT',... 
'Motor DBSOFF V1RTV2RT', 'Motor DBSOFF V1TV2T', 'Motor DBSOFF V1RTV2T', 'Motor DBSOFF V1TV2RT', ...
'Motor DBSON V1RTV2RT', 'Motor DBSON V1TV2T', 'Motor DBSON V1RTV2T', 'Motor DBSON V1TV2RT',... 
'Limbic DBSOFF V1RTV2RT', 'Limbic DBSOFF V1TV2T', 'Limbic DBSOFF V1RTV2T', 'Limbic DBSOFF V1TV2RT', ...
'Limbic DBSON V1RTV2RT', 'Limbic DBSON V1TV2T', 'Limbic DBSON V1RTV2T', 'Limbic DBSON V1TV2RT',... 
'Associative DBSOFF V1RTV2RT', 'Associative DBSOFF V1TV2T', 'Associative DBSOFF V1RTV2T',... 
'Associative DBSOFF V1TV2RT', 'Associative DBSON V1RTV2RT', 'Associative DBSON V1TV2T', ...
'Associative DBSON V1RTV2T', 'Associative DBSON V1TV2RT'}; % filename_3

% initialize arrays to hold the data
within_ON = [];
between_ON = [];
within_OFF = [];
between_OFF = [];
 
% read within ON values from both file filename_1 and filename_2 and combine
for i = 1:length(sheetnames_1_2)
    % check if the sheet corresponds to 'ON'
    if contains(sheetnames_1_2{i}, 'ON')
        % read data from filename_1
        data_ON_filename_1 = readmatrix(filename_1, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % read data from filename_2
        data_ON_filename_2 = readmatrix(filename_2, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % combine values into within_ON array
        within_ON = [within_ON; data_ON_filename_1; data_ON_filename_2]; % makes within_ON 24x1 double
    end
end

% read within OFF values from both file filename_1 and filename_2 and combine
for i = 1:length(sheetnames_1_2)
    % check if the sheet corresponds to 'OFF'
    if contains(sheetnames_1_2{i}, 'OFF')
        % read data from filename_1
        data_OFF_filename_1 = readmatrix(filename_1, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % read data from filename_2
        data_OFF_filename_2 = readmatrix(filename_2, 'Sheet', sheetnames_1_2{i}, 'Range', 'A2:A4');
        % combine values into within_OFF array
        within_OFF = [within_OFF; data_OFF_filename_1; data_OFF_filename_2]; % makes within_OFF 24x1 double
    end
end

% read between ON values
for i = 1:length(sheetnames_3)
    % check for 'ON' in sheet names
    if contains(sheetnames_3{i}, 'ON')
        data_BETWEEN_ON = readmatrix(filename_3, 'Sheet', sheetnames_3{i}, 'Range', 'A2:A4');
        between_ON = [between_ON; data_BETWEEN_ON]; % makes between_ON 48x1 double
    end
end

% read between OFF values
for i = 1:length(sheetnames_3)
    % check for 'OFF' in sheet names
    if contains(sheetnames_3{i}, 'OFF')
        data_BETWEEN_OFF = readmatrix(filename_3, 'Sheet', sheetnames_3{i}, 'Range', 'A2:A4');
        between_OFF = [between_OFF; data_BETWEEN_OFF];  % makes between_OFF 48x1 double
    end
end

% determine maximum length and pad arrays with NaNs
max_length = max([length(within_ON), length(between_ON), length(within_OFF), length(between_OFF)]);
within_ON = [within_ON; nan(max_length - length(within_ON), 1)];
between_ON = [between_ON; nan(max_length - length(between_ON), 1)];
within_OFF = [within_OFF; nan(max_length - length(within_OFF), 1)];
between_OFF = [between_OFF; nan(max_length - length(between_OFF), 1)];

% combine data into a single cell array for boxplot
iccdata_mat = {within_ON, within_OFF, between_ON, between_OFF};
labels = {'ON', 'OFF', 'ON', 'OFF'};
group_index = [1 1 2 2];

% create the boxplot
figure('Position', [1, 160, 400, 600]);
daboxplot(iccdata_mat,'groups', group_index, 'colors', [0.5, 0.8, 1; 0, 0.45, 0.7; 0.5, 0.8, 1; 0, 0.45, 0.7],'whiskers', 0, 'scatter', 1, 'scattersize', 200, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5, 'outsymbol', 'kx', 'boxwidth', 2.5, 'boxspacing', 0.1);
set(gca, 'XTickLabel', labels, 'FontSize', 50); % set x-tick labels to your specified labels
ylabel('ICC','FontSize', 50, 'Color', 'k');
% yticks(0.3:0.1:0.8); % edit depending on axis size
% ylim([0.3 0.8]);

% run Wilxcoxon signed rank test with uneven samples
[p_within,h_within] = signrank(within_ON, within_OFF); 
fprintf('P-value: %.7e\n', p_within); 
fprintf('H-value: %.0f\n', h_within); 
[p_between,h_between] = signrank(between_ON, between_OFF);
fprintf('P-value: %.7e\n', p_between); 
fprintf('H-value: %.0f\n', h_between); 

% get the current axis limits for positioning the p-value text
ylimits = ylim;

% set the y position for the p-values slightly below the top of the plot
y_offset = 0.02 * (ylimits(2) - ylimits(1));  % offset for positioning

% display p-values on the plot just below within and between labels
text(1.5, ylimits(2) - y_offset, sprintf('p = %.1e', p_within), 'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', 'k');
text(3.5, ylimits(2) - y_offset, sprintf('p = %.1e', p_between), 'HorizontalAlignment', 'center', 'FontSize', 25, 'Color', 'k');

% set other properties
set(gca, 'FontSize', 30); % Increase axis font size
set(gcf, 'Color', 'w'); % Optional: set figure background color

effect_within = meanEffectSize(within_ON, within_OFF, Effect="cliff", paired=true, ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(p_within);
disp(effect_within);

effect_between = meanEffectSize(between_ON, between_OFF, Effect="cliff", paired=true, ConfidenceIntervalType="bootstrap", ...
    BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(p_between);
disp(effect_between);

%% Figure 4B - Connectomes of Patient in longitudinal cohort with outlier ICC (visual inspection)

pt_ID_outliermotorICC = 3; % Patient ID in longitudinal cohort

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

for i=1:4 
    figure('Position',[0 0 900 750]);
    set(gcf,'color','w');
    fc=squeeze(fcs_motor_2(pt_ID_outliermotorICC,:,:,:)); %avg across selected patients
    imagesc(squeeze(fc(i,:,:))); %plots only one run at a time
    set(gca,'XTickLabel',corrname,'XTickLabelRotation',50);
    set(gca,'YTickLabel',corrname);
    set(gca, 'XTick', (1:length(corrname)), 'FontSize', 12);
    set(gca, 'YTick', (1:length(corrname)), 'FontSize',12);
    grid off
    box off
    colormap("parula");
     if i==1
    title('DBS OFF (Test)' ,'FontSize', 25);
    elseif i==2
    title('DBS OFF (Retest)' ,'FontSize', 25);
    elseif i==3
    title('DBS ON (Test)' ,'FontSize', 25);
    elseif i==4
    title('DBS ON (Retest)' ,'FontSize', 25);
    end
    a = colorbar; caxis([-1, 1]);
    ylabel(a,'Functional Connectivity','FontSize',25,'Rotation',90);

end