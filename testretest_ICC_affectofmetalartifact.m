%% Calculating ICC Values for Functional Connectivity and Brain Variability in Affected 
%% vs. Unaffected Brain ROIs by Creating a Map of Electrode Spatial Trajectory, 
%% and Analysis of Affected vs. Unaffected ICCs (Figure 5)

% This script identifies which ROIs are affected and unaffected by the DBS
% leads, calculates respective ICCs (affected vs. unaffected), and corresponding analyses.
% This script is meant to be run section by section in order and used AFTER running the
% testretest_ICC.m script.

%% ------------------------------------------------------------------------------------- %

%PART 1: Create a 3D reconstruction of the electrodes using paraview [OUTSIDE OF MATLAB]

%% ------------------------------------------------------------------------------------- %%

% First install ITK Snap and Paraview softwares
% The file glassbrain.vtk is the base brain image in MNI atlas space.
% The overlay files is the mask you created which needs to be converted to VTK format before display. Note: Both these files are masks. 
% ITK Snap is used to convert the binary electrode masks into 3D surface mesh in the VTK format.
% In ITKSnap, File> Open> electrode_mask.nii  > Next > Finish
% Segmentation > Open segmentation > electrode_mask.nii
% Segmentation > Export as surface mesh > Export meshes for all labels as a single scene > Save as VTKPolyData file
% Open Paraview. File > Open > glassbrain.vtk
% Again File > Open > electrode_mask.vtk
% In the pipeline browser on the left side, click on the eye symbol next to each file to display the image. 
% For the glass brain, change the opacity and color. 
% Now you can see the ROIs overlaid on glass brain.

%% ------------------------------------------------------------------------------------- %

%PART 2: Identify ROIs that are affected and unaffected by the DBS lead for
%each patient

%% ------------------------------------------------------------------------------------- %%
% Repeat this part until you have run it for each patient using both the qsm and the 
% qsm_dgm atlases. 

% At the end, you will have two output text files containing the affected
% ROIs for all patients in the R and L hemispheres, respectively.

%% Set paths (modify to where the data is on your computer and change for each individual patient)


addpath('/LOCAL/PATH/HERE/spm12.local/'); % path to SPM 12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
addpath('/LOCAL/PATH/HERE/matlab_snippets/handle_nifti/nifti_tools/'); % path to the "Tools for NIfTI and ANALYZE image" function (Jimmy Shen (2025). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange.)

mnirois = '/YOUR/PATH/HERE/'; % path to where atlas NIfTI files and text files are stored
mnirois_txt = 'qsm.txt'; % qsm atlas text file - change to 'qsm_dgm.txt' as needed
mnirois_nii = 'qsm.nii'; % qsm atlas NIfTI file - change to 'qsm_dgm.nii' as needed
pathmnirois_nii = [mnirois mnirois_nii]; 
pathmnirois_txt = [mnirois mnirois_txt];  

dir_data = '/YOUR/PATH/HERE/'; % path to where the output patient CONN files are stored.
pt_id = 'PT_ID_HERE'; % the individual patient's ID

dir_patient = [dir_data pt_id]; % path to individual patient CONN files and where output trajectory NIfTI files will be saved

dir_mask = ['/YOUR/PATH/HERE/' pt_id]; % path to where patient level binary electrode masks are stored
mask_name = 'qsm_lead-seg-otsu_cleaned.nii'; % name of binary electrode mask

%% Align (reslice) masks to atlases
% This step allows you to view the mask overlayed with the patient's
% fMRI NiFTI file.

% load reference NIfTI (the one with the desired voxel size & dimensions)
refNii = spm_vol(pathmnirois_nii);

% load the target NIfTI (the one that needs resampling)
targetNii = spm_vol(dir_mask);

% configure reslicing
flags = struct('interp', 1, 'mask', 0, 'mean', 0, 'which', 1, 'wrap', [0 0 0], 'prefix', 'qsm_');

% reslice the target file to match the atlas
spm_reslice({refNii.fname, targetNii.fname}, flags);

%% Load in mnirois brain atlas in MNI space 

FileStruct = load_untouch_nii(pathmnirois_nii); 
rois = FileStruct.img; % original mask
imagesc(rois(:,:,90)); colormap (colorcube); % test slice
size(rois)

figure;
set(gca,'xtick',[], 'ytick',[]);
set(gcf,'color','white')
montage(rois, 'DisplayRange', [min(rois(:)), max(rois(:))]);
colormap(colorcube);

%% Multiple electrode mask by brain atlas 

cd(dir_mask)

FileStruct=load_untouch_nii(mask_name); % load lead mask
electrodes=FileStruct.img;

roi_dims=size(rois);

% calculate the padding needed for each dimension
pad_x_pre = floor((roi_dims(1) - size(electrodes, 1)) / 2);
pad_x_post = ceil((roi_dims(1) - size(electrodes, 1)) / 2);
pad_y_pre = floor((roi_dims(2) - size(electrodes, 2)) / 2);
pad_y_post = ceil((roi_dims(2) - size(electrodes, 2)) / 2);
pad_z_pre = floor((roi_dims(3) - size(electrodes, 3)) / 2);
pad_z_post = ceil((roi_dims(3) - size(electrodes, 3)) / 2);

% pad the electrode array to match mni_rois dimensions
electrodes_padded = padarray(electrodes, [pad_x_pre, pad_y_pre, pad_z_pre], 0, 'pre');
electrodes_padded = padarray(electrodes_padded, [pad_x_post, pad_y_post, pad_z_post], 0, 'post');

% voxel size check
info_electrodes = niftiinfo(mask_name);
info_mnirois = niftiinfo(pathmnirois_nii);

if isequal(info_electrodes.PixelDimensions, info_mnirois.PixelDimensions)
    disp('Voxel sizes match.');
else
    disp('Voxel sizes do not match.');
end

% spatial orientation check
disp('Affine transformation matrix for electrode(s):');
disp(info_electrodes.Transform.T); % equivalent to ImageOrientationPatient and ImagePositionPatient
disp('Affine transformation matrix for mni_rois:');
disp(info_mnirois.Transform.T); % equivalent to ImageOrientationPatient and ImagePositionPatient

% multiply electrode by atlas (to quantify)
electrodes_labeled = (electrodes_padded).*rois;

% show electrode over atlas (visual only)
electrodes_labeled_2 = uint16(1000*electrodes_padded)+rois; 
montage(electrodes_labeled_2(:,:,65:120)); colormap (colorcube); 

%% Identify brain regions electrode and IPG passes through

cd(pathmnirois_txt); % text files for all atlases stored here

values = readtable(mnirois_txt, 'Delimiter', '\n', 'ReadVariableNames', false); % change for each mniroi 
disp(values);

% analyze electrodes
histogram(nonzeros(electrodes_labeled)) % excluding zeros, plot all the values
set(gcf,'color','white')
xlabel('brain atlas label number')
ylabel('count')

k_electrodes = unique(nonzeros(electrodes_labeled));
fprintf('This patient DBS lead is passing through brain regions: %s\n', strjoin(values.Var1(k_electrodes), ', '));

cd(dir_data);

% output files
right_output_file = 'right_electrode_brain_regions_output_new.txt';
left_output_file = 'left_electrode_brain_regions_output_new.txt';

% append regions for right electrode
if ~isempty(k_electrodes)
    existing_regions_R = string.empty; 
    if exist(right_output_file, 'file')
        existing_content = readcell(right_output_file, 'Delimiter', '\n');
        
        % initialize an empty array for existing regions
        for line = existing_content'
            % check if the line contains the expected format
            if contains(line{1}, 'Region ') && contains(line{1}, ':')
                region = extractBetween(line{1}, 'Region ', ':');
                existing_regions_R = [existing_regions_R; string(strtrim(region))];  
            end
        end
    end
    
    % open file for appending
    fid_R = fopen(right_output_file, 'a');
    if fid_R ~= -1
        fprintf(fid_R, 'Patient: %s\n', pt_id);  % replace with actual patient ID
        for i = 1:length(k_electrodes)
            region = string(values.Var1(k_electrodes(i)));
            if contains(region,  ' r', 'IgnoreCase', true)
                if ~ismember(region, existing_regions_R)
                    fprintf(fid_R, 'Region %d: %s\n', k_electrodes(i), region);
                end
            end
        end
        fprintf(fid_R, '\n');
        fclose(fid_R);
        fprintf('Appended output to %s\n', right_output_file);
    else
        error('Could not open file for writing: %s', right_output_file);
    end
else
    disp('No valid indices found in right.');
end

% append regions for left electrode
if ~isempty(k_electrodes)
    existing_regions_L = string.empty; 
    if exist(left_output_file, 'file')
        existing_content = readcell(left_output_file, 'Delimiter', '\n');
        
        % initialize an empty array for existing regions
        for line = existing_content'
            % check if the line contains the expected format
            if contains(line{1}, 'Region ') && contains(line{1}, ':')
                region = extractBetween(line{1}, 'Region ', ':');
                existing_regions_L = [existing_regions_L; string(strtrim(region))];  
            end
        end
    end
    
    % open file for appending
    fid_L = fopen(left_output_file, 'a');
    if fid_L ~= -1
        fprintf(fid_L, 'Patient: %s\n', pt_id);
        for i = 1:length(k_electrodes)
            region = string(values.Var1(k_electrodes(i)));
            if contains(region,  ' l', 'IgnoreCase', true)
                if ~ismember(region, existing_regions_L)
                    fprintf(fid_L, 'Region %d: %s\n', k_electrodes(i), region);
                end
            end
        end
        fprintf(fid_L, '\n');
        fclose(fid_L);
        fprintf('Appended output to %s\n', left_output_file);
    else
        error('Could not open file for writing: %s', left_output_file);
    end
else
    disp('No valid indices found in left.');
end

%% Visualize brain regions electrode and IPG passes through

% electrodes
mni_rois_trajectory = zeros(size(rois)); 

% map both side brain regions
for i = 1:length(k_electrodes) 
    % for every brain region in the lead trajectory
    mni_rois_trajectory(rois == k_electrodes(i)) = k_electrodes(i); 
end

% set non-passing regions to zero
mni_rois_trajectory(mni_rois_trajectory == 0) = 0; 

% display unique values to confirm both side regions
unique_regions = unique(mni_rois_trajectory);
disp('Unique regions in electrode trajectories:');
disp(unique_regions); % should only show the brain regions and 0

figure;
cmap = lines(101); % use 'lines' to generate distinct colors
cmap(1, :) = [0 0 0]; % set color for index 0
for idx = 2:101
cmap(idx, :) = rand(1, 3); % random color for each index
end
imagesc(mni_rois_trajectory(:,:,76)); colormap (cmap);

figure;
montage(mni_rois_trajectory*5000); colormap colorcube; 

cd (dir_patient)

% save the entire trajectory as a NIfTI file
FileStruct.img = mni_rois_trajectory; 
FileStruct.hdr.dime.dim = [3, size(FileStruct.img, 1), size(FileStruct.img, 2), size(FileStruct.img, 3), 1, 1, 1, 1]; 

cd(dir_patient);

save_untouch_nii(FileStruct, 'mni_rois_trajectory_new_qsm.nii.gz');
disp('Saved mni_rois_trajectory.... as NIfTI file.'); 


%% ------------------------------------------------------------------------ %%

          %PART 3: ICC Calculation - Effect of Metal Artifact

%% ------------------------------------------------------------------------ %%

% IMPORTANT: Part 2 and part 3 do NOT run in conjunction. You must finish
% running Part 2 for ALL patients of interest before moving on to part 3, as part 3
% analysis requires output from all patients using BOTH the qsm and qsm_dgm
% atlases.

%% Add data paths
lead_seg_data = 'YOUR/PATH/HERE'; % path to where outputs text files from 'testretest_dbs_spatialtrajectorymap.m' are stored.

left_electrode_ROIs = 'left_electrode_brain_regions_output_new.txt'; % 
right_electrode_ROIs = 'right_electrode_brain_regions_output_new.txt'; 

% Add paths to functions.
addpath('YOUR/PATH/HERE/ICC.m'); % Add a path to where the ICC function is stored. (Arash Salarian (2025). Intraclass Correlation Coefficient (ICC) (https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc), MATLAB Central File Exchange.)
addpath('YOUR/PATH/HERE/daboxplot/'); % Add a path to where the daboxplot function is stored. (Povilas Karvelis (2025). daboxplot (https://github.com/frank-pk/DataViz/releases/tag/v3.2.3), GitHub.)

patientIDs = {'LIST ALL PT IDs HERE'}; % A list of all the patient IDs

dataset_name = 'YOUR_DATASET_NAME_HERE'; % The name of the dataset which will be used to name the output Excel file.

%% Load Lead Segmentation Data
% Read in output data indicating ROIs are affected for each patient.
cd(lead_seg_data)

left_electrode = fileread(left_electrode_ROIs);
right_electrode = fileread(right_electrode_ROIs);

% Regular expression to extract all patients with their conditions (ROIs)
patientPattern = 'Patient: (\S+) \(([^)]+)\)';
patients = {regexp(left_electrode, patientPattern, 'tokens'), ...
            regexp(right_electrode, patientPattern, 'tokens')};

% Regular expression to extract regions (those starting with "Region")
regionPattern = 'Region \d+: ([^()]+)';

% Initialize variables to store patients and their regions
affected_regions = {cell(numel(patients{1}), 3), cell(numel(patients{2}), 3)};

% Loop through both left and right electrodes
for side = 1:2
    electrode_data = {left_electrode, right_electrode};
    startIdx = 1;
    
    for i = 1:numel(patients{side})
        if i < numel(patients{side})
            nextStartIdx = regexp(electrode_data{side}(startIdx:end), ['Patient: ' patients{side}{i+1}{1}], 'once');
            if ~isempty(nextStartIdx)
                endIdx = startIdx + nextStartIdx - 1;
            else
                endIdx = numel(electrode_data{side});
            end
        else
            endIdx = numel(electrode_data{side});
        end
        
        currentSection = electrode_data{side}(startIdx:endIdx);
        
        % Store patient ID and brain target
        affected_regions{side}{i, 1} = patients{side}{i}{1};
        affected_regions{side}{i, 2} = patients{side}{i}{2};
        
        % Extract and store the affected ROIs
        patientRegions = regexp(currentSection, regionPattern, 'tokens');
        patientRegions = strtrim(patientRegions);
        % Flatten the nested cell structure and store the regions
        affected_regions{side}{i, 3} = [patientRegions{:}]; 

        
        % Update start index
        startIdx = endIdx + 1;
    end
end

% Display results and check that data was extracted properly.
disp('Affected Regions for Left Electrode:');
disp(affected_regions{1});

disp('Affected Regions for Right Electrode:');
disp(affected_regions{2});

% Initialize the matrix to store combined regions
combined_regions = cell(numel(patientIDs), 2);

% Loop through each patient ID
for p = 1:numel(patientIDs)
    patientID = patientIDs{p};
    combinedROIs = {};  % Initialize empty cell for combined ROIs
    
    % Check both left and right affected regions
    for side = 1:2
        for i = 1:size(affected_regions{side}, 1)
            if strcmp(affected_regions{side}{i, 1}, patientID)
                % Add the regions to the combined list
                combinedROIs = [combinedROIs; affected_regions{side}{i, 3}(:)];
            end
        end
    end
    
    % Remove duplicate regions
    combinedROIs = unique(combinedROIs);
    
    % Store the result in the new matrix
    combined_regions{p, 1} = patientID;
    combined_regions{p, 2} = combinedROIs;
end

% Display the combined regions matrix and check that ROIs are correct
disp(combined_regions);

corrname_whole_unaffected = cell(size(combined_regions, 1), 2);

for i = 1:size(combined_regions, 1)
    patientID = combined_regions{i, 1};
    affectedROIs = combined_regions{i, 2};

    % Replace affected ROIs with blank spaces while preserving order
    updatedROIs = corrname_whole;  % Start with the full list

    % Loop through and replace affected ROIs with blank spaces
    for j = 1:length(updatedROIs)
        if ismember(updatedROIs{j}, affectedROIs)
            updatedROIs{j} = '';  % Replace with blank space
        end
    end

    % Store the results
    corrname_whole_unaffected{i, 1} = patientID;
    corrname_whole_unaffected{i, 2} = updatedROIs;
end

% Check that affected regions were properly removed.
disp(corrname_whole_unaffected);

% Create indexes for unaffected regions and affected regions
unaffected_idx = cell(size(corrname_whole_unaffected, 1), 1); % to store the unaffected indices for each patient
affected_idx = cell(size(corrname_whole_unaffected, 1),1);

for i = 1:size(corrname_whole_unaffected, 1)
    % Get the list of unaffected ROIs for this patient
    unaffectedROIs = corrname_whole_unaffected{i, 2}; 
    
    % Create an empty vector to store indices of unaffected regions
    patient_unaffected_index = [];
    patient_affected_index = [];
    
    % Loop through each of the 142 Whole Brain ROIs
    for j = 1:142
        if ~isempty(unaffectedROIs{j}) % Check if the ROI is unaffected (non-blank)
            patient_unaffected_index = [patient_unaffected_index, j]; % Add the index of unaffected ROIs
        else
            patient_affected_index = [patient_affected_index, j]; % Add the index of affected ROIs
        end
    end
    
    % Store the unaffected indices for the current patient
    unaffected_idx{i} = patient_unaffected_index;
    affected_idx{i} = patient_affected_index;
    
end

%% calculate ICC values for functional connectivity and brain variability.

% Initialize matrices to store ICC results for each connectome
% FC = functional connectivity
% BV = brain variability
% 7 = 7 ICC output variables
FC_ICC_DBSOFF_whole_unaffected = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
FC_ICC_DBSON_whole_unaffected = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON
BV_ICC_DBSOFF_whole_unaffected = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
BV_ICC_DBSON_whole_unaffected = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON

FC_ICC_DBSOFF_whole_affected = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
FC_ICC_DBSON_whole_affected = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON
BV_ICC_DBSOFF_whole_affected = zeros(num_subj, 7); % ICC results for 'whole' connectome, DBS OFF
BV_ICC_DBSON_whole_affected = zeros(num_subj, 7);  % ICC results for 'whole' connectome, DBS ON

% Arrange rows (measures) and columns (rater or test-retest)  & perform ICC
for i = 1:num_subj

    % Get the unaffected indices for the current patient
    unaffected_idx_i = unaffected_idx(i);
    affected_idx_i = affected_idx(i);

    % Ensure unaffected_indices is numeric (convert if necessary)
    unaffected_idx_i = cell2mat(unaffected_idx_i); % Convert cell to numeric array if needed
    affected_idx_i = cell2mat(affected_idx_i);

    % Functional connectivity - whole brain unaffected
    M_OFF = [nonzeros(fcs_whole_tril(i, 1, unaffected_idx_i, :)) nonzeros(fcs_whole_tril(i, 2, unaffected_idx_i, :))];
    disp(M_OFF)
    M_ON = [nonzeros(fcs_whole_tril(i, 3, unaffected_idx_i, :)) nonzeros(fcs_whole_tril(i, 4, unaffected_idx_i, :))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0); %calc ICC
    FC_ICC_DBSOFF_whole_unaffected(i, :) = [r_off, LB, UB, F, df1, df2, p]; %assign to results vector
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    FC_ICC_DBSON_whole_unaffected(i, :) = [r_on, LB, UB, F, df1, df2, p];

    % Brain variability - whole brain unaffected
    M_OFF = [nonzeros(brain_var_whole(i, 1, unaffected_idx_i)) nonzeros(brain_var_whole(i, 2,unaffected_idx_i))];
    disp(M_OFF)
    M_ON = [nonzeros(brain_var_whole(i, 3, unaffected_idx_i)) nonzeros(brain_var_whole(i, 4, unaffected_idx_i))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    BV_ICC_DBSOFF_whole_unaffected(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    BV_ICC_DBSON_whole_unaffected(i, :) = [r_on, LB, UB, F, df1, df2, p];

    % Functional connectivity - whole brain affected
    M_OFF = [nonzeros(fcs_whole_tril(i, 1, affected_idx_i,:)) nonzeros(fcs_whole_tril(i, 2, affected_idx_i, :))];
    disp(M_OFF)
    M_ON = [nonzeros(fcs_whole_tril(i, 3, affected_idx_i, :)) nonzeros(fcs_whole_tril(i, 4, affected_idx_i, :))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0); %calc ICC
    FC_ICC_DBSOFF_whole_affected(i, :) = [r_off, LB, UB, F, df1, df2, p]; %assign to results vector
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    FC_ICC_DBSON_whole_affected(i, :) = [r_on, LB, UB, F, df1, df2, p];

    % Brain variability - whole brain affected
    M_OFF = [nonzeros(brain_var_whole(i, 1, affected_idx_i)) nonzeros(brain_var_whole(i, 2,affected_idx_i))];
    disp(M_OFF)
    M_ON = [nonzeros(brain_var_whole(i, 3, affected_idx_i)) nonzeros(brain_var_whole(i, 4, affected_idx_i))];
    disp(M_ON)
    [r_off, LB, UB, F, df1, df2, p] = ICC(M_OFF, ICC_type, alpha, r0);
    BV_ICC_DBSOFF_whole_affected(i, :) = [r_off, LB, UB, F, df1, df2, p];
    [r_on, LB, UB, F, df1, df2, p] = ICC(M_ON, ICC_type, alpha, r0);
    BV_ICC_DBSON_whole_affected(i, :) = [r_on, LB, UB, F, df1, df2, p];
end

%% Save ICC Values to Excel

% Functional Connectivity ICC Values

% Define header for ICC table
fc_icc_header = {'FC ICC', 'Lower Bound', 'Upper Bound', 'F', 'df1', 'df2', 'p-value'};

% Create tables for each connectome type and DBS condition
table_fc_whole_dbs_off_affected = array2table(FC_ICC_DBSOFF_whole_affected, 'VariableNames', fc_icc_header);
table_fc_whole_dbs_on_affected = array2table(FC_ICC_DBSON_whole_affected, 'VariableNames', fc_icc_header);

table_fc_whole_dbs_off_unaffected = array2table(FC_ICC_DBSOFF_whole_unaffected, 'VariableNames', fc_icc_header);
table_fc_whole_dbs_on_unaffected = array2table(FC_ICC_DBSON_whole_unaffected, 'VariableNames', fc_icc_header);


% Define the filename for the Excel file
filename = [dataset_name '_FC_ICC_Results_affectedunaffected.xlsx'];

% Write each table to a separate sheet in the Excel file
writetable(table_fc_whole_dbs_off_affected, filename, 'Sheet', 'Whole_DBS_OFF_affected');
writetable(table_fc_whole_dbs_on_affected, filename, 'Sheet', 'Whole_DBS_ON_affected');

writetable(table_fc_whole_dbs_off_unaffected, filename, 'Sheet', 'Whole_DBS_OFF_unaffected');
writetable(table_fc_whole_dbs_on_unaffected, filename, 'Sheet', 'Whole_DBS_ON_unaffected');

%Brain Variability

% Define header for ICC table
bv_icc_header = {'BV ICC', 'Lower Bound', 'Upper Bound', 'F', 'df1', 'df2', 'p-value'};

% Create tables for each connectome type and DBS condition
table_bv_whole_dbs_off_affected = array2table(BV_ICC_DBSOFF_whole_affected, 'VariableNames', bv_icc_header);
table_bv_whole_dbs_on_affected = array2table(BV_ICC_DBSON_whole_affected, 'VariableNames', bv_icc_header);

table_bv_whole_dbs_off_unaffected = array2table(BV_ICC_DBSOFF_whole_unaffected, 'VariableNames', bv_icc_header);
table_bv_whole_dbs_on_unaffected = array2table(BV_ICC_DBSON_whole_unaffected, 'VariableNames', bv_icc_header);

% Define the filename for the Excel file
filename = [dataset_name '_BV_ICC_Results_affectedunaffected.xlsx'];

% Write each table to a separate sheet in the Excel file
writetable(table_bv_whole_dbs_off_affected, filename, 'Sheet', 'Whole_DBS_OFF_affected');
writetable(table_bv_whole_dbs_on_affected, filename, 'Sheet', 'Whole_DBS_ON_affected');

writetable(table_bv_whole_dbs_off_unaffected, filename, 'Sheet', 'Whole_DBS_OFF_unaffected');
writetable(table_bv_whole_dbs_on_unaffected, filename, 'Sheet', 'Whole_DBS_ON_unaffected');


%% ------------------------------------------------------------------------ %%

                   % PART 4: Analysis - Figure 5

%% ------------------------------------------------------------------------ %%

%% Figure 5C - Box Plot of Affected vs. Unaffected vs Whole Brain ICC values 

ICC_values = 'YOUR/PATH/HERE/'; % path to where ICC values from the testretest_ICC.m script are stored.
FC_ICC_file = 'FC_ICC_Results.xlsx'; % name of the FC ICC output file from the testretest_ICC.m script.
BV_ICC_file = 'BV_ICC_Results.xlsx'; % % name of the BV ICC output file from the testretest_ICC.m script. 

FC_ICC_DBSON_whole = readtable("FC_ICC_Results.xlsx", 'Sheet', 'Whole_DBS_ON');
FC_ICC_DBSOFF_whole = readtable("FC_ICC_Results.xlsx", 'Sheet', 'Whole_DBS_OFF');
BV_ICC_DBSON_whole = readtable("BV_ICC_Results.xlsx", 'Sheet', 'Whole_DBS_ON');
BV_ICC_DBSOFF_whole = readtable("BV_ICC_Results.xlsx", 'Sheet', 'Whole_DBS_OFF');

FC_ICC_DBSON_r_whole = FC_ICC_DBSON_whole(:,1);
FC_ICC_DBSOFF_r_whole = FC_ICC_DBSOFF_whole(:,1);
BV_ICC_DBSON_r_whole = BV_ICC_DBSON_whole(:,1);
BV_ICC_DBSOFF_r_whole = BV_ICC_DBSOFF_whole(:,1);

FC_ICC_DBSON_r_whole_unaffected = FC_ICC_DBSON_whole_unaffected(:,1);
FC_ICC_DBSOFF_r_whole_unaffected = FC_ICC_DBSOFF_whole_unaffected(:,1);
BV_ICC_DBSON_r_whole_unaffected = BV_ICC_DBSON_whole_unaffected(:,1);
BV_ICC_DBSOFF_r_whole_unaffected = BV_ICC_DBSOFF_whole_unaffected(:,1);

FC_ICC_DBSON_r_whole_affected = FC_ICC_DBSON_whole_affected(:,1);
FC_ICC_DBSOFF_r_whole_affected = FC_ICC_DBSOFF_whole_affected(:,1);
BV_ICC_DBSON_r_whole_affected = BV_ICC_DBSON_whole_affected(:,1);
BV_ICC_DBSOFF_r_whole_affected = BV_ICC_DBSOFF_whole_affected(:,1);

M_FC_ICC_dbsoff_leadanalysis =zeros(num_subj,3,1); % create empty reshaped matrices for functional connectivity ICC values.
M_FC_ICC_dbson_leadanalysis =zeros(num_subj,3,1);
M_BV_ICC_dbsoff_leadanalysis =zeros(num_subj,3,1); % create empty reshaped matrices for brain variability ICC values.
M_BV_ICC_dbson_leadanalysis=zeros(num_subj,3,1);

categories_leadanalysis = {'WB', 'U', 'A'}; % WB, U, A = whole brain, unaffected, affected

M_FC_ICC_dbsoff_leadanalysis(:, 1, 1) = FC_ICC_DBSOFF_r_whole;
M_FC_ICC_dbsoff_leadanalysis(:, 2, 1) = FC_ICC_DBSOFF_r_whole_unaffected;
M_FC_ICC_dbsoff_leadanalysis(:, 3, 1) = FC_ICC_DBSOFF_r_whole_affected;

M_FC_ICC_dbson_leadanalysis(:, 1, 1) = FC_ICC_DBSON_r_whole;
M_FC_ICC_dbson_leadanalysis(:, 2, 1) = FC_ICC_DBSON_r_whole_unaffected;
M_FC_ICC_dbson_leadanalysis(:, 3, 1) = FC_ICC_DBSON_r_whole_affected;

M_BV_ICC_dbsoff_leadanalysis(:, 1, 1) = BV_ICC_DBSOFF_r_whole;
M_BV_ICC_dbsoff_leadanalysis(:, 2, 1) = BV_ICC_DBSOFF_r_whole_unaffected;
M_BV_ICC_dbsoff_leadanalysis(:, 3, 1) = BV_ICC_DBSOFF_r_whole_affected;

M_BV_ICC_dbson_leadanalysis(:, 1, 1) = BV_ICC_DBSON_r_whole;
M_BV_ICC_dbson_leadanalysis(:, 2, 1) = BV_ICC_DBSON_r_whole_unaffected;
M_BV_ICC_dbson_leadanalysis(:, 3, 1) = BV_ICC_DBSON_r_whole_affected;


on_color = [0.5, 0.8, 1];
off_color = [0, 0.45, 0.7];

% FC DBS ON

figure('Color', 'w','Position', [100, 100, 300, 500]);
daboxplot(M_FC_ICC_dbson_leadanalysis,'color', on_color,'whiskers' , 0, 'scatter', 1, 'scattersize', 150, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories_leadanalysis,'fill',1,'boxspacing',0.1, 'boxwidth', 2);
ylabel('ICC', 'FontSize',25);
ylim([0.2,0.9]);
yticks(0.2:0.1:0.9);
xlim([0.5, 4.5]);
set(gca, 'XTickLabel', categories_leadanalysis, 'FontSize', 25, 'XTickLabelRotation', 0);

[p_au, h_au] = signrank(squeeze(M_FC_ICC_dbson_leadanalysis(:,2)), squeeze(M_FC_ICC_dbson_leadanalysis(:,3)));
[p_wu, h_wu] = signrank(squeeze(M_FC_ICC_dbson_leadanalysis(:,1)), squeeze(M_FC_ICC_dbson_leadanalysis(:,2)));
[p_wa, h_wa] = signrank(squeeze(M_FC_ICC_dbson_leadanalysis(:,1)), squeeze(M_FC_ICC_dbson_leadanalysis(:,3)));

effect_au = meanEffectSize(squeeze(M_FC_ICC_dbson_leadanalysis(:,2)), squeeze(M_FC_ICC_dbson_leadanalysis(:,3)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wu = meanEffectSize(squeeze(M_FC_ICC_dbson_leadanalysis(:,1)), squeeze(M_FC_ICC_dbson_leadanalysis(:,2)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wa = meanEffectSize(squeeze(M_FC_ICC_dbson_leadanalysis(:,1)), squeeze(M_FC_ICC_dbson_leadanalysis(:,3)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);

disp([p_au, p_wu, p_wa])
disp(effect_au)
disp(effect_wu)
disp(effect_wa)

% FC DBS OFF
figure('Color', 'w','Position', [100, 100, 300, 500]);
daboxplot(M_FC_ICC_dbsoff_leadanalysis,'color', off_color,'whiskers' , 0, 'scatter', 1, 'scattersize', 150, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories_leadanalysis,'fill',1,'boxspacing',0.1, 'boxwidth', 2);
ylabel('ICC', 'FontSize',25);
ylim([0.2,0.9]);
yticks(0.2:0.1:0.9);
xlim([0.5, 4.5]);
set(gca, 'XTickLabel', categories_leadanalysis, 'FontSize', 25, 'XTickLabelRotation', 0);

[p_au, h_au] = signrank(squeeze(M_FC_ICC_dbsoff_leadanalysis(:,2)), squeeze(M_FC_ICC_dbsoff_leadanalysis(:,3)));
[p_wu, h_wu] = signrank(squeeze(M_FC_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_FC_ICC_dbsoff_leadanalysis(:,2)));
[p_wa, h_wa] = signrank(squeeze(M_FC_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_FC_ICC_dbsoff_leadanalysis(:,3)));
disp([p_au, p_wu, p_wa])

effect_au = meanEffectSize(squeeze(M_FC_ICC_dbsoff_leadanalysis(:,2)), squeeze(M_FC_ICC_dbsoff_leadanalysis(:,3)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wu = meanEffectSize(squeeze(M_FC_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_FC_ICC_dbsoff_leadanalysis(:,2)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wa = meanEffectSize(squeeze(M_FC_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_FC_ICC_dbsoff_leadanalysis(:,3)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(effect_au)
disp(effect_wu)
disp(effect_wa)

% BV DBS ON
figure('Color', 'w','Position', [100, 100, 300, 500]);
daboxplot(M_BV_ICC_dbson_leadanalysis,'color', on_color,'whiskers' , 0, 'scatter', 1, 'scattersize',150, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories_leadanalysis,'fill',1,'boxspacing',0.1, 'boxwidth', 2);
ylabel('ICC', 'FontSize',25);
ylim([0.6,1]);
yticks(0.6:0.1:1);
xlim([0.5, 4.5]);
set(gca, 'XTickLabel', categories_leadanalysis, 'FontSize', 25, 'XTickLabelRotation', 0);

[p_au, h_au] = signrank(squeeze(M_BV_ICC_dbson_leadanalysis(:,2)), squeeze(M_BV_ICC_dbson_leadanalysis(:,3)));
[p_wu, h_wu] = signrank(squeeze(M_BV_ICC_dbson_leadanalysis(:,1)), squeeze(M_BV_ICC_dbson_leadanalysis(:,2)));
[p_wa, h_wa] = signrank(squeeze(M_BV_ICC_dbson_leadanalysis(:,1)), squeeze(M_BV_ICC_dbson_leadanalysis(:,3)));
disp([p_au, p_wu, p_wa])

effect_au = meanEffectSize(squeeze(M_BV_ICC_dbson_leadanalysis(:,2)), squeeze(M_BV_ICC_dbson_leadanalysis(:,3)), Effect="cliff",paired=true, ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wu = meanEffectSize(squeeze(M_BV_ICC_dbson_leadanalysis(:,1)), squeeze(M_BV_ICC_dbson_leadanalysis(:,2)), Effect="cliff",paired=true, ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wa = meanEffectSize(squeeze(M_BV_ICC_dbson_leadanalysis(:,1)), squeeze(M_BV_ICC_dbson_leadanalysis(:,3)), Effect="cliff",paired=true, ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(effect_au)
disp(effect_wu)
disp(effect_wa)

% BV DBS OFF
figure('Color', 'w','Position', [100, 100, 300, 500]);
daboxplot(M_BV_ICC_dbsoff_leadanalysis,'color', off_color,'whiskers' , 0, 'scatter', 1, 'scattersize', 150, 'scattercolors', {'k', 'k'}, 'scatteralpha', 0.5,'outsymbol','kx','xtlabels', categories_leadanalysis,'fill',1,'boxspacing',0.1, 'boxwidth', 2);
ylabel('ICC', 'FontSize',25);
ylim([0.6,1]);
yticks(0.6:0.1:1);
xlim([0.5, 4.5]);
set(gca, 'XTickLabel', categories_leadanalysis, 'FontSize', 25, 'XTickLabelRotation', 0);

[p_au, h_au] = signrank(squeeze(M_BV_ICC_dbsoff_leadanalysis(:,2)), squeeze(M_BV_ICC_dbsoff_leadanalysis(:,3)));
[p_wu, h_wu] = signrank(squeeze(M_BV_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_BV_ICC_dbsoff_leadanalysis(:,2)));
[p_wa, h_wa] = signrank(squeeze(M_BV_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_BV_ICC_dbsoff_leadanalysis(:,3)));
disp([p_au, p_wu, p_wa])

effect_au = meanEffectSize(squeeze(M_BV_ICC_dbsoff_leadanalysis(:,2)), squeeze(M_BV_ICC_dbsoff_leadanalysis(:,3)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wu = meanEffectSize(squeeze(M_BV_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_BV_ICC_dbsoff_leadanalysis(:,2)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
effect_wa = meanEffectSize(squeeze(M_BV_ICC_dbsoff_leadanalysis(:,1)), squeeze(M_BV_ICC_dbsoff_leadanalysis(:,3)), Effect="cliff",paired=true,ConfidenceIntervalType="bootstrap", ...
      BootstrapOptions=statset(UseParallel=true),NumBootstraps=10000);
disp(effect_au)
disp(effect_wu)
disp(effect_wa)