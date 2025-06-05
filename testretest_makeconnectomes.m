%% Make Connectomes - Functional Connectivity and Brain Variability Calculations

% This script calculates functional connectivity and brain variability values 

%% Step 1: Set Data Paths

dir = '/YOUR/PATH/HERE/'; % Replace with the path to where your CONN files are stored.
dir_postv1 = [dir '/YOUR/PATH/HERE/'] % Replace with the path to the session 01 CONN file 
dir_postv2 = [dir '/YOUR/PATH/HERE/'] % Replace with the path to the session 02 CONN file 

output = '/YOUR/PATH/HERE/';

%% Step 2: Initalize Values

% Input variable corresponding to the preprocessed fMRI dataset.
results_dir= dir_postv1; % change here to choose session 01 or 02
cd (results_dir)
 
% Input the number of subjects.
subj = 1:16; 

prefix='qsm_atlas'; % MAKE SURE THIS MATCHES THE PREFIX OF YOUR ATLAS's ROIs

% Input the conditions.
conditions = load('_list_conditions.mat');
conditions.allnames %conditions 5-8 are DBS OFF RUN 1, DBS OFF RUN 2, DBS ON RUN 1, DBS ON RUN 2
condition = 5:8;

% Input the brain network.
brain_network = 1; % 1=whole-brain, 2=motor, 3=limbic, 4=associative

% Load one example dataset to extract ROI names.
ROI = load('ROI_Subject001_Condition000.mat'); 

% Create a list of all possible ROIs.
names = string(ROI.names);
disp(names); % Print all ROI names to verify content.

%% Step 3: Choose network ROIs by Name

% Initialize the ROI names.

roi_labels = {[prefix '.OFC'], [prefix '.Ins'], [prefix '.Cingulate-Ant'], ...
    [prefix '.Temporal-Pole'], [prefix '.Caud'], [prefix '.Put'], ...
    [prefix '.Sub'], [prefix '.TL'], [prefix '.Amy'], [prefix '.Hip'], ...
    [prefix '.ParaH'], [prefix '.Front'], ...
    [prefix '.Rol'], [prefix '.Ang'], [prefix '.Supr'], ...
    [prefix '.Precent'], [prefix '.Postcent'], [prefix '.Paracent'], ...
    [prefix '.Supp'], [prefix '.Cerebel'], [prefix '.Verm'], 'GPe', 'GPi', 'STN',}; % These labels are based on our combined QSM and QSM DGM atlas

% Initialize the ROI index.
roi_idx = cell(size(roi_labels));

% Fill the ROI index with indices from the preprocessed data that match the ROI names.
for i = 1:length(roi_labels)
    roi_idx{i} = find(contains(names, roi_labels{i}));
end

%Assign each ROI index to a corresponding variable.
[ofc, insula, cing_ant, temp_pole, caud, put, sub_nig, thal, amyg, hip, parahip, ...
    front, rol_op, angular, supramarg, precent, postcent, ...
    paracent, supp, cereb, vermis, gpe, gpi, stn] = roi_idx{:};

% Choose motor network ROIs.
allrois_idx_motor=[precent postcent paracent supp caud put sub_nig ... 
    thal gpi gpe stn cereb vermis]; %tot=60

disp(allrois_idx_motor) %test names - should return a 1x60 double.

% Choose limbic network ROIs.
allrois_idx_limb=[ofc insula cing_ant temp_pole caud put sub_nig ...
    thal amyg hip parahip gpi gpe stn]; %tot=48

disp(allrois_idx_limb) %test names - should return a 1x48 double.

% Choose associative network ROIs.
allrois_idx_assoc=[front rol_op angular supramarg caud put sub_nig ...
    thal gpi gpe stn]; %tot=46

disp(allrois_idx_assoc) %test names - should return a 1x46 double.

% Choose wholebrain ROIs.
whole=find(contains(names,prefix) | contains(names,'GPi')| contains(names,'GPe')| contains(names,'STN')); % This should call 212 ROIs: 206 are the entire QSM atlas, and 6 are from the QSM DGM atlas (GPi L/R, GPe L/R, and STN L/R)
temp_idx1=find(contains(names(whole),[prefix '.Vermis-10'])); % mark beginning of WM tracts

disp(whole(temp_idx1)) % test for vermis-10

temp_idx2=find(contains(names(whole),[prefix '.Caudate-Nucleus'])); % mark end of WM tracts

disp(whole(temp_idx2)) % test for caudate-nuc

whole_excludeWM=whole([1:temp_idx1 temp_idx2:length(whole)]); % patch markers together

disp(whole_excludeWM) % test for exclusion of WM tracts %144 ROIs

temp_idx3=find(contains(names(whole_excludeWM),[prefix '.Globus-Pallidus'])); % mark GP

disp(whole_excludeWM(temp_idx3)) % test for GP

whole_excludeWM_excludeGP = whole_excludeWM([1:temp_idx3-1 temp_idx3+2:length(whole_excludeWM)]); %Remove GP from QSM, since we're already pulling it from QSM DGM

allrois_idx_whole=whole_excludeWM_excludeGP;
disp(allrois_idx_whole); %test names - should return a 1x142 double.

% Choose network based on user input.
if brain_network==1 
    allrois_idx = allrois_idx_whole;
elseif brain_network==2
    allrois_idx = allrois_idx_motor;
elseif brain_network==3
    allrois_idx = allrois_idx_limb;
elseif brain_network==4
    allrois_idx = allrois_idx_assoc;
end

% Check that ROI names are correct for specified brain network.
allrois_names=names(allrois_idx);
disp(allrois_names); 

% Clean up names by removing atlas prefixes and check that they are correct. 
corrname=string(allrois_names);
corrname = erase(corrname, {'qsm_atlas.','qsm_dgm.'}); 
match="(" + wildcardPattern + ")";
corrname=erase(corrname,match);
corrname = strip(corrname);
disp(corrname);

%% Step 4: Compute functional connectivity & brain variability matrices

% Build an empty matrix and vector.
fcs = zeros(length(subj),length(condition), length(allrois_idx),length(allrois_idx));
brain_var = zeros(length(subj), length(condition),length(allrois_idx));

% Load subjects ROI files.
cd(results_dir)

for i = 1:length(subj)
    for l=1:length(condition)
        if subj(i)<10
            j=num2str(subj(i),'%02d');
            ROI = load(['ROI_Subject0' j '_Condition000.mat']);
        else
            ROI = load(['ROI_Subject0' num2str(subj(i)) '_Condition000.mat']);
        end

        w = max(0,ROI.conditionsweights{condition(l)}{1}); %extract condition timeseries (see https://www.nitrc.org/forum/message.php?msg_id=17373)
        idx = find( w>0 );
        ROI.data = cellfun(@(x)x(idx,:), ROI.data, 'uni',0); %now this should be 168 TP

        % Generate functional matrices.
        for n=1:length(allrois_idx)
            for m=1:length(allrois_idx)
                fcs(i,l,m,n)=corr(ROI.data{allrois_idx(n)}(:),ROI.data{allrois_idx(m)}(:));
            end
        end

        % Generate brain variability values.
        for n=1:length(allrois_idx)
            brain_var(i,l,n)=std(ROI.data{allrois_idx(n)}(:));
        end
    end
end

% Check that fcs & brainvar contains all the matrices.
size(fcs) 
size(brain_var)

%% Step 5: Save connectomes & brain var

% Make an output directory.
cd(output)

if ~exist('testretest_connectomes', 'dir')
mkdir('testretest_connectomes')
end

cd ('/YOUR/PATH/HERE/testretest_connectomes/')

% Determine the filename based on the directory and brain network.
if strcmp(results_dir, dir_postv1)
    if brain_network == 1
        filename = 'results_whole_v1.mat';
    elseif brain_network == 2
        filename = 'results_motor_v1.mat';
    elseif brain_network == 3
        filename = 'results_limbic_v1.mat';
    elseif brain_network == 4
        filename = 'results_assoc_v1.mat';
    end
elseif strcmp(results_dir, dir_postv2)
    if brain_network == 1
        filename = 'results_whole_v2.mat';
    elseif brain_network == 2
        filename = 'results_motor_v2.mat';
    elseif brain_network == 3
        filename = 'results_limbic_v2.mat';
    elseif brain_network == 4
        filename = 'results_assoc_v2.mat';
    end
end

% Save the file.
save(filename, 'fcs', 'brain_var', 'names', 'allrois_idx', 'subj', 'corrname');

