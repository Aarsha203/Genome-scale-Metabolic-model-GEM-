%%% after changing glu-D upper bound to -0.1 and the conducting F value
% Define the directory where the MATLAB files are located
folderPath = 'D:/MATLAB/models_wo_outliers';

% List all MATLAB files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Loop through each MATLAB file
for i = 1:length(mat_files)
    % Load the model from the MATLAB file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Find the index of 'EX_glc(e)' reaction
    rxnIndex = find(ismember(model.rxns, 'EX_glc(e)'));
    if ~isempty(rxnIndex)
        % Set the upper bound of 'EX_glc(e)' reaction to -0.1
        newUpperBoundValue = -0.1;
        model.ub(rxnIndex) = newUpperBoundValue;
        
        % Check if 'EX_lac_L[e]' is present in the model
        lac_index = find(strcmpi('EX_lac_L[e]', model.rxns));
        if ~isempty(lac_index)
            % Change objective to 'EX_lac_L[e]'
            model = changeObjective(model, 'EX_lac_L[e]');
            
            % Optimize the model to find the optimal objective value
            sol = optimizeCbModel(model);
            
            % Print the optimal objective value (f)
            disp(['Optimal objective value (f) for ', mat_files(i).name, ': ', num2str(sol.f)]);
        else
            disp(['EX_lac_L[e] reaction not found in the model for ', mat_files(i).name, '.']);
        end
    else
        disp(['EX_glc(e) reaction not found in the model for ', mat_files(i).name, '.']);
    end
end

%% 
% Find the index of the biomass reaction(for individual model)
biomassIdx = find(strcmp(model.rxns, 'biomass_reaction'));

if ~isempty(biomassIdx)
 % Set the lower bound of the biomass reaction to 0.01
    model.lb(biomassIdx) = 0.01;
else
    fprintf('Biomass reaction not found in the model.\n');
end
% Change the objective to 'EX_lac_L[e]' (individual)
model = changeObjective(model, 'EX_lac_L[e]');
%change of biomass reaction lb into 0.01Then, it checks if 'EX_lac_L[e]' is
%present in the model, and if not, it adds it as an exchange reaction. then
%it optimises the model (for the folder)
% Define the folder path
folderPath = 'D:/MATLAB/models_wo_outliers';

% Get a list of all .mat files in the folder
modelFiles = dir(fullfile(folderPath, '*.mat'));

% Initialize a cell array to store the results
results = cell(length(modelFiles), 2);

% Loop through each .mat file
for i = 1:length(modelFiles)
    % Load the model from the .mat file
    modelFilePath = fullfile(folderPath, modelFiles(i).name);
    load(modelFilePath, 'model');
    
    % Find the index of the biomass reaction
    biomassIdx = find(strcmp(model.rxns, 'biomass_reaction'));

    % Check if the biomass reaction exists
    if ~isempty(biomassIdx)
        % Set the lower bound of the biomass reaction to 0.01
        model.lb(biomassIdx) = 0.01;

        % Check if 'EX_lac_L[e]' is present in the model, if not, add it
        if ~ismember('EX_lac_L[e]', model.rxns)
            model = addExchangeRxn(model, {'lac_L[e]'}, 0, 1000);
        end

        % Check if 'EX_lac_L[e]' is now present in the model
        if ismember('EX_lac_L[e]', model.rxns)
            % Change the objective to 'EX_lac_L[e]'
            model = changeObjective(model, 'EX_lac_L[e]');

            % Optimize the model
            sol = optimizeCbModel(model);
            
            % Store the result (f value) in the results cell array
            results{i, 1} = modelFiles(i).name;
            results{i, 2} = sol.f;
        else
            fprintf('Objective reaction not found in model: %s\n', modelFiles(i).name);
        end
    else
        % Biomass reaction not found
        fprintf('Biomass reaction not found in model: %s\n', modelFiles(i).name);
    end
end

% Display the results
disp('Results:');
disp(results);

% Write the results to an Excel file
resultFilePath = fullfile(folderPath, 'result.xlsx');
xlswrite(resultFilePath, results);
%% this recognise if the model have already EX_lac[e]if not it will add
%%% code used for adding lac_L[e] exchange reaction to all the models
% Define the folder path
folderPath = 'D:/MATLAB/models_wo_outliers';

% List all MATLAB files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Loop through each MATLAB file
for i = 1:length(mat_files)
    % Load the model from the MATLAB file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Check if the model already contains the exchange reaction for 'pyr(e)'
    if ~any(contains(model.rxns, 'EX_pyr(e)'))
        % Add exchange reaction for 'pyr(e)' to the model
        model = addExchangeRxn(model, {'pyr(e)'}, 0, 1000);
        
        % Save the modified model
        save(fullfile(folderPath, ['modified_', mat_files(i).name]), 'model');
    else
        disp(['EX_pyr(e) already present in ', mat_files(i).name]);
    end
end

%% change metabolite co2(e) to co2[e]
% Define the directory where the MATLAB files are located
folderPath = 'D:/MATLAB/models_wo_outliers';

% List all MATLAB files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Initialize arrays to store model names
models_with_co2_e = {};
models_without_co2_e = {};

% Loop through each MATLAB file
for i = 1:length(mat_files)
    % Load the model from the MATLAB file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Check if 'co2(e)' is present in the metabolites
    if any(contains(model.mets, 'co2(e)'))
        % Iterate through each metabolite
        for j = 1:length(model.mets)
            % Replace 'co2(e)' with 'co2[e]'
            model.mets{j} = strrep(model.mets{j}, 'co2(e)', 'co2[e]');
        end
        
        % Save the modified model
        save(fullfile(folderPath, ['modified_', mat_files(i).name]), 'model');
        
        % Store the model name in the array
        models_with_co2_e = [models_with_co2_e, mat_files(i).name];
    else
        % Store the model name in the array
        models_without_co2_e = [models_without_co2_e, mat_files(i).name];
    end
end

% Display the models with co2(e) and those without co2(e)
disp('Models with co2(e):');
disp(models_with_co2_e);
disp('Models without co2(e):');
disp(models_without_co2_e);

%%%%%%%%%%%
%% 
% Find the index of the 'EX_glc(e)' reaction in the model
rxnIndex = find(ismember(model.rxns, 'EX_glc(e)'));

% Set the new upper bound value
newUpperBoundValue = -0.1;

% Update the upper bound of the reaction with the new value
model.ub(rxnIndex) = newUpperBoundValue;
%% 
%%% % Find the index of the reaction ADK1 in model.rxns
reactionIndex = find(strcmp(model.rxns, 'CO2tm'));
% Check if the reaction ADK1 exists in the model
if ~isempty(reactionIndex)
% Get the gene associated with the reaction ADK1 from model.gene
associatedGene = model.genes(reactionIndex);
% Print the associated gene
fprintf('The gene associated with the reaction ADK1 is: %s\n', associatedGene{:});
else
fprintf('The reaction ADK1 does not exist in the model.\n');
end
%% rxns to gpr 
% Load the Excel file
filename = 'C:\Users\hp\OneDrive\Desktop\IIT_M\collateral lethal and subsytems.xlsx';
sheet = 1;
range = 'A:A';
rxns = readcell(filename, 'Sheet', sheet, 'Range', range);

% Initialize cell array to store grRules
grRules = cell(size(rxns));
%% find the rxns in BL 21.03
folder_path = "D:\MATLAB\models_wo_outliers";

% Get a list of all MAT files in the folder
mat_files = dir(fullfile(folder_path, 'BN*.mat'));

% Initialize a structure to store the counts for each model
model_counts = struct();

% Loop through each MAT file
for i = 1:numel(mat_files)
    % Load the MAT file
    mat_data = load(fullfile(folder_path, mat_files(i).name));
    
    % Extract the model name from the MAT file name (assuming the file name format is 'modelName.mat')
    [~, model_name, ~] = fileparts(mat_files(i).name);
    
    % Check if the 'model' variable exists in the MAT file
    if isfield(mat_data, 'model') && isfield(mat_data.model, 'rxns')
        % Get the count of entries in 'model.rxns'
        rxns_count = numel(mat_data.model.rxns);
        
        % Store the count for the current model
        model_counts.(model_name) = rxns_count;
    else
        % If the 'model' variable doesn't exist or 'rxns' field is missing, store 0
        model_counts.(model_name) = 0;
    end
end

% Display the counts for each model
disp('Number of model.rxns for each individual model:');
disp(model_counts);
%% % Define the folder containing the MAT files
folder_path = "D:\MATLAB\models_wo_outliers";

% Define the models of interest
models_of_interest = {'BN16', 'BN26'};

% Initialize a structure to store the model.rxns for each model
model_rxns = struct();

% Loop through each model of interest
for i = 1:numel(models_of_interest)
    model_name = models_of_interest{i};
    
    % Load the MAT file for the current model
    mat_file_path = fullfile(folder_path, [model_name '.mat']);
    mat_data = load(mat_file_path);
    
    % Check if the 'model' variable exists in the MAT file
    if isfield(mat_data, 'model') && isfield(mat_data.model, 'rxns')
        % Store the model.rxns in the structure
        model_rxns.(model_name) = mat_data.model.rxns;
    else
        % Display a message if 'model' or 'model.rxns' is missing
        disp(['Model ' model_name ' or model.rxns is missing.']);
    end
end

% Find missing reactions for each model
for i = 1:numel(models_of_interest)
    model_name = models_of_interest{i};
    
    % Find missing reactions by comparing with the first model
    if i > 1
        missing_rxns = setdiff(model_rxns.(models_of_interest{1}), model_rxns.(model_name));
    else
        missing_rxns = [];
    end
    
    % Print out the missing reactions for the current model
    disp(['Missing reactions for ' model_name ':']);
    disp(missing_rxns);
end
%% code to find which reactions abbreviations are present in the model.
% Define the path to the directory containing the models
models_path = 'D:\MATLAB\models_wo_outliers';

% Get a list of all MAT files starting with 'BN' in the directory
model_files = dir(fullfile(models_path, 'BN*.mat'));

% Initialize variables to store model names
models_with_co2 = {};
models_without_co2 = {};

% Loop through each MAT file
for i = 1:numel(model_files)
    % Load the MAT file
    model_data = load(fullfile(models_path, model_files(i).name));
    
    % Check if the model contains the reaction EX_co2(e)
    if isfield(model_data.model, 'rxns') && any(strcmp(model_data.model.rxns, 'EX_pyr(e)'))
        % Store the model name if the reaction is present
        models_with_co2 = [models_with_co2, model_files(i).name];
    else
        % Store the model name if the reaction is not present
        models_without_co2 = [models_without_co2, model_files(i).name];
    end
end

% Display the models with and without EX_co2(e)
disp('Models with EX_co2(e):');
disp(models_with_co2);

disp('Models without EX_co2(e):');
disp(models_without_co2);

%%%
%% adding an exchnage reaction to the whole matfile
% Define the paths for input and output folders
input_folder_path = 'D:\MATLAB\models_wo_outliers';
output_folder_path = 'D:\MATLAB\models_with_co2';

% Get a list of all .mat files starting with 'BN' in the input folder
mat_files = dir(fullfile(input_folder_path, 'BN*.mat'));

% Create the output folder if it doesn't exist
if ~exist(output_folder_path, 'dir')
    mkdir(output_folder_path);
end

% Loop through each .mat file
for i = 1:numel(mat_files)
    % Load the matfile
    mat_data = load(fullfile(input_folder_path, mat_files(i).name));
    
    % Check if 'model' is a field in the loaded matfile
    if isfield(mat_data, 'model')
        model = mat_data.model;
        
        % Check if the model already has 'pyr[e]' exchange reaction
        if ~any(contains(model.rxns, 'EX_co2(e)'))
            % Add exchange reaction for 'pyr[e]' to the model
            model = addExchangeRxn(model, {'co2(e)'}, 0, 1000);
            
            % Save the modified model to the output folder
            output_file_path = fullfile(output_folder_path, mat_files(i).name);
            save(output_file_path, 'model');
        end
    end
end
%% changimg pyruvate from EX_lac_L[e] to EX_lac_L(e)

%%changimg pyruvate 
% Define the path to the directory containing the MAT files
folder_path = 'D:/MATLAB/models_wo_outliers/modified';

% Get a list of all MAT files in the folder
mat_files = dir(fullfile(folder_path, '*.mat'));

% Loop through each MAT file
for i = 1:numel(mat_files)
    % Load the MAT file
    mat_data = load(fullfile(folder_path, mat_files(i).name));
    
    % Check if 'model' is a field in the loaded matfile
    if isfield(mat_data, 'model')
        % Replace 'EX_pyr[e]' with 'EX_pyr(e)' in model.rxns
        if isfield(mat_data.model, 'rxns')
            mat_data.model.rxns = regexprep(mat_data.model.rxns, 'EX_pyr\[e\]', 'EX_pyr(e)');
            
            % Overwrite the existing MAT file with the modified model
            save(fullfile(folder_path, mat_files(i).name), '-struct', 'mat_data');
        else
            warning(['No rxns field found in ' mat_files(i).name]);
        end
    else
        warning(['No model field found in ' mat_files(i).name]);
    end
end
%% joined optimisation code 
folderPath = 'D:/MATLAB/models_wo_outliers/modified';
 modelFiles = dir(fullfile(folderPath, '*.mat'));
 for i = 1:length(modelFiles)

modelFilePath = fullfile(folderPath, modelFiles(i).name);
    load(modelFilePath, 'model');

sol = optimizeCbModel(model);

fprintf('Model: %s\n', modelFiles(i).name);
    fprintf('Optimal objective value: %f\n', sol.f);
 end

% Write the results to an Excel file
filename = fullfile(folderPath, 'model_results.xlsx');
xlswrite(filename, results);
%% check which reactions is present in all the models
% Define the path to the directory containing the MAT files
folder_path = 'D:\MATLAB\models_wo_outliers';

% Get a list of all MAT files in the folder
mat_files = dir(fullfile(folder_path, 'BN*.mat'));

% Initialize a cell array to store model.rxns containing 'EX_lac_L[e]'
rxns_with_ex_co2 = {};

% Loop through each MAT file
for i = 1:numel(mat_files)
    % Load the MAT file
    mat_data = load(fullfile(folder_path, mat_files(i).name));
    
    % Check if 'model' is a field in the loaded matfile
    if isfield(mat_data, 'model')
        % Check if 'model.rxns' is a field in the 'model' structure
        if isfield(mat_data.model, 'rxns')
            % Check if 'EX_lac_L[e]' is present in 'model.rxns'
            if any(contains(mat_data.model.rxns, 'EX_co2(e)'))
                % Add the file name to the list
                rxns_with_ex_co2 = [rxns_with_ex_co2, mat_files(i).name];
            end
        else
            warning(['No rxns field found in ' mat_files(i).name]);
        end
    else
        warning(['No model field found in ' mat_files(i).name]);
    end
end

% Display the model files containing 'EX_lac_L[e]' in model.rxns
disp('MAT files containing EX_co2(e) in model.rxns:');
disp(rxns_with_ex_25hvitd3);
%% code to analyse single and double deletion
% Define the directory where the MATLAB files are located
folderPath = 'D:/MATLAB/kannan/models'; %change accordingly 

% List all MATLAB files in the directory will run for all the file in the
% folder which has mat extension
mat_files = dir(fullfile(folderPath, '*.mat'));

% Initialize results cell array, define how you want to give your result
% output, it will be shown like that you can also opt for table format 
results = cell(length(mat_files), 5);

% Loop through each MATLAB file
for i = 1:length(mat_files)
    % Load the model from the MATLAB file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');

    % Extract model.rxns, already loaded
    model_rxns = model.rxns;

    % Find the intersection of 'media_rxn' and 'model_rxns'
    int_sect = intersect(media_rxn, model_rxns);

    % Find the indices of 'media_rxn' in 'int_sect'
    index_num = find(ismember(media_rxn, int_sect));

    % Set the lower bounds based on 'media_flux' which is already loaed in
    % work space 
    for j = index_num
        a = media_flux(j);
        low_rxn = media_rxn(j);
        low_bound = -(a);
        model = changeRxnBounds(model, low_rxn, low_bound, 'l');
    end

    % Perform single gene deletion analysis
    [grRatio, ~, ~, ~, ~] = singleGeneDeletion(model);

    % Perform double gene deletion analysis
    [grRatioDble, ~, ~] = doubleGeneDeletion(model);

    % Calculate non-scaled epistasis
    [m, ~] = size(grRatio);
    b = zeros(m);
    for j = 1:m
        for k = 1:m
            b(j, k) = grRatio(j) * grRatio(k);
        end
    end
    dia_grRatio = triu(b, 1);
    dia_grRatiodble = triu(grRatioDble, 1);
    d = dia_grRatiodble - dia_grRatio;
    positive_epistasis = sum(d(:) > 1e-6);
    negative_epistasis = sum(d(:) < -1e-6);

    % Calculate scaled epistasis
    e = abs(d);
    f = d ./ e;
    F_positive_epistasis = sum(f(:) > 1e-6);
    F_negative_epistasis = sum(f(:) < -1e-6);

    % Store the results
    results{i, 1} = mat_files(i).name;
    results{i, 2} = positive_epistasis;
    results{i, 3} = negative_epistasis;
    results{i, 4} = F_positive_epistasis;
    results{i, 5} = F_negative_epistasis;
end

% Display the results
disp('Results:');
disp(results);

% Write the results to an Excel file
resultFilePath = fullfile(folderPath, 'result.xlsx');
xlswrite(resultFilePath, results);
%% 
% Define the directory where the MATLAB files are located
folderPath = 'D:\MATLAB\models_wo_outliers';

% List all MATLAB files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Loop through each MATLAB file
for i = 1:length(mat_files)
    % Load the model from the MATLAB file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Iterate through each metabolite
    for j = 1:length(model.mets)
        % Check if the metabolite is 'co2(e)'
        if strcmp(model.mets{j}, 'co2(e)')
            % 'co2(e)' is already present, no need to make any changes
            disp(['co2[e] is already present in ', mat_files(i).name]);
            break; % Exit the loop
        elseif strcmp(model.mets{j}, 'co2[e]')
            % 'co2[e]' is already present, no need to make any changes
            disp(['co2[e] is already present in ', mat_files(i).name]);
            break; % Exit the loop
        elseif contains(model.mets{j}, 'co2(e)')
            % Replace 'co2(e)' with 'co2[e]'
            model.mets{j} = strrep(model.mets{j}, 'co2(e)', 'co2[e]');
            disp(['Replaced co2(e) with co2[e] in ', mat_files(i).name]);
        end
    end
    
    % Save the modified model
    save(fullfile(folderPath, ['modified_', mat_files(i).name]), 'model');
end
%% FVA and calculating jaccard index
% Define the directory where the MATLAB files are located
folderPath = 'D:/MATLAB/models_wo_outliers';

% List all MATLAB files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Loop through each MATLAB file
for i = 1:length(mat_files)
    % Load the model from the MATLAB file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Calculate flux variability
    [minFlux, maxFlux] = fluxVariability(model);

    % Calculate Jaccard index
    J = fvaJaccardIndex(minFlux, maxFlux);

    % Display the Jaccard index along with the filename
    disp(['File: ', mat_files(i).name, ', Jaccard index (J) = ', num2str(J)]);
end
%% adding a exchange rtn to all the mat files and changing the og matfile 
% Define the folder path
folderPath = 'D:/MATLAB/models_wo_outliers/modified';
% List all MAT files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Loop through each MAT file
for i = 1:length(mat_files)
    % Load the model from the MAT file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Check if 'pyr[e]' exchange reaction is already present
    if ~any(strcmp(model.rxns, 'EX_tcynt(e)'))
        % Add 'pyr[e]' exchange reaction to the model
        model = addExchangeRxn(model, {'tcynt[e]'}, 0, 1000);
        
        % Save the modified model back to the same MAT file
        save(filename, 'model');
        
        disp(['Added tcynt[e] exchange reaction to ', mat_files(i).name]);
    else
        disp(['EX_tcynt(e) already present in ', mat_files(i).name]);
    end
end
%% 
%%optimisiation for all models 
% Folder containing the models
folderPath = 'D:/MATLAB/models_wo_outliers/modified';

% Get list of model files
modelFiles = dir(fullfile(folderPath, '*.mat'));

% Initialize results cell array
results = cell(length(modelFiles), 2);

% Loop through each model file
for i = 1:length(modelFiles)
    % Load the model
    modelFilePath = fullfile(folderPath, modelFiles(i).name);
    load(modelFilePath, 'model');

    % Optimize the model
    sol = optimizeCbModel(model);

    % Store results
    results{i, 1} = modelFiles(i).name;
    results{i, 2} = sol.f;

    % Print the results
    fprintf('Model: %s\n', modelFiles(i).name);
    fprintf('Optimal objective value: %f\n', sol.f);
end

% Write the results to an Excel file
filename = fullfile(folderPath, 'model_results.xlsx');
xlswrite(filename, results);





