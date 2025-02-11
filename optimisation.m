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

%%% code used for adding lac_L[e] exchange reaction to all the models
% Define the folder path
folder_path = 'D:/MATLAB/models_wo_outliers';

% Get a list of all .mat files in the folder
mat_files = dir(fullfile(folder_path, '*.mat'));

% Loop through each .mat file
for i = 1:numel(mat_files)
    % Load the matfile
    mat_data = load(fullfile(folder_path, mat_files(i).name));
    
    % Check if 'model' is a field in the loaded matfile
    if isfield(mat_data, 'model')
        model = mat_data.model;
        
        % Add exchange reaction for 'lac_L[e]' to the model
        model = addExchangeRxn(model, {'lac_L[e]'}, 0, 1000);
        
        % Save the modified model back to the same file
        save(fullfile(folder_path, mat_files(i).name), 'model');
    end
end
%%%%%%%%%%%
% Find the index of the 'EX_glc(e)' reaction in the model
rxnIndex = find(ismember(model.rxns, 'EX_glc(e)'));

% Set the new upper bound value
newUpperBoundValue = -0.1;

% Update the upper bound of the reaction with the new value
model.ub(rxnIndex) = newUpperBoundValue;
%%%%% % Find the index of the reaction ADK1 in model.rxns
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
