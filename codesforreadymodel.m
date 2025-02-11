%converting mets () to []
% Define the folder path
folderPath = "D:\MATLAB\models_wo_outliers";

% List all MAT files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Initialize counter for models containing 'co2[e]'
models_with_co2_e = 0;

% Loop through each MAT file
for i = 1:length(mat_files)
    % Load the model from the MAT file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Check if 'co2[e]' is present in the metabolites
    if any(contains(model.mets, 'co2[e]'))
        models_with_co2_e = models_with_co2_e + 1;
        
        % Iterate through each metabolite
        for j = 1:length(model.mets)
            % Replace 'co2(e)' with 'co2[e]'
            model.mets{j} = strrep(model.mets{j}, 'co2(e)', 'co2[e]');
        end
        
        % Save the modified model
        save(fullfile(folderPath, ['modified_', mat_files(i).name]), 'model');
    end
end

disp(['Number of models with co2[e]: ', num2str(models_with_co2_e)]);
%% finding if one rtn is present in the whole model 
% Define the folder path
folderPath = 'D:/MATLAB/models_wo_outliers/modified';

% List all MAT files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Initialize a counter for models containing 'EX_lac_L' in model.rxns
lac_L_present_count = 0;

% Loop through each MAT file
for i = 1:length(mat_files)
    % Load the model from the MAT file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Check if 'EX_lac_L' is present in the reaction IDs
    if any(strcmp(model.mets, 'oh1[c]'))
        lac_L_present_count = lac_L_present_count + 1;
        disp(['oh1[c] found in model ', mat_files(i).name]);
    end
end

disp(['Total number of models containing oh1[c]: ', num2str(lac_L_present_count)]);

%% finding duplicate metabolite, if there are same metabolite more than once
% Set the path to the model directory
modelDir = 'D:/MATLAB/models_wo_outliers/modified';

% List all MAT files in the directory
mat_files = dir(fullfile(folderPath, '*.mat'));

% Loop through each MAT file
for i = 1:length(mat_files)
    % Load the model from the MAT file
    filename = fullfile(folderPath, mat_files(i).name);
    load(filename, 'model');
    
    % Find repeated metabolites
    duplicate_mets = unique(model.mets([diff(sort(model.mets)) == 0; false]));
    
    % Display model name and repeated metabolites, if any
    if ~isempty(duplicate_mets)
        disp(['Model: ', mat_files(i).name, ' - Repeated Metabolites: ', strjoin(duplicate_mets, ', ')]);
    end
end
%% convert rxns fron [e] to (e)

folder_path = 'D:/MATLAB/models_wo_outliers/modified';

% Get a list of all files in the folder
model_files = dir(fullfile(folder_path, '*.mat'));

% Loop through each model file
for file_idx = 1:numel(model_files)
    % Load the model
    file_path = fullfile(folder_path, model_files(file_idx).name);
    load(file_path); % Load the model from the file
    
    % Replace '[' with '(' and ']' with ')' in model.rxns
    for i = 1:numel(model.rxns)
        model.rxns{i} = strrep(model.rxns{i}, '[', '(');
        model.rxns{i} = strrep(model.rxns{i}, ']', ')');
    end
    
    % Save the modified model back to the same file
    save(file_path, 'model');
end
%% code for FBA then biomass to 50% then FVA then 
% Folder path
folder_path = 'D:\MATLAB\models_wo_outliers\modified';

% Get a list of all .mat files in the specified folder starting with "TT"
mat_files = dir(fullfile(folder_path, '*CN*.mat'));

% Loop through each .mat file
for i = 1:numel(mat_files)
    % Load the model
    filename = fullfile(folder_path, mat_files(i).name);
    disp(['Processing model: ' filename]);
    load(filename);
    
    % Code 1: Set biomass flux to 50%
    solution = optimizeCbModel(model);
    biomass_index = find(contains(model.rxns, 'biomass_reaction'));
    biomass_flux = solution.f * 0.5;
    model.lb(biomass_index) = biomass_flux;
    model.ub(biomass_index) = biomass_flux;

    % Code 2: Perform Flux Variability Analysis (FVA)
    [minFlux, maxFlux] = fluxVariability(model);

    % Compute flux span
    if numel(maxFlux) == numel(minFlux)
        flux_span = maxFlux - minFlux;

        % Prepare data for Excel
        data_flux_span = [model.rxns, num2cell(flux_span)];

        % Prepare Excel file
        excel_filename = fullfile(folder_path, [mat_files(i).name(1:end-4) '_flux_span.xlsx']);

        % Write results to Excel
        xlswrite(excel_filename, {'Reaction', 'Flux Span'}, 'Sheet1', 'A1');
        xlswrite(excel_filename, data_flux_span, 'Sheet1', 'A2');

        disp(['Flux span computed for ' num2str(length(flux_span)) ' reactions.']);
        disp(['Results saved to ' excel_filename]);
    else
        disp(['Error: The dimensions of maxFlux and minFlux in ' filename ' are not compatible.']);
    end
end
%% % Define the path to the Excel file
excelFile = 'C:/Users/hp/Downloads/biomass 50% FBA & FVA .xlsx';

% Read all sheet names from the Excel file
[~, sheetNames] = xlsfinfo(excelFile);

% Initialize a structure to store matrices
matrices = struct();

% Loop through each sheet
for sheetIdx = 1:numel(sheetNames)
    % Read data from the current sheet
    [~,~,rawData] = xlsread(excelFile, sheetNames{sheetIdx});
    
    % Get the reaction names from column A
    reactions = rawData(2:end, 1);
    
    % Get the flux values from column B onwards
    fluxValues = rawData(2:end, 2:end);
    
    % Convert numeric values to double, and leave NaNs as they are
    fluxValues_numeric = cellfun(@(x) isnumeric(x) && ~isnan(x), fluxValues);
    fluxValues_numeric(~fluxValues_numeric) = 0; % Convert logicals to 0 for non-numeric values
    fluxValues_numeric = cellfun(@(x) x, fluxValues, 'UniformOutput', false); % Convert NaN back to NaN
    
    % Store the matrix in the structure
    matrices.(genvarname(sheetNames{sheetIdx})) = fluxValues_numeric;
end

% Access matrices using field names
% For example, to access the matrix from the third sheet:
% thirdMatrix = matrices.SheetName3;