clear; % Clear the workspace
rng('default');
% Add necessary directories to the path
addpath('data', 'src');

folder_path = 'data'; % Set the folder path for the datasets

% Retrieve information about MAT files in the data directory
mat_files = dir(fullfile(folder_path, '*.mat'));

% Extract and save file names into a cell array
file_names = {mat_files.name};

lambda_TopoLa = 0.05; 

% Initialize a results table to store dataset performances
results = table([], [], [], 'VariableNames', {'Dataset', 'NMI', 'ARI'});

% Iterate through the datasets
for i = 1:min(30, length(file_names)) % Loop through the first 30 datasets or total available
    load(fullfile(folder_path, file_names{i})); % Load the current dataset
    datasetID = file_names{i}(1:end-4); % Extract the dataset ID
    C = length(unique(labs)); % Identify the number of clusters

    [y, S_raw, S_TopoLa, ydata, alphaK, timeOurs, converge, LF] = SIMLR_TopoLa(in_X', C, 10, 0, 0, lambda_TopoLa); % Run SIMLR with the current alpha
    NMI_i = Cal_NMI(y, labs); % Calculate NMI
    ARI_i = Cal_ARI(y, labs, 'adjusted'); % Calculate ARI
    
    

    % Append the new row of results
    results = [results; {datasetID, NMI_i, ARI_i}];
end

% Save the final results to a CSV file
writetable(results, 'Result.csv');
