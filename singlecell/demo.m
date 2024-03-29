clear; % Clear the workspace

% Add necessary directories to the path
addpath('data', 'src', 'core');

folder_path = 'data'; % Set the folder path for the datasets

% Retrieve information about MAT files in the data directory
mat_files = dir(fullfile(folder_path, '*.mat'));

% Extract and save file names into a cell array
file_names = {mat_files.name};

alpha_NR = [100,10,1,0.1, 0.01, 0.001, 0.0001]; % Define alpha values for experimentation

% Initialize a results table to store dataset performances
results = table([], [], [], [], 'VariableNames', {'Dataset', 'NMI', 'ARI', 'BestAlpha'});

% Iterate through the datasets, adjust the range as necessary
for i = 1:min(30, length(file_names)) % Loop through first 30 datasets or total available
    load(fullfile(folder_path, file_names{i})); % Load the current dataset
    datasetID = file_names{i}(1:end-4); % Extract the dataset ID
    C = length(unique(labs)); % Identify the number of clusters
    rng(i, 'twister'); % Set random seed for reproducibility
    
    bestNMI = 0; bestARI = 0; bestalpha = 0; % Initialize best performance metrics
    
    % Iterate through each alpha value
    for alpha = alpha_NR
        [y, ~, ~, ~, ~] = SIMLR(in_X', C, 10, 0, 0, alpha); % Run SIMLR with the current alpha
        NMI_i = Cal_NMI(y, labs); % Calculate NMI
        ARI_i = Cal_ARI(y, labs, 'adjusted'); % Calculate ARI
        
        % Update best values if current results are better
        if NMI_i > bestNMI
            bestNMI = NMI_i; bestARI = ARI_i; bestalpha = alpha;
        end
    end
    
    % Display the best NMI value for the current dataset
    fprintf('The best NMI for dataset %s is %f with alpha %f\n', datasetID, bestNMI, bestalpha);
    
    % Append the new row of results
    results = [results; {datasetID, bestNMI, bestARI, bestalpha}];
end

% Save the results to a .mat file for future use
save('results_NR.mat', 'results');
