addpath('src');

% Get all MAT files in the Data folder
matFiles = dir('data/*.mat');
lambda = 1e-4; % Fix alpha value

% Initialize result storage
NMI_results = zeros(length(matFiles), 1);
ARI_results = zeros(length(matFiles), 1);

% Loop through each MAT file
for i = 1:length(matFiles)
    % Get the full file path
    filePath = fullfile(matFiles(i).folder, matFiles(i).name);
    
    % Load the MAT file and display variables
    data = load(filePath);
    disp(['Processing file: ' matFiles(i).name]);
    disp('Loaded variables:');
    disp(fieldnames(data));
    
    % Check if variables A and lab exist
    if isfield(data, 'A') && isfield(data, 'lab')
        A = data.A;
        labs = data.lab;
    else
        disp(['File ' matFiles(i).name ' does not contain variable "A" or "lab"']);
        continue;
    end
    
    % Convert table type to array (if applicable)
    if istable(A)
        A = table2array(A);
    end

    % Apply TopoLa and compute results
    A_mod = TopoLa(A, lambda); % Apply NR with fixed alpha
    y = louvain(A_mod); % Use existing function to compute y
    
    % Calculate NMI and ARI
    NMI_results(i) = Cal_NMI(y, labs); % Calculate NMI
    ARI_results(i) = Cal_ARI(y, labs, 'adjusted'); % Calculate ARI
end

% Create a table to store results
resultsTable = table();
for i = 1:length(matFiles)
    % Get the dataset name
    [~, dbname, ~] = fileparts(matFiles(i).name);
    
    % Store results in the table
    resultsTable = [resultsTable; table({dbname}, NMI_results(i), ARI_results(i), 'VariableNames', {'Dataset', 'NMI', 'ARI'})];
end

% Save the table to a CSV file
writetable(resultsTable, 'Results.csv');
fprintf('Results saved to Results.csv\n');

disp('All results have been saved to Results.csv.');
