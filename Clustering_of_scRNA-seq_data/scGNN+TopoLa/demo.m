clear
addpath('result_scGNN') % Add path for result_scGNN folder
addpath('data') % Add path for data folder
addpath('src') % Add path for function folder

% Get all .mat files in the result_scGNN folder
mat_files = dir('result_scGNN/*.mat');
selected_files = {mat_files.name};

dataset = selected_files; % Store dataset file names


NMI_values = []; % Create an empty array to store NMI values
ARI_values = []; % Create an empty array to store ARI values

% Initialize a results table
results = table([], [], [], 'VariableNames', {'Dataset', 'NMI', 'ARI'});

for i = 1:length(dataset)
    % Perform analysis on the current dataset
    load(['data/', dataset{i}]); % Load data
    datasetID = dataset{i}(1:end-4); % Extract dataset ID
    load(['result_scGNN/', dataset{i}]); % Load results for the current dataset

    % Generate the graph using TopoLa
    X_graph = TopoLa(U, S, V);
    C = length(unique(labs)); % Number of clusters
    if length(X_graph) < length(labs)
        X_graph_ = zeros(length(labs), length(labs)); % Initialize an empty graph
        X_graph_(1:length(X_graph), 1:length(X_graph)) = X_graph; % Fill in the graph
        X_graph = X_graph_; % Update the graph
    end
    
    % Perform community detection using Louvain algorithm
    com = louvain(X_graph);
    NMI = nmi(com, labs); % Calculate NMI
    ARI = Cal_ARI(com, labs); % Calculate ARI
    NMI_values = [NMI_values; NMI]; % Add NMI value to the array
    ARI_values = [ARI_values; ARI]; % Add ARI value to the array

    % Append the new row of results
    results = [results; {datasetID, NMI, ARI}];
end

% Save the results to a CSV file
writetable(results, 'Result.csv');
