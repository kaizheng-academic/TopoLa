% Define the range of ii values
ii_values = 1:22;

% Define the alpha values for iteration
alpha_NR = [100, 1, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9];

% Initialize a table to store the dataset indices, best NMI values, and corresponding alpha values
results = table([], [], [], 'VariableNames', {'DatasetIndex', 'BestNMI', 'BestAlpha'});

% Iterate over the range of ii values
for ii = ii_values
    % Initialize variables to keep track of the best NMI and alpha for the current ii
    best_NMI_for_ii = -inf; % Start with a very small number
    best_alpha_for_ii = NaN; % Placeholder for the best alpha
    
    % Iterate over the range of alpha values
    for alpha = alpha_NR
        % Call your function with the current ii and alpha, and get the NMI value
        current_NMI = run_hiC_network(ii, alpha);
        
        % Check if the current NMI is better than the best NMI for the current ii
        if current_NMI > best_NMI_for_ii
            best_NMI_for_ii = current_NMI; % Update the best NMI
            best_alpha_for_ii = alpha; % Update the best alpha
        end
    end
    
    % Add the results for the current ii to the results table
    newRow = {ii, best_NMI_for_ii, best_alpha_for_ii};
    results = [results; newRow];
end

% After collecting all results, save them to a MAT file
save('NRND.mat', 'results');
