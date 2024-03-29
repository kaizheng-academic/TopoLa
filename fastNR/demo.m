% Define the range of database IDs
databaseIDs = 1:22;

% Initialize a table to store the database IDs and corresponding NMI values
results = table([], [], 'VariableNames', {'DatabaseID', 'NMI'});

% Iterate over the range of database IDs
for dbID = databaseIDs
    % Call your function with the current database ID and record the NMI value
    current_NMI = run_hiC_network(dbID);
    
    % Add the results for the current database ID to the results table
    newRow = {dbID, current_NMI};
    results = [results; newRow];
end

% After collecting all results, save them to a MAT file
save('fastNR.mat', 'results');
