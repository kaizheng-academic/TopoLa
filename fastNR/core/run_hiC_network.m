function [NMI] = run_hiC_network(databaseID)
    % Add all subdirectories of the current directory to the path
    addpath(genpath(pwd));

    % Load the raw Hi-C contact network and corresponding labels for a given resolution and database ID
    HiC_raw_resolution1K = load(['mat_', num2str(databaseID), '_1000.txt']); % Raw HiC network with 1K resolution
    label_resolution1K = load(['class_', num2str(databaseID), '_1000.txt']); % Node labels in HiC network with 1K resolution

    % Apply fastNR and NE
    HiC_NR = fastNR(HiC_raw_resolution1K,databaseID); % Apply fast Network Regularization
    HiC_NR = Network_Enhancement(HiC_NR);  % Apply Network Enhancement

    % Detect communities in the network using Louvain algorithm
    communities_NR = louvain(HiC_NR);

    % Calculate and return the Normalized Mutual Information (NMI) between the detected communities and true labels
    NMI = nmi(communities_NR, label_resolution1K);
end
