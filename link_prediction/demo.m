% Initialize parameters and dataset paths
n = 1;
alpha = 0.01;
beta = 0.1;
Path = {'USAir', 'NS', 'PB', 'Yeast', 'Celegans', 'Power', 'Ecoli', 'Fdataset', 'Cdataset'};
num_simulation = 1;
K = 10;  % K-fold cross-validation

% Setup
addpath('dataset', 'core');
rng('default');
load(Path{n});  % Load the dataset

% Prepare the network and index matrices
network = full(net);  % Ensure the network is a full matrix
M = triu(ones(size(network)), 1);  % Upper triangular matrix, exclude diagonal
all_index = find(M);  % Indices of the upper triangle

% Initialize arrays for metrics
aucArray = zeros(num_simulation, K);
auprArray = zeros(num_simulation, K);

% Simulation starts
for sim = 1:num_simulation
    fprintf('Simulation %d of %d\n', sim, num_simulation);
    
    % Cross-validation indices
    crossval_idx = crossvalind('Kfold', zeros(size(all_index)), K);
    
    % Loop over each fold for cross-validation
    for fold = 1:K
        fprintf('Fold %d of %d\n', fold, K);

        % Split data into training and test sets based on fold
        test_idx = all_index(crossval_idx == fold);
        train = network;
        train(test_idx) = 0;  % Set test positions to 0 in training data
        train = triu(train, 1) + triu(train, 1)';  % Ensure symmetry

        % Apply the predictive model
        predictedScores = TRWR(train, alpha, beta);  % Replace 'TRWR' with your model function

        % True labels and predicted scores for the test set
        testActual = network(test_idx);

        % Compute AUC and AUPR, ensuring all are vectors
        [X, Y, T, AUC] = perfcurve(testActual, predictedScores(test_idx), 1);
        [Xpr, Ypr, Tpr, AUPR] = perfcurve(testActual, predictedScores(test_idx), 1, 'xCrit', 'reca', 'yCrit', 'prec');
        
        % Store metrics
        aucArray(sim, fold) = AUC;
        auprArray(sim, fold) = AUPR;
    end
end

% Compute final average results across all simulations and folds
finalAUC = mean(aucArray, 'all');
finalAUPR = mean(auprArray, 'all');

% Display the results
fprintf('\nFinal Results:\n');
fprintf('Average AUC: %.4f\n', finalAUC);
fprintf('Average AUPR: %.4f\n', finalAUPR);
