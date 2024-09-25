clear;

% Initialize parameters
N = 500; % Number of nodes
E = 10000; % Number of edges

% Initialize an empty sparse adjacency matrix
A = sparse(N, N);

% Generate random edges without duplicates or self-loops
edges = randperm(N^2, 2*E); % Select 2E potential edge locations randomly
[i, j] = ind2sub([N, N], edges); % Convert linear indices to subscript indices
valid = i ~= j; % Exclude self-loops by ensuring no edge connects a node to itself
i = i(valid);
j = j(valid);
A(sub2ind([N, N], i, j)) = 1; % Add edges to the adjacency matrix
A(sub2ind([N, N], j, i)) = 1; % Ensure the graph is undirected by mirroring edges

% Alternatively, load a precomputed adjacency matrix
load("A.mat");

% Calculate the degree of each node
degrees = sum(A);

% Compute the number of common neighbors (BJ) for each pair of nodes
BJ = A * A; 

% Define the regularization parameter
alpha = 1e-3;

% Compute the topological dissimilarity (D_topola) using the regularization parameter
D_topola = (1/alpha*eye(size(A'*A)) + A' * A) \ (A' * A); 

% Calculate the denominator for Jaccard similarity
AI = repmat(degrees, N, 1) + repmat(degrees', 1, N) - BJ; 

% Calculate the Jaccard similarity for each pair of nodes
x_values = BJ ./ (AI - BJ); 

% Identify pairs of nodes with non-zero common neighbors
nonzero_common = BJ > 0;  

% Filter the calculated values based on having non-zero common neighbors
filtered_x_values = x_values(nonzero_common); 
filtered_D_topola = D_topola(nonzero_common); 

% Plot the relationship between topological similarity and D_topola for pairs with non-zero common neighbors
figure;
scatter(filtered_x_values, filtered_D_topola, 'filled');
xlabel('Topological similarity');
ylabel('D_{topo}');
title('Scatter plot with D_{topola} values for non-zero common neighbors');
grid on;
