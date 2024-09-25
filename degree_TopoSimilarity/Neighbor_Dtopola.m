clear;

% Define the number of nodes and edges
N = 500; % Number of nodes
E = 10000; % Number of edges

% Initialize an empty sparse adjacency matrix
A = sparse(N, N);

% Generate random edges ensuring no duplicates and no self-loops
edges = randperm(N^2, 2*E); % Randomly select 2E potential edge locations
[i, j] = ind2sub([N, N], edges); % Convert linear indices to subscripts
valid = i ~= j; % Ensure no self-loops by excluding indices where i equals j
i = i(valid);
j = j(valid);
A(sub2ind([N, N], i, j)) = 1;
A(sub2ind([N, N], j, i)) = 1; % Make the graph undirected by mirroring edges

% Alternatively, load a precomputed adjacency matrix
load("A.mat");

% Calculate the degree of each node
degrees = sum(A);

% Calculate Jaccard similarity for each pair of nodes
BJ = A * A; % Number of common neighbors
AI = repmat(degrees, N, 1) + repmat(degrees', 1, N) - BJ; % Total neighbors minus common ones for Jaccard denominator
similarity = BJ ./ (AI - BJ); % Calculate Jaccard similarity
similarity(1:N+1:end) = 0; % Zero out diagonal to ignore self-similarity
similarity(isnan(similarity)) = 0; % Handle division by zero by setting NaNs to zero
similarity = full(similarity); % Convert from sparse to full matrix for visualization

alpha = 1e-3; % Regularization parameter
D_topola = (1/alpha*eye(size(A'*A)) + A' * A) \ (A' * A); % Calculate topological dissimilarity

% Plot the distribution of similarity values
figure;
histogram(similarity(:), 50, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
title('Distribution of Similarity Elements');
xlabel('Similarity Value');
ylabel('Frequency');
grid on;

% Calculate the sum of degrees minus common neighbors for each node pair
degree_sum_minus_common = bsxfun(@plus, degrees, degrees') - (A * A');

% Identify indices of edges within specified similarity ranges
indices_1 = similarity > 0.045 & similarity <= 0.055; % Range 1
indices_2 = similarity > 0.095 & similarity <= 0.105; % Range 2
indices_3 = similarity > 0.145 & similarity <= 0.155; % Range 3

% Extract topological dissimilarity and degree_sum_minus_common values for the identified ranges
D_topola_1 = D_topola(indices_1);
degree_sum_minus_common_1 = degree_sum_minus_common(indices_1);
D_topola_2 = D_topola(indices_2);
degree_sum_minus_common_2 = degree_sum_minus_common(indices_2);
D_topola_3 = D_topola(indices_3);
degree_sum_minus_common_3 = degree_sum_minus_common(indices_3);

% Plot the relationship between sum of degrees minus common neighbors and topological dissimilarity
figure;
hold on;
scatter(degree_sum_minus_common_1, D_topola_1, 'b.');
scatter(degree_sum_minus_common_2, D_topola_2, 'r.');
scatter(degree_sum_minus_common_3, D_topola_3, 'g.');
ylabel('D_{topo}');
xlabel('Degrees');
legend({'Similarity 5%', 'Similarity 10%', 'Similarity 15%'}, 'Location', 'Best');
hold off;
