  function [ A] = Network_Enhancement(X,k)      

        % Assuming X is the matrix containing data points

        % Use NearestNeighbors to generate the neighborhood graph
        neigh = createns(X, 'nsmethod', 'kdtree');  % Create nearest neighbor search object
        k = 5;  % Set the number of neighbors
        [neighbors, dists] = knnsearch(neigh, X, 'K', k+1);  % Find the k+1 nearest neighbors for each data point (the first neighbor is the point itself)
        
        % Build the adjacency matrix (data type: double)
        n = size(X, 1);
        A = sparse(repmat((1:n)', [k, 1]), neighbors(:, 2:end), 1, n, n);
        A = max(A, A');  % Ensure the adjacency matrix is symmetric
        A = full(A);  % Convert the sparse matrix to a full matrix
        A = double(A);  % Convert the data type to double

  end