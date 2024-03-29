function [Matrix] = fastNR(R_Omega, databaseID)
    % Perform a randomized SVD on the input matrix R_Omega using the randQB_fp algorithm.
    [U, S, V] = svdsketch(R_Omega, 0.001); % The second parameter here controls the accuracy of the approximation.
    
    % Set the default parameter to 0.001.
    parameter = 0.001; 
    if databaseID == 2
        parameter = 0.01; % Adjust the parameter for database ID 2
    end

    % Extract the diagonal singular values from S and determine the lambda based on the parameter.
    S_values = diag(S);
    k = round(parameter * length(S_values)); % Calculate the index of the singular value.
    lambda = S_values(k)^2; % Lambda is determined by squaring the k-th singular value.
    
    % Apply NR to the singular values.
    for i = 1:length(S_values)
        if S_values(i) ~= 0   
            S_values(i) = S_values(i)^3 / (S_values(i)^2 + lambda); % Regularize each non-zero singular value.
        end
    end
    
    % Update the S matrix with the new singular values.
    S = diag(S_values);
    
    % Reconstruct the matrix.
    Matrix = U * S * V';
end