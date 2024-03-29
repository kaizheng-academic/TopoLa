function [M_recover, M, N, iter] = FNNM(R_Omega, X, Y, lambda1, lambda2, mu1, mu2, maxiter, tol,dc,dg,alpha)
%% FNNM: Feature and nuclear norm minimization for matrix completion.
% Usage:  [M_recover, M, N, iter] = FNNM(R_Omega, X, Y, lambda1, lambda2, mu1, mu2, maxiter, tol)
%
% Inputs:
%        R_Omega                 - a matrix recording the observed positions in the target matrix.
%        lambda1, lambda2        - parameters needed to give.
%        X, Y                    - they are features or side information.
%        mu1, mu2                - they both are set to be 0.1.
%        maxiter                 - maximum number of iterations.
%        tol                     - tolerance of termination conditions.
%
% Outputs:
%        M_recover               - the completed matrix.
%        M, N                    - M is a low-rank matrix and N is a sparse matrix.
%        iter                    - the number of iterations.


%% Step1: Initialization
[n, m] = size(R_Omega);
r_a = size(X, 2);
r_b = size(Y, 2);
M = zeros(n, m);
N = zeros(r_a, r_b);


%% Step2: Iteration
i = 1;
stop1 = 1;
matr = M + X * N * Y';
while(stop1 > tol )
    stop1_1 = stop1;
    matr_1 = matr;
    % the process of getting U
    U = M - mu1 * (M + X * N * Y' - R_Omega) .* double(R_Omega ~= 0);
    % the process of getting M
    M_1 = fsvt(U, lambda1 * mu1);
    % the process of getting V
    V = N - mu2 * X' * ((M_1 + X * N * Y' - R_Omega) .* double(R_Omega ~= 0)) * Y;
    % the process of getting N
    N_1 = fL1(V, lambda2 * mu2);
    
    matr = M_1 + X * N_1 * Y';
    stop1 = norm(matr - matr_1, 'fro') / norm(matr_1, 'fro');
    M = M_1;
    N = N_1;
    i=i+1;
    
    if i < maxiter
        iter = i - 1;
    else
        iter = maxiter;
        break
    end
end
%% Step1: Initialization



%% Step3: The completed matrix
M_recover = X * N * Y' + M;

temp=[dg,M_recover;M_recover',dc]; % Integration matrix \mathbf{A}

W = lsqminnorm( (1/alpha*eye(size(temp'*temp)) + temp' * temp),temp')* temp;
%ind_nan = find(isnan(W));
%W = (1/alpha*eye(size(temp'*temp)) + temp' * temp)\temp' * temp; % Computing the weighting matrix
Matrix = temp * W;
M_recover = Matrix(1:n,n+1:end); 
end
