function [Matrix] = WRWR_EH(A)
    
    [U,S,V] = svd(A,'econ');
    s = diag(S);
    s = s.^2;
    k = round(0.001 * length(s));
    threshold = s(k);
    
    if size(A,1) > size(A,2)
        C = (threshold*eye(size(A'*A)) + A' * A) \ A' * A;
        Matrix = C * A;
    else
        C = (threshold*eye(size(A*A')) + A * A') \ A * A';
        Matrix = A' * C;
    end
    
end