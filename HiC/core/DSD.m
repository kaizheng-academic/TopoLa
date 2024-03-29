function [Matrix] = DSD(A)


    n = length(A);
    degree = sum(A,1);
    p = A / degree;
    pi = degree/sum(degree);
    
    p_=repmat(p, 1, n);
    pi_=repmat(pi', 1, n);
    Matrix = squareform(pdist(inv(eye(n) - p_ - pi_),'cityblock'));
    Matrix = Rbf_kernel(Matrix);
    





