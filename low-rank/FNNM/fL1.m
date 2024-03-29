function Y = fL1(X, a)
%% This is a L1 operator for matrix 'X' by thretholding parameter 'a'.
[n1, n2] = size(X);
for i = 1 : n1
    for j = 1 : n2
        if X(i, j) > a
            X(i, j) = X(i, j) - a;
        elseif X(i, j) <-a
            X(i, j) = X(i, j) + a;
        else
            X(i, j) = 0;
        end
    end
end
Y = X;
end