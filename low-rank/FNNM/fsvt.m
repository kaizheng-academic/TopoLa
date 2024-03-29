function Y = fsvt(X, a)
%% This is a singular value thresholding operator for matrix 'X 'by thretholding parameter 'a'.
[S, V, D] = svd(X, 'econ');
v = diag(V);
[V_row, V_col] = size(V);
a = a * ones(size(v));
v_new = zeros(size(v));
nonZero = v > a;
v_new(nonZero) = v(nonZero) - a(nonZero);
if V_row < V_col
    Y = S * [diag(v_new), zeros(V_row, V_col - V_row)] * D';
else Y = S * [diag(v_new); zeros(V_row - V_col, V_col)] * D';
end
end