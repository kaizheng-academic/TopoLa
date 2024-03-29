function [B, C] = fsampling(A, sr)
% A: Data
% sr: Sampling ratio
% B: Missing data
% C: The position matrix of missing entries

posA = double(A ~= 0);
num_nonzero = sum(sum(A ~= 0));
num_remove = floor(num_nonzero * (1 - sr));
pos_nonzero = find(A ~= 0);
remove = datasample(pos_nonzero, num_remove, 'Replace', false);
A(remove) = 0;
B = A;
posB = double(B ~= 0);
C = posA - posB;
end