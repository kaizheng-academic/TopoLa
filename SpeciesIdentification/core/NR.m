function [Matrix] = LO(R_Omega,alpha)

[n, m] = size(R_Omega);
M = zeros(n, m);
M_recover =R_Omega;
temp=R_Omega;

W =(1/alpha*eye(size(temp'*temp))+temp' * temp)\temp' * temp;
Matrix = temp * W;

end
