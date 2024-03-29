function [Matrix] = WRWR(R_Omega,alpha,beta)

[n, m] = size(R_Omega);
M = zeros(n, m);
temp=R_Omega;

C =(1/alpha*eye(size(temp'*temp))+temp' * temp)\temp' * temp;
W= R_Omega ./ sum(R_Omega, 2); 
W(isnan(W)) = 0;

Q =(1-beta)*inv(eye(size(W))-beta*W);

Matrix = Q * C * R_Omega;

end
