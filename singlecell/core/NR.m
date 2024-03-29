function [Matrix] = NR_EH(A,alpha)

[n, m] = size(A);

temp=A;

C =(1/alpha*eye(size(temp'*temp))+temp' * temp)\temp' * temp;


Matrix =  C * A;

end
