function [Matrix] = TopoLa(A,lambda)

[n, m] = size(A);

temp=A;

C =(1/lambda*eye(size(temp'*temp))+temp' * temp)\temp' * temp;


Matrix =  C * A;

end
