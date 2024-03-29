function [Matrix] = NR(A,alpha)
   

   
    [n, m] = size(A);
    M = zeros(n, m);
    temp=A;
    
    C =(1/alpha*eye(size(temp'*temp))+temp' * temp)\temp' * temp;
    
    Matrix =  C * A;

end


