function [Matrix] = DRN(A,lambda)

    [U,S,V] = svdsketch(A,1e-1);
    
    
    
    
    
    S_=diag(S);
    for i=1:length(S_)
         if S_(i)~=0   
            S_(i)=S_(i)^3/(S_(i)^2+1/lambda);
         end
    end
    S(logical(eye(size(S))))=S_;
    
    Matrix=U*S*V';

   
    
    end
    