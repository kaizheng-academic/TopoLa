function [Matrix] = DRN(R_Omega)

    [U,S,V] = svdsketch(R_Omega,0.001);
    
    
    
    
    
    S_=diag(S);
    k = round(0.001 * length(S));
    lambda = S_(k)*S_(k);
    
    for i=1:length(S_)
         if S_(i)~=0   
            S_(i)=S_(i)^3/(S_(i)^2+lambda);
         end
    end
    S(logical(eye(size(S))))=S_;
    
    Matrix=U*S*V';
    
    end
    