function [Matrix] = TopoLa_self(A,p)

    [U,S,V] = svd(A,'econ');
    s = diag(S);
    s = s.^2;
    k = round(p * length(s))+1;
    
    threshold = s(k);
    S_=diag(S);
    
    for i=1:length(S_)
         if S_(i)~=0   
            S_(i)=S_(i)^3/(S_(i)^2+threshold);
         end
    end
    S(logical(eye(size(S))))=S_;
    
    Matrix=U*S*V';

   
    
    end