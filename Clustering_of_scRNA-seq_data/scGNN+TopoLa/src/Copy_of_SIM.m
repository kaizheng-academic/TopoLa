  function [S] = SIM(X)      

        D_Kernels = multipleK(X);
        distX = mean(D_Kernels,3);
        S0 = max(max(distX))-distX;
        S0 = Network_Diffusion(S0,k);
        S0 = NE_dn(S0,'ave');
        S= (S0 + S0')/2;

  end

function D_Kernels = multipleK(x)


N = size(x,1);
KK = 0;
sigma = [2:-0.25:1];
Diff = (dist2(x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
allk = 10:2:30;
t=1;
for l = 1:length(allk)
    if allk(l) < (size(x,1)-1)
        TT=mean(T(:,2:(allk(l)+1)),2)+eps;
        Sig=(repmat(TT,1,n)+repmat(TT',n,1))/2;
        Sig=Sig.*(Sig>eps)+eps;
        for j = 1:length(sigma)
            W=normpdf(Diff,0,sigma(j)*Sig);
            Kernels(:,:,KK+t) = (W + W')/2;
            t = t+1;
        end
    end
end

for i = 1:size(Kernels,3)
    K = Kernels(:,:,i);
    k = 1./sqrt(diag(K)+1);
    %G = K.*(k*k');
    G = K;
    D_Kernels(:,:,i) = (repmat(diag(G),1,length(G)) +repmat(diag(G)',length(G),1) - 2*G)/2;
    D_Kernels(:,:,i) = D_Kernels(:,:,i) - diag(diag(D_Kernels(:,:,i)));
  
end

end

function W = Network_Diffusion(A, K)
%K = min(2*K, round(length(A)/10));
A = A-diag(diag(A));
P = (dominateset(double(abs(A)),min(K,length(A)-1))).*sign(A);
DD = sum(abs(P'));
P = P + (eye(length(P))+diag(sum(abs(P'))));
P = (TransitionFields(P));
[U,D] = eig(P);
d = real((diag(D))+eps);
alpha = 0.8;
beta = 2;
d = (1-alpha)*d./(1-alpha*d.^beta);


D = diag(real(d));
W = U*D*U';

W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),DD);
W=D*(W);
W = (W+W')/2;

end


function wn=NE_dn(w,type);
w = w*length(w);
w = double(w);
D=sum(abs(w),2)+eps;

if type == 'ave'
    D=1./D;
    D=sparse(1:length(D),1:length(D),D);
    wn=D*w;
elseif type == 'gph'
    D=1./sqrt(D);
    D=sparse(1:length(D),1:length(D),D);
    wn=D*(w*D);
end
end