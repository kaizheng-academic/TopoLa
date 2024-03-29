%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is a demo of the FNNM method on movie recommendation.
% MovieLens 100K is downloaded from [1].
%
% [1] F. M. Harper, J. A. Konstan. The movielens datasets: History and 
% context. Transactions on Interactive Intelligent Systems, 5(4):19, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [] = Demo_Movielens_FNNM_Topola()
%% STEP 1: Loading Data (MovieLens 100K)
clear
rng('default')
addpath('FNNM');
addpath('functions');
addpath('movielens_dataset');
load movielens_data
REPEAT = 5;
sr0 = 1;
for sr = 0.1 : 0.2 : 0.9;
    [Tdata, Tpos] = fsampling(E_org, sr); 
    [XX, R_A] = qr(X', 0);
    [YY, R_B] = qr(Y', 0);
    S2=corrcoef(XX');
    S1=corrcoef(YY');
    tol1 = 1e-4;
    maxiter = 3000;
    mu1 = 0.1;
    mu2 = mu1;
    
%% STEP 2: Choosing Parameters
    seq_lambda1 = [10];
    seq_lambda2 = [1e-1];
    seq_alpha =[1e-5 1e-4 1e-3  1e-2 1e-1  1 10];
    training_error_matrix = zeros(length(seq_alpha), length(seq_lambda1), length(seq_lambda2));
    for k = 1 :length(seq_alpha)
        for i = 1 : length(seq_lambda1 )
            for j = 1: length(seq_lambda2)
                [data_after_dump, dump_position] = selectdata(Tdata, 0.1);
                [M_ResultMat, A, B, iter] = FNNM_Topola(data_after_dump, XX, YY, seq_lambda1(i), seq_lambda2(j), mu1, mu2, maxiter, tol1,S1,S2,seq_alpha(k));
                [RMSE0, RelativeErr0] = fRMSE_new(E_org, M_ResultMat, Tpos);
                training_error_matrix(k,i, j) = RMSE0;
            end
        end
    end
    index = find(training_error_matrix == min(min(training_error_matrix)));
    [Ia,I_row, I_col] = ind2sub(size(training_error_matrix), index);
    lambda1 = seq_lambda1(I_row(1));
    lambda2 = seq_lambda2(I_col(1));
    alpha = seq_alpha(Ia(1))
%% STEP 3: Running for 5 times
    ts = 1;
    for times = 1 : REPEAT
        [M_ResultMat, A, B, iter] = FNNM_Topola(Tdata, XX, YY, lambda1, lambda2, mu1, mu2, maxiter, tol1,S1,S2,alpha);
        [RMSE, RelativeErr] = fRMSE_new(E_org, M_ResultMat, Tpos);
        result(ts, :) = [RMSE, RelativeErr, sr, lambda1, lambda2,alpha]
        ts = ts + 1;
        
        if times == REPEAT
            STmean = mean(result, 1)
            STstd = std(result, 1);
            VFmean(sr0, :) = [sr, STmean(1,1), STmean(1,2), lambda1, lambda2, mu1, mu2, iter, tol1, maxiter,alpha];
            VFstd(sr0, :) = [STstd(1,1), STstd(1,2)];
            break
        end
        
        [Tdata, Tpos] = fsampling(E_org, sr);
    end
    sr0 = sr0 + 1;
end

%% STEP 4: Saving Results
save MovieLens100K_FNNM_Topola VFmean VFstd
end
