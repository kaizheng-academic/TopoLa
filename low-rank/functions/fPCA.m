function  U_dr = fPCA(A, ratio)
pca_l = 0;
singular_sum = svd(A);
for i = 1 : rank(A)
    if sum(singular_sum(1 : i)) > ratio * sum(singular_sum)
        pca_l = i;
        break;
    end
end
[U_dr, S_dr, V_dr] = svds(A, pca_l);
end