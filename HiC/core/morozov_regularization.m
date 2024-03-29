function [alpha_opt] = morozov_regularization(A, alpha_range)
% 计算A的奇异值分解


% 矩阵A的列数
n = size(A, 2);

% 初始化 eta
eta = zeros(length(alpha_range), 1);

% 遍历 alpha_range，计算 eta
for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    x_alpha = (A'*A + alpha^2 * eye(n)) \ (A'*A);  % 使用Tikhonov正则化求解反问题
    eta(i) = norm(alpha * x_alpha) / norm(A * x_alpha - A);  % 计算η(α)
end

% 选择最接近指定范围内的α
[~, idx] = min(abs(eta - 1));
alpha_opt = alpha_range(idx);

end