
function Gaussian_Dis = calculate_gaussian(A)
    n = size(A, 1); % 获取矩阵的行数
    Gaussian_Dis = zeros(n, n); % 创建大小为 n × n 的全零矩阵
   
    sum = 0; % 记录总数
    

    % 按行来求（有按列来的）
    for i = 1:n
        temp = norm(A(i,:));
        sum = sum + temp^2;
    end
    pare_a = 1 / (sum / n);

    for i = 1:n
        for j = 1:n
            Gaussian_Dis(i,j) = exp(-pare_a * (norm(A(i,:) - A(j,:))^2));
        end
    end
end


