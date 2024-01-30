function P_omp = DOA_OMP(u, scan_a, omp_toler)
% 基于正交匹配追踪法实现来波的DOA估计
%   u           对回波自相关矩阵做酉对角化后，最大特征值对应的酉向量
%   scan_a      DOA估计的栅格，在稀疏恢复理论中也称为"超完备字典"
%   omp_toler   容许的最小误差，低于此误差时结束迭代
%   P_omp       基于OMP法得到的不同方向来波功率估计

% 计算阵元数量和栅格数量
[rec_len,scan_len] = size(scan_a);
% 初始化标签集合
omp_omiga_mask = zeros(1,scan_len,'logical');
% 初始化初始向量
omp_r = u;
for i = 1:scan_len
    % 求字典矩阵中与残差向量rk-1最强相关的原子
    foo = abs(omp_r' * scan_a);
    [~,idx] = max(foo);
    % 往omiga中添加该列索引
    omp_omiga_mask(idx) = true;
    % 依据索引从字典中选出需要的列向量
    omp_omiga_foo = scan_a(:,omp_omiga_mask);
    % 奇异判断
    rms = det(omp_omiga_foo' * omp_omiga_foo);
    % 退出条件
    if rms < omp_toler
        break;
    end
    % 最小化min(x) ||y-fa*x|| L2范数的平方最小
    omp_xk = (omp_omiga_foo' * omp_omiga_foo) \ omp_omiga_foo' * u;
    % 更新残差
    omp_r = (eye(rec_len) - omp_omiga_foo / (omp_omiga_foo.'*omp_omiga_foo) * omp_omiga_foo.') * u;
end
P_omp = zeros(1,scan_len);
% 提取结果
foo = find(omp_omiga_mask==true);
P_omp(foo(1:length(omp_xk))) = abs(omp_xk)';
% 归一化
P_omp = P_omp ./ max(P_omp);
P_omp(P_omp < 1e-4) = 1e-4;
end