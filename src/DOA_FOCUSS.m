function P_focuss = DOA_FOCUSS(scan_a, u, lamda_spe, lamda_reg, lamda_err)
% DOA_FOCUSS 基于欠定系统局灶解法(Focal Under determinedSystem Solver)
% 实现稀疏恢复获得DOA估计结果
%   scan_a      DOA估计的栅格，在稀疏恢复理论中也称为"超完备字典"
%   u           对回波自相关矩阵做酉对角化后，最大特征值对应的酉向量
%   lamda_spe   稀疏因子，效果类似于结果的范数约束
%   lamda_reg   正则化因子，过大会趋于0解，过小结果发散
%   lamda_err   迭代结束误差
%   P_focuss    通过FOCUSS法得到的归一化来波方向功率估计

% 计算阵元数
rec_len = size(u,1);
% 确定s的初始值
Dg = scan_a;
s0 = Dg' / (Dg * Dg') * u;
for i = 1:1000
    W = diag(s0.^(1-lamda_spe/2));
    s = W * W' * Dg' / (Dg * (W * W') * Dg' + lamda_reg .* eye(rec_len)) * u;
    if norm(s - s0,2) / norm(s0,2) < lamda_err
        break;
    end
    s0 = s;
end
P_focuss = abs(s);
% 数据饱和钳制
P_focuss = P_focuss ./ max(P_focuss);
P_focuss(P_focuss < 1e-4) = 1e-4;
end