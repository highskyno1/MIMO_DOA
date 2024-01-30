function s_pinv = DOA_PINV(u, scan_a)
% DOA_PINV 基于伪逆法实现稀疏恢复，效果等同于最大似然估计法
%   u       对回波自相关矩阵做酉对角化后，最大特征值对应的酉向量
%   scan_a  DOA估计的栅格，在稀疏恢复理论中也称为"超完备字典"
%   s_pinv  基于PINV法得到的不同来波方向的归一化功率

s_pinv = u' * pinv(scan_a)';
s_pinv = abs(s_pinv);
s_pinv = s_pinv ./ max(s_pinv);
end