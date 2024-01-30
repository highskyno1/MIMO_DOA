function P_ml = DOA_ML(scan_a, R_sig)
% DOA_ML 基于最大似然估计法得到DOA估计
%   scan_a  DOA估计栅格矩阵
%   R_sig   接收信号的自相关矩阵，阵元数*阵元数
%   P_ml    通过ML法得到的归一化估计结果

% 计算阵元数和栅格数
[rec_len,scan_len] = size(scan_a);
P_ml = zeros(1,scan_len);
for i = 1:scan_len
    scan = scan_a(:,i);
    Pa = scan / (scan'*scan) * scan';
    P_ml(i) = trace(Pa*R_sig) / rec_len;
end
P_ml = abs(P_ml);
P_ml = P_ml ./ max(P_ml);
end