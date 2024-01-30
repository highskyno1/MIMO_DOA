function s_sbl = DOA_EM_SBL(sigma, scan_a, R_sig, shot_len, err_lim, times_lim)
% 基于期望最大化-稀疏贝叶斯学习方法实现DOA估计
%   sigma       估计的噪声方差
%   scan_a      DOA估计的栅格，在稀疏恢复理论中也称为"超完备字典"
%   R_sig       回波自相关矩阵
%   shot_len    快拍数
%   err_lim     迭代误差限，迭代退出的条件之一
%   times_lim   迭代次数限制，迭代退出的条件之二

% 计算DOA栅格数量
scan_len = size(scan_a,2);
times_cnt = 0;
Gamma = eye(scan_len)*0.1;  % 初始Gamma矩阵
while true 
    times_cnt = times_cnt + 1;
    % E-step
    Sigma_x = pinv(sigma * (scan_a'*scan_a) + pinv(Gamma));
    Mu_x = Sigma_x / sigma * scan_a' * R_sig;
    % M-step
    Gamma_new = Gamma;
    for i = 1:scan_len
        mu_xn = Mu_x(i,:); 
        Gamma_new(i,i) = mu_xn*mu_xn'/shot_len + Sigma_x(i,i);
    end
    if sum(abs(diag(Gamma_new - Gamma))) < err_lim || times_cnt > times_lim
        break;
    end
    Gamma = Gamma_new;
end
Gamma_new = abs(diag(Gamma_new));
s_sbl = Gamma_new ./ max(Gamma_new);
end