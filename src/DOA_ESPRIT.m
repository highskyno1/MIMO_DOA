function [DOA_esp_ml, DOA_esp_tls] = DOA_ESPRIT(x_sig, target_len, lamda, d)
% DOA_ESPRIT 基于旋转不变子空间法实现DOA
%   x_sig       每个阵元接收到的信号矩阵，阵元数*快拍数
%   target_len  目标数量
%   lamda       载波波长
%   d           阵元间隔
%   DOA_esp_ml  基于最大似然估计准则得到的估计结果
%   DOA_esp_tls 基于最小二乘准则得到的估计结果

% 计算阵元数
rec_len = size(x_sig,1);
% 回波子阵列合并
x_esp = [x_sig(1:rec_len-1,:);x_sig(2:rec_len,:)];
% 计算协方差
R_esp = cov(x_esp');
% 特征分解
[~,~,W] = eig(R_esp);
% 获取信号子空间
W = fliplr(W);
U_s = W(:,1:target_len);
% 拆分
U_s1 = U_s(1:rec_len-1,:);
U_s2 = U_s(rec_len:end,:);

%% LS-ESPRIT法
mat_esp_ml = pinv(U_s1) * U_s2;
% 获取对角线元素并解算来向角
DOA_esp_ml = -angle(eig(mat_esp_ml));
DOA_esp_ml = asin(DOA_esp_ml .* lamda ./ 2 ./ pi ./ d);
DOA_esp_ml = rad2deg(DOA_esp_ml);

%% TLS-ESPRIT
V = [U_s1,U_s2];
[~,~,V] = svd(V);
% 提取E12和E22
E12 = V(1:target_len,target_len+1:end);
E22 = V(target_len+1:end,target_len+1:end);
mat_esp_tls = - E12 / E22;
% 获取对角线元素并解算来向角
DOA_esp_tls = -angle(eig(mat_esp_tls));
DOA_esp_tls = asin(DOA_esp_tls .* lamda ./ 2 ./ pi ./ d);
DOA_esp_tls = rad2deg(DOA_esp_tls);
end