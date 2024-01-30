function P_music = DOA_MUSIC(U, target_len, scan_a)
% DOA_MUSIC 基于多重信号分类法实现DOA估计
%   U           对接收信号的自相关矩阵做酉对角化分解后的左酉矩阵
%   target_len  目标数量
%   rec_len     阵元数量
%   scan_a      DOA估计栅格

% 计算阵元数和栅格数
[rec_len,scan_len] = size(scan_a);
W = fliplr(U);
% 已知只有两个信源，用信号子空间法
U_s = W(:,1:target_len);
U_music = eye(rec_len) - U_s * U_s';
P_music = zeros(1,scan_len);
for i = 1:scan_len
    foo = scan_a(:,i);
    P_music(i) = 1 / (foo' * U_music * foo);
end
P_music = abs(P_music);
P_music = P_music ./ max(P_music);
end