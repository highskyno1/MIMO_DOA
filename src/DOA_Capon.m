function P_capon = DOA_Capon(scan_a, R)
% DOA_Capon 基于Capon法实现DOA估计
%   scan_a  需要估计的来波方向
%   R       快拍的协方差矩阵，L*L维，L为快拍数
%   P_capon 各个scan_a栅格上的来波归一化强度

scan_len = size(scan_a,2);
P_capon = zeros(1,scan_len);
for i = 1:scan_len
    foo = scan_a(:,i);
    P_capon(i) = 1 / (foo' * pinv(R) * foo);
end
P_capon = abs(P_capon);
P_capon = P_capon ./ max(P_capon);
end