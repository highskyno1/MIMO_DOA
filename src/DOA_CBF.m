function P_cbf = DOA_CBF(scan_a, R)
% DOA_CBF 基于常规波束形成法实现DOA估计
%   scan_a  需要估计的来波方向
%   R       快拍的协方差矩阵，L*L维，L为快拍数
%   P_cbf   各个scan_a栅格上的来波归一化强度

scan_len = size(scan_a,2);
P_cbf = zeros(1,scan_len);
for i = 1:scan_len
    foo = scan_a(:,i);
    P_cbf(i) = foo' * R * foo / (foo' * foo)^2;
end
P_cbf = abs(P_cbf);
P_cbf = P_cbf ./ max(P_cbf);
end