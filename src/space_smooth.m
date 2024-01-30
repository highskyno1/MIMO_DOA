function R_sig = space_smooth(rec_len, equ_l, x_sig)
% 基于空间平滑法解决相干信源问题
%   rec_len 真实阵元数量
%   equ_l   等效阵元数量
%   x_sig   回波矩阵，真实阵元数量*快拍数

% 计算平滑阵元数量
rec_len_equ = rec_len - equ_l + 1;
% 计算子阵的自相关
R_sig = zeros(equ_l,equ_l);
for i = 1:rec_len_equ
    % 第i个子阵
    foo = x_sig(i:equ_l+i-1,:);
    % 计算协方差
    R_sig = R_sig + (foo * foo')./rec_len;
end
R_sig = R_sig ./ rec_len_equ;
end