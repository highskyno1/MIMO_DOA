function [x_sig, sigma, R_sig] = echo_generate(rec_len, d, lamda, src_angle, shot_len, f0, BW, Fg, snr, is_coherent)
% echo_generate 按照条件基于随机数滤波法生成回波矩阵，并计算噪声方差和回波的自相关矩阵
%   rec_len     阵元数量
%   d           阵元间隔
%   lamda       载波波长
%   src_angle   来波方向方向
%   shot_len    快拍数量
%   f0          载波中心频率
%   BW          信号带宽
%   Fg          系统全局采样率
%   snr         回波信号信噪比
%   x_sig       回波信号矩阵，阵元数*快拍数
%   sigma       估计的噪声方差
%   R_sig       回波的自相关系数矩阵
%   is_coherent 是否产生相干信号

% 计算来波目标数量
tar_len = length(src_angle);
% 计算目标导向矩阵
i = 0:rec_len-1;
A = exp(-1i*2*pi*d/lamda.*i'.*sin(src_angle));
% 生成信源
S = rand(tar_len,shot_len);
if is_coherent
    % 相干模式下，所有信源使用同一个频段
    BP = Filter_BP_IIR(Fg,f0-BW,f0-BW/2,f0+BW/2,f0+BW);
    S(1,:) = filter(BP,S(1,:),2);
    for i = 2:tar_len
        S(i,:) = S(1,:) * (abs(randn()) + 1);
    end
else
    % 非相干模式下，根据中心频率对信源做带通滤波(取代上变频)
    % 针对每个信源做IIR滤波，保证每个信源频谱基本不相交
    % 以保证信号信源间的正交性
    bw_each = linspace(f0-BW,f0+BW,tar_len);
    for i = 1:tar_len
        center = bw_each(i);
        % 获取IIR带通滤波器
        BP = Filter_BP_IIR(Fg,center-BW,center-BW/2,center+BW/2,center+BW);
        S(i,:) = filter(BP,S(i,:),2);
    end
end
% 计算回波
x_sig = A * S;
% 添加噪声
[x_sig,sigma] = awgn(x_sig,snr,"measured");
% 计算回波的自相关矩阵
R_sig = (x_sig * x_sig')./shot_len;
end