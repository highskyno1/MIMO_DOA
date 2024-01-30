%{
    CBF、Capon、MUSIC、LS-ESPRIT、TLS-ESPRIT、ML法、FOCUSS、OMP、CVX、伪逆和SBL结果比较，
    2024/01/25：第一版
%}
close all;

%% 参数定义区
c = physconst('LightSpeed');    % 光速 m/s
f0 = 10e9;                      % 载波频率 Hz
fg = 50e9;                      % 系统全局采样率
BW = 1e8;                       % 信号带宽
src_angle = deg2rad([-10,-20]); % 目标来向角 °
M = 32;                         % 阵元数量
L = 128;                        % 快拍数
snr = 20;                       % 信噪比 dB
scan_angle = deg2rad(-60:0.1:60);   % 扫描角度
is_space_smooth = false;        % 是否使用空间平滑法
is_coherent = false;            % 是否产生相干的信源

%% 参数计算区
lambda = c / f0;                % 载波波长 m*
d = lambda / 2;                 % 阵元间隔 m
P = length(src_angle);          % 目标数

%% 回波生成
[x_sig, sigma, R_sig] = echo_generate(M, d, lambda, src_angle, L, f0, BW, fg, snr, is_coherent);

%% 空间平滑法
if is_space_smooth
    equ_l = 14; % 等效阵元数量，必须小于真实阵元数 @rec_len
    % 进行空间平滑滤波以破坏信源的相干性
    R_sig = space_smooth(M, equ_l, x_sig);
    % 将阵元数改变为等效阵元数
    M = equ_l;
end

%% 扫描算法通用区
% 计算波束形成矩阵（DOA栅格矩阵，在稀疏恢复中也称为"完备字典"）
i = 0:M-1;
a = exp(-1i*2*pi*d/lambda.*i'.*sin(scan_angle));
% 计算Rxx的主特征矢量
[U,~] = eig(R_sig);
u = U(:,end);
% 判断用户是否安装CVX工具箱
is_cvx_access = canIUseCVX();

%% CBF法
P_cbf = DOA_CBF(a, R_sig);

%% Capon法
P_capon = DOA_Capon(a, R_sig);

%% MUSIC法
P_music = DOA_MUSIC(U, P, a);

%% ESPRIT法
% 分别返回基于ML准则与TLS准则的ESPRIT估计结果
% 实测这两个准则下的结果相差不大
[DOA_esp_ml, DOA_esp_tls] = DOA_ESPRIT(x_sig, P, lambda, d);

%% ML法
P_ml = DOA_ML(a, R_sig);

%% FOCUSS法
% 设定迭代条件
lamda_reg = 1e-4;   % 正则化因子
lamda_err = 1e-4;   % 迭代结束误差
lamda_spe = 0;      % 稀疏因子
P_focuss = DOA_FOCUSS(a, u, lamda_spe, lamda_reg, lamda_err);

%% 尝试OMP法
% 设置退出条件
omp_toler = 1e-4;
P_omp = DOA_OMP(u, a, omp_toler);

%% 尝试伪逆法
P_pinv = DOA_PINV(u, a);

%% 尝试EM-SBL法
err_lim = 1e-3;             % 误差限
times_lim = 30;             % 迭代次数限制
s_sbl = DOA_EM_SBL(sigma, a, R_sig, L, err_lim, times_lim);

%% 尝试CVX法
if is_cvx_access
    tor_lim = 1e-1; % CVX的误差容许限，结果发散时请增大该值
    % 分别尝试P1范数、P1.5范数与P2范数约束下的结果
    P_cvx_1 = DOA_CVX(u, a, 1, tor_lim);
    P_cvx_1_5 = DOA_CVX(u, a, 1.5, tor_lim);
    P_cvx_2 = DOA_CVX(u, a, 2, tor_lim);
end

%% 绘制方位角对比图
plot_x = rad2deg(scan_angle);
figure;
if is_cvx_access
    subplot(311);
else
    subplot(211);
end
plot(plot_x,10*log10(P_cbf),'LineWidth',1);
hold on
plot(plot_x,10*log10(P_capon),'LineWidth',1);
plot(plot_x,10*log10(P_music),'LineWidth',1);
plot(plot_x,10*log10(P_focuss),'LineWidth',1.5);
plot(plot_x,10*log10(s_sbl),'LineWidth',1.5);
scatter(DOA_esp_ml,zeros(1,length(DOA_esp_ml)),60,'filled');
xline(rad2deg(src_angle),'--','真实来波方向');
legend('CBF','Capon','MUSIC','FOCUSS','EM-SBL','LS-ESPRIT',Location='southeast');
xlabel('目标信号来向角度/°');
ylabel('归一化强度/dB');
if is_cvx_access
    subplot(312);
else
    subplot(212);
end
plot(plot_x,10*log10(P_omp),'LineWidth',1.5);
hold on
plot(plot_x,10*log10(P_pinv),'LineWidth',1.5);
plot(plot_x,10*log10(P_ml),'LineWidth',1.5);
scatter(DOA_esp_tls,zeros(1,length(DOA_esp_tls)),60,'filled');
xline(rad2deg(src_angle),'--','真实来波方向');
legend('OMP','PINV','ML','TLS-ESPRIT',Location='southeast');
xlabel('目标信号来向角度/°');
ylabel('归一化强度/dB');
if is_cvx_access
    subplot(313);
    plot(plot_x,10*log10(P_cvx_1),'LineWidth',1);
    hold on;
    plot(plot_x,10*log10(P_cvx_1_5),'LineWidth',1);
    plot(plot_x,10*log10(P_cvx_2),'LineWidth',1);
    xline(rad2deg(src_angle),'--','真实来波方向');
    legend('CSV_{P1}','CSV_{P1.5}','CSV_{P2}',Location='southeast');
    xlabel('目标信号来向角度/°');
    ylabel('归一化强度/dB');
end
