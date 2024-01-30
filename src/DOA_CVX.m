function s_cvx = DOA_CVX(u, scan_a, p_norm, tor_lim)
% DOA_CVX 尝试利用凸优化方法实现稀疏恢复，获得DOA估计
%   !!使用前必须先安装CVX凸优化工具箱，下载地址为：http://cvxr.com/cvx/download/
%   u       对回波自相关矩阵做酉对角化后，最大特征值对应的酉向量
%   scan_a  DOA估计的栅格，在稀疏恢复理论中也称为"超完备字典"
%   s_cvx   基于CVX法得到的不同方向来波功率估计
%   p_norm  约束结果的P范数
%   tor_lim 对稀疏恢复整体(u - scan_a * s_cvx)的容许限度
%           如果结果为+inf，说明发散了，请增大该值
% 注意在CVX工具箱的scope中
% 尽量调用工具箱自己的函数（定义在CVX解压目录/functions 文件夹下）
% 使用MATLAB定义的函数可能会有问题

% 计算DOA栅格数量
scan_len = size(scan_a,2);
%% CVX工具箱调用语法~开始
cvx_begin
variable s_cvx(scan_len,1)
minimize(  sum(pow_abs(s_cvx, p_norm)) );
subject to
    sum(pow_abs(u - scan_a * s_cvx,2)) <= tor_lim;
cvx_end
% CVX工具箱调用语法~结束
s_cvx = abs(s_cvx);
s_cvx = s_cvx ./ max(s_cvx);
end