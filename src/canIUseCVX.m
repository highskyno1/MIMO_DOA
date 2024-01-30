function res = canIUseCVX()
try 
    cvx_where;
    res = true;
catch
    warning('!您的电脑没有安装CVX工具箱，请前往：http://cvxr.com/cvx/download/ 下载安装!');
    warning('或者你可以注释掉DOA_CVX与canIUseCVX两个函数的调用以关闭本警告');
    res = false;
end
end