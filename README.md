# MIMO阵列信号来向DOA估计实现~含FOCUSS、OMP、贝叶斯学习(SBL)等稀疏重构法和常规、子空间法、空间平滑滤波法
这是一篇CSDN博文的附录代码，文章地址：https://blog.csdn.net/weixin_41476562/article/details/135833238

仓库地址：https://github.com/highskyno1/MIMO_DOA

基于凸优化法（CVX）实现DOA时，==**依赖CVX工具箱**==，如果你的MATLAB没有安装，请前往[这里](http://cvxr.com/cvx/download/)下载，解压后在MATLAB命令行，cd到解压目录并执行其中的“cvx_setup.m”文件进行安装。如果不做CVX部分的仿真，可忽略这一步。
# 前言
波达方向估计(Direction Of Arrival, DOA)也称为测向、空间谱估计，为利用电磁波来获取目标或信源对天线阵列的角度信息，主要应用于雷达、通信、电子侦察与对抗等领域。
本文利用MIMO天线阵列实现DOA相关算法的总结，主要仿真实现了常规波束形成(CBF)、Capon和最大似然估计(ML)三种常规方法，多重信号分类法(MUSIC)、LS-ESPRIT和TLS-ESPRIT三种子空间方法，欠定系统聚焦法(FOCUSS)、正交匹配追踪法(OMP)、凸优化法(CVX)、伪逆法(PINV)和期望最大化-稀疏贝叶斯学习法(EM-SBL)等稀疏恢复方法。对比了上述方法在常规、低信噪比、低快拍以及信源相干情况下的性能，并研究了空间平滑算法在处理相干信源问题上的表现。
