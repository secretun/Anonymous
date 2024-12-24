1.主函数为fastica,函数定义为：
function [Out1, Out2, Out3] = fastica(mixedsig, varargin)

- mixedsig: 观测混合信号，M*K列表示有 M路观测信号，共有K个采样点。
- varargin: 可调参数，具体可以参照《Matlab中FastICA工具箱的使用说明.pdf》

- Out1: icasig,解混合后的信号，每一行表示估计独立分量（IC）。
- Out2：A，计算出来的混合矩阵
- Out3: W，计算出来的解混合矩阵


2.使用样例（MixedS.mat为测试数据）：
eg1: [icasig, A, W] = fastica(MixedS)

eg2: [icasig] = fastica(MixedS)

eg3: [A, W] = fastica(MixedS)

eg4: % 'numOfIC'为独立分量数，默认和输入混合信号行数相同，'displayMode'设置绘图模式。

     [icasig, A, W] = fastica(MixedS, 'numOfIC', 3, 'displayMode', 'signals')


3.具体说明可以参见《Matlab中FastICA工具箱的使用说明.pdf》


4.fasticag提供可视化操作，但是个人实践发现无法导入数据，所以无法正常使用。
