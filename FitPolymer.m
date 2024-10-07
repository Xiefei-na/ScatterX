function [fitresult, gof] = FitPolymer(qc, logIc1)
%CREATEFIT(QC,LOGIC1)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      X 输入: qc
%      Y 输出: logIc1
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 20-May-2024 18:55:07 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData] = prepareCurveData( qc, logIc1 );

% 设置 fittype 和选项。
ft = fittype( 'log(  I0*(1./(1+c1.*q.^2)  )', 'independent', 'q', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';


% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% % 绘制数据拟合图。
% figure( 'Name', '无标题拟合 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'logIc1 vs. qc', '无标题拟合 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % 为坐标区加标签
% xlabel( 'qc', 'Interpreter', 'none' );
% ylabel( 'logIc1', 'Interpreter', 'none' );
% grid on


