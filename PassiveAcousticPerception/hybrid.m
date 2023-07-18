figure
hold on
axis equal

% 设置图形字体和网格线
set(gca, 'FontName', 'Times New Roman');
grid on;
set(gca,'GridColor',[0.8 0.8 0.8],'GridAlpha',0.5, 'LineWidth', 1.5)

% 画前10个椭圆，a的间距是0.1
for i = 0:10
    a = 2 + i*0.2;  % a从1开始，增长0.1
    b = a/3;  % 保持a/b的比例为2:1
    r = b;
    rectangle('Position',[-r, -r, 2*r, 2*r], 'Curvature',[1 1], 'EdgeColor', 'r')  % 圆心在(0,0)

    % 创建椭圆的参数方程
    t = linspace(0, 2*pi, 1000);
    x = a*cos(t) + 1;  % 椭圆中心点均位移2单位
    y = b*sin(t);

    plot(x, y, 'b')  % 'b'表示黑色
end
% 
% % 画后10个椭圆，a的间距是0.2
% for i = 0:9
%     a = 2.45 + i*0.06;  % a从2开始，增长0.2
%     b = a/3;  % 保持a/b的比例为2:1
%     r = b;
%     rectangle('Position',[-r, -r, 2*r, 2*r], 'Curvature',[1 1], 'EdgeColor', 'r')  % 圆心在(0,0)
%     % 创建椭圆的参数方程
%     t = linspace(0, 2*pi, 1000);
%     x = a*cos(t) + 1;  % 椭圆中心点均位移2单位
%     y = b*sin(t);
% 
%     plot(x, y, 'b')  % 'b'表示黑色
% end
% 
% for i = 0:9
%     a = 2.99 + i*0.07;  % a从2开始，增长0.2
%     b = a/3;  % 保持a/b的比例为2:1
%     r = b;
%     rectangle('Position',[-r, -r, 2*r, 2*r], 'Curvature',[1 1], 'EdgeColor', 'r')  % 圆心在(0,0)
%     % 创建椭圆的参数方程
%     t = linspace(0, 2*pi, 1000);
%     x = a*cos(t) + 1;  % 椭圆中心点均位移2单位
%     y = b*sin(t);
% 
%     plot(x, y, 'b')  % 'b'表示黑色
% end

xlim([-3,5]);
ylim([-3,3]);
hold off
