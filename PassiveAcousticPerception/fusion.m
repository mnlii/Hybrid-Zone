clear all;
%timeWifi;
%timeAcou;
load("tr1.mat");
acouDistance = toSave;
load("tr2.mat");
acouAngle = toSave;
load("log7.mat")
wifiDistance = toSave;
acouDistance(end-25:end,:)=[];
acouAngle(end-25:end,:)=[];
wifiDistance(end-20:end,:)=[];
subplot(2,1,1);
plot(wifiDistance);
subplot(2,1,2);
plot(acouDistance);
    



% 确定插值后的数据点数
num_points = length(acouDistance);
num_interp_points = length(wifiDistance);

% 生成插值的时间点
x = linspace(1, num_interp_points, num_points);

% 插值处理
wifiDistance_interp = interp1(wifiDistance, x);

% 绘制插值后的图形
%figure;
%subplot(2, 1, 1);
%plot(wifiDistance_interp);
%title('Interpolated WiFi Distance');
%subplot(2, 1, 2);
%plot(acouDistance);
%title('Acoustic Distance and Angle');

% 更新wifiDistance为插值后的结果
wifiDistance = wifiDistance_interp;
delta_t = 0.04;  % ∆t为0.04秒
acouDistance_shifted = [acouDistance(2:end); 0];  % 将acouDistance数组向前移动一个位置
acouVelocity = (acouDistance_shifted - acouDistance) / delta_t;

% 长度
N = length(wifiDistance);

% 转换角度为弧度
acouAngle = acouAngle * pi / 180; 

% 初始化变量
P_sound = zeros(1, N); %声音源位置
P_wifi = ones(1, N);   %Wifi接收器位置
P_u = zeros(1, N);     %目标物体位置
V_u = zeros(1, N);     %目标物体速度

t_interval = 0.04;
% 计算初始位置
P_u(1) = acouDistance(1) * exp(1i * acouAngle(1));

for i = 2:N
    % 计算速度
    V_u(i) = acouVelocity(i) * exp(1i * acouAngle(i));
    
    % 用速度更新位置
    P_u(i) = P_u(i-1) + V_u(i) * t_interval;
end

% 计算距离
D_us = abs(P_u - P_sound);
D_uw = abs(P_u - P_wifi);

% 计算速度的分量
m1 = real((P_u - P_sound)./abs(P_u - P_sound));
m2 = imag((P_u - P_sound)./abs(P_u - P_sound));
m3 = real(((P_u - P_sound)./abs(P_u - P_sound) + (P_u - P_wifi)./abs(P_u - P_wifi)));
m4 = imag(((P_u - P_sound)./abs(P_u - P_sound) + (P_u - P_wifi)./abs(P_u - P_wifi)));

V_us = - real(V_u .* (P_u - P_sound)./abs(P_u - P_sound));
V_uw = - real(V_u .* ((P_u - P_sound)./abs(P_u - P_sound) + (P_u - P_wifi)./abs(P_u - P_wifi)));

% 融合wifi和声学信息计算速度和距离
vx = (m4 .* V_us - m2 .* V_uw)./(m2 .* m3 - m1 .* m4);
vy = (m3 .* V_us - m1 .* V_uw)./(m2 .* m3 - m1 .* m4);
P_infer = vx + 1i * vy;

% 绘制速度随时间的变化图
figure;
plot(1:N, abs(V_u));
title('Speed over Time');
xlabel('Time');
ylabel('Speed');

% 绘制距离随时间的变化图
figure;
plot(1:N, abs(P_u));
title('Distance over Time');
xlabel('Time');
ylabel('Distance');
