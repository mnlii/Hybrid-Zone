tic
clc;clear;
close all;
addpath('.\Functions')
name = 'person2a';
FigureParam = InitFigParam();
SignalParam = InitSignalParam();
FilterParam = InitFilterParam(SignalParam);
MicroParam = InitMicroParam(SignalParam,FigureParam);
chirp = GenerateFMCWSignal(SignalParam);
[BPF_b,BPF_a] = DesignBPF(FilterParam.F_Center,FilterParam.BPF_Filter_Pass,FilterParam.BPF_Filter_Stop,SignalParam.SampleFrequency);
[HPF_b,HPF_a] = DesignHPF(FilterParam.HPF_Filter_Pass,FilterParam.HPF_Filter_Stop,SignalParam.SampleFrequency);
[LPF_b,LPF_a] = DesignLPF(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,SignalParam.SampleFrequency);
[Signal, mic_index] = Load6MicData(SignalParam,FilterParam,name);
[Signal, BGSignal] = align(SignalParam,Signal);
[data,Signal]=BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,false);

shift = 0;
[P_up, P_down, v, d, shift,tmp] = Estimation(SignalParam,FilterParam,MicroParam,data,true, shift);
x = tmp;
trajectory_up = P_up;
trajectory_down = P_down;

trajectory = [P_up; abs(P_down)];
velocity = v;
distance = d;
%shift = 0;
while size(Signal,1)>SignalParam.ChirpSize
    [data, Signal]=BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,false);
    [P_up, P_down, v, d, shift, tmp] = Estimation(SignalParam,FilterParam,MicroParam,data, false, shift);
    x = [x; tmp];
    trajectory_up=[trajectory_up;P_up];
    trajectory_down=[trajectory_down;P_down];
    trajectory = [trajectory; P_up; abs(P_down)];
    velocity = [velocity; v];
    distance = [distance; d];
end
delta = 6.8 * shift / 1920;
trajectory(:,1) = trajectory(:,1) + delta;

%position = zeros(size(trajectory, 1), 2); % Initialize position array

% Convert angle from degrees to radians
theta_rad = deg2rad(trajectory(:,2));
theta_rad = hampel(theta_rad,20,0.05);
theta_rad = smooth(theta_rad,20);
% Compute the coordinates



trajectory = hampel(trajectory, 10, 1);
trajectory(:,1) = smooth(trajectory(:,1), 10);
trajectory = hampel(trajectory, 5, 0.5);

trajectory_up = hampel(trajectory_up, 10, 1);
trajectory_up(:,1) = smooth(trajectory_up(:,1), 10);
trajectory_up = hampel(trajectory_up, 5, 0.5);

z_scores = abs(zscore(velocity)); % 计算z-scores
TF = z_scores > 2; % 创建一个逻辑向量，标记所有z-scores绝对值大于2的数据点

% 创建一个新的速度向量，其中离群值被替换为NaN
velocity_cleaned = velocity;
velocity_cleaned(TF) = NaN;

% 创建一个时间向量，这将用于插值
t = 1:length(velocity);

% 使用'pchip'插值方法填补离群值
velocity_cleaned = interp1(t(~TF), velocity(~TF), t, 'pchip');

% 进行Hampel滤波和移动平均平滑
velocity_smoothed = hampel(velocity_cleaned, 10, 1);
velocity_smoothed = smooth(velocity_smoothed, 10);
velocity_smoothed = hampel(velocity_smoothed, 5, 0.5);


distance = hampel(distance, 10, 1);
distance = smooth(distance, 10);
distance = hampel(distance, 5, 0.5);


d2v = diff(distance) / 0.08;
d2v = [d2v(1,:);d2v];
d2v_mean = (d2v(:, 1) + velocity_smoothed(:, 1))/2;

d2v_up = diff(trajectory_up(:, 1))/0.08;
d2v_up = [d2v_up(1,:);d2v_up];

d2v_mean = hampel(d2v_mean, 10, 1);
d2v_mean = smooth(d2v_mean, 10);
d2v_mean = hampel(d2v_mean, 5, 0.5);


d2v_up = hampel(d2v_up, 10, 1);
d2v_up = smooth(d2v_up, 10);
d2v_up = hampel(d2v_up, 5, 0.5);

figure; 
subplot(4,1,1);
plot(trajectory_up(:, 1));
subplot(4,1,2);
plot(distance(:, 1));
subplot(4, 1, 3);
plot(d2v_up(:,1))
subplot(4,1,4)
plot(d2v_mean(:,1))





theta_rad_avg = mean(reshape(theta_rad, [2, size(theta_rad,1)/2]), 1)';
savedData = [distance(:,1), theta_rad_avg, d2v_mean(:,1)];
FileName = ['acouRes', name, '.mat'];
save(FileName, 'savedData');

position(:,1) = distance(:,1) .* cos(theta_rad_avg);
position(:,2) = distance(:,1) .* sin(theta_rad_avg);

figure;
position = hampel(position, 10, 0.1);
position(:,1) = smooth(position(:,1), 10);
position(:,2) = smooth(position(:,2), 10);
position = hampel(position, 5, 0.05);
plot(position(:,1), position(:,2), 'r-'); % Plot y against x using dots
xlim([0,1]);
ylim([0,1]);