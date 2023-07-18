tic
clc;clear;
close all;
addpath('.\Functions')
FigureParam = InitFigParam();
SignalParam = InitSignalParam();
FilterParam = InitFilterParam(SignalParam);
MicroParam = InitMicroParam(SignalParam,FigureParam);
chirp = GenerateFMCWSignal(SignalParam);
[BPF_b,BPF_a] = DesignBPF(FilterParam.F_Center,FilterParam.BPF_Filter_Pass,FilterParam.BPF_Filter_Stop,SignalParam.SampleFrequency,FigureParam);
% High pass filter (Derive FMCW waves)
[HPF_b,HPF_a] = DesignHPF(FilterParam.HPF_Filter_Pass,FilterParam.HPF_Filter_Stop,SignalParam.SampleFrequency,FigureParam);
% Low pass filter (Derive FMCW waves)
[LPF_b,LPF_a] = DesignLPF(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,SignalParam.SampleFrequency,FigureParam);
[Signal, mic_index] = Load6MicData(SignalParam,FilterParam,'/3');
[Signal, BGSignal] = align(SignalParam,Signal);
%[data,Signal]=BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,true);
[data,Signal]=BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,true);
shift = 0;
[P_up, P_down, v, d, shift,tmp] = Estimation(SignalParam,FilterParam,MicroParam,data,true, shift);
x = tmp;
trajectory_up = P_up;
trajectory_down = P_down;

trajectory = [P_up; abs(P_down)];
velocity = v;
distance = d;
% figure
% plot(P(1),P(2),'o');
% axis([0.3 0.9 -60 60])
while size(Signal,1)>SignalParam.ChirpSize
    [data, Signal]=BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,false);
    [P_up, P_down, v, d, shift, tmp] = Estimation(SignalParam,FilterParam,MicroParam,data,false, shift);
    x = [x; tmp];
    %fprintf('%d    %d\n',P_up(:,2), P_down(:,2));
    trajectory_up=[trajectory_up;P_up];
    trajectory_down=[trajectory_down;P_down];
    trajectory = [trajectory; P_up; abs(P_down)];
    velocity = [velocity; v];
    distance = [distance; d];
%     hold on;
%     plot(P(1),P(2),'o');
%     axis([0.3 0.9 -60 60]) 

%     pause(SignalParam.ChirpT)

end





load('FullAngleMat.mat');  % 这应该给你一个变量，比如叫做FullAngleMat
num_columns = size(x, 2);
%
segment_length = round(size(x, 1) / 10) - 1;
x = x((1:segment_length * 10), :);
% 将x分成10段
x_segments = mat2cell(x, repmat(segment_length, 1, 10), num_columns);

% 循环处理每段数据
for seg_idx = 1:10
    % 获取当前段的数据
    x_segment = x_segments{seg_idx};
    
    [c, lags] = xcorr(x_segment, 'coeff');
    maxLag = zeros(num_columns, num_columns); % 存储每个子图中最大频率对应的lag
    
    %figure; % 每段数据创建一个新的图形
    for i = 1:num_columns
        for j = 1:num_columns
            %subplot(num_columns, num_columns, (i-1)*num_columns+j);
            idx = find(abs(lags)<=8);
            %stem(lags(idx),c(idx,(i-1)*num_columns+j));
            % title(sprintf('Cross-correlation between column %d and %d', i, j));
            % xlim([-8 8]);  
            [~, maxIdx] = max(abs(c(idx,(i-1)*num_columns+j)));
            maxLag(i,j) = lags(idx(maxIdx));
        end
    end
    %disp(maxLag);
    
    distances = zeros(360, 1);
    for angle = 1:360
        distances(angle) = norm(FullAngleMat(:, :, angle) - maxLag, 'fro');
    end
    
    min_distance = min(distances);
    estimated_angles = find(distances == min_distance);
    fprintf('%d ', estimated_angles);
    %display(estimated_angles);
end


delta = 6.8 * shift / 1920;

trajectory(:,1) = trajectory(:,1) + delta;


position = zeros(size(trajectory, 1), 2); % Initialize position array

% Convert angle from degrees to radians
theta_rad = deg2rad(trajectory(:,2));

% Compute the coordinates
position(:,1) = trajectory(:,1) .* cos(theta_rad);
position(:,2) = trajectory(:,1) .* sin(theta_rad);


figure
subplot(4,1,1)
position = hampel(position, 10, 1);
position(:,1) = smooth(position(:,1), 10);
position = hampel(position, 5, 0.5);
plot(position(:,1), position(:,2), '.'); % Plot y against x using dots
xlim([0,1]);
ylim([0,1]);

title('Plot of y versus x');
xlabel('x');
ylabel('y');
title('merge-chirp position')

trajectory = hampel(trajectory, 10, 1);
trajectory(:,1) = smooth(trajectory(:,1), 10);
trajectory = hampel(trajectory, 5, 0.5);

subplot(4,1,2)
plot(trajectory(:,1))
ylim([0 1])

title('merge-chirp distance')

subplot(4,1,3)
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

% 绘图
plot(velocity_smoothed(:,1))
%plot(velocity(:,1))


title('fv -> v')
subplot(4,1,4)
distance = hampel(distance, 10, 1);
distance = smooth(distance, 10);
distance = hampel(distance, 5, 0.5);
plot(distance(:,1))

title('fd -> d')
% % FileName = 'acouRes3.mat';
% save(FileName,'trajectory');