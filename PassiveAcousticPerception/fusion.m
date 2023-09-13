function [] = fusion(filname)
acou = load(['acouRes',filname,'.mat']);
acou = acou.savedData;
acouDistance = acou(:, 1);
acouAngle = acou(:, 2);
acouVelocity = acou(:, 3);
wifi = load(['wifiRes', filname, '.mat']);
wifi = wifi.savedData;
wifiVelocity = wifi(:,1);
timeWifi = wifi(:,2)';

time = min(length(acouVelocity) * 0.08, timeWifi(:,end));

% Define the original time vectors
timeAcou = linspace(0, time, length(acouVelocity));

% Define the new time vector
timeNew = linspace(timeWifi(:,1), time, (time - timeWifi(:,1)) * 101);

% Perform the interpolation
wifiVelocityInterp = interp1(timeWifi, wifiVelocity, timeNew);
acouDistanceInterp = interp1(timeAcou, acouDistance, timeNew);
acouAngleInterp = interp1(timeAcou, acouAngle, timeNew);
acouVelocityInterp = interp1(timeAcou, acouVelocity, timeNew);
acouPosition(:,1) = acouDistanceInterp .* cos(acouAngleInterp);
acouPosition(:,2) = acouDistanceInterp .* sin(acouAngleInterp);
% figure
% subplot(2,1,1)
% plot(acouVelocityInterp)
% title('v_s after interp')
% subplot(2,1,2)
% plot(wifiVelocityInterp)
% title('v_w after interp')

% Initialize velocity result vector
v_result = zeros(length(timeNew), 2);

d = acouDistanceInterp(1);
theta = acouAngleInterp(1);

for idx = 1:length(timeNew)
    % Define known variables.
    xs = 0;
    ys = 0;
    xt = 0;
    yt = 0;
    xr = 0;
    yr = 1.6;
    vs = acouVelocityInterp(idx);
    vw = wifiVelocityInterp(idx);

    % Calculate xh and yh

    xh = d * cos(theta);
    yh = d * sin(theta);
    % Calculate lh
    lh = sqrt(xh^2 + yh^2);
    mt = sqrt((xh - xt)^2 + (yh - yt)^2);
    mr = sqrt((xh - xr)^2 + (yh - yr)^2);

    % Calculate alpha_x and alpha_y
    alpha_x = (xh - xs) / lh;
    alpha_y = (yh - ys) / lh;

    % Calculate beta_x and beta_y
    beta_x = ((xh - xt) / mt) + ((xh - xr) / mr);
    beta_y = ((yh - yt) / mt) + ((yh - yr) / mr);

    % Define L matrix and v_obs vector
    L = [alpha_x, alpha_y; beta_x, beta_y];
    v_obs = [vs; vw];

    % Calculate v
    v = inv(L' * L) * L' * v_obs;
    
    % Store the result
    v_result(idx, :) = v;

    % xh = xh + v_result(idx, 1) * 0.005;
    % yh = yh + v_result(idx, 2) * 0.005;
end

% Plot the velocity
% figure;
% subplot(2,1,1);
% plot(timeNew, v_result);
% xlabel('Time');
% ylabel('Velocity');
% legend('vx', 'vy');
% title('Velocity versus Time');


% v_t_res = v_result(:,1) + 1i * v_result(:,2);
% v_t_res = v_t_res - mean(v_t_res);
% v_f_res1 = abs(fft(v_t_res));
% v_f_res2 = abs(fft(v_result(:,2)));

tmp = abs(v_result(:,1) + 1i * v_result(:,2));
% subplot(2,1,2);
% plot(timeNew, tmp);
% xlabel('Time');
% ylabel('Velocity');

v = v_result(:,1) + 1i * v_result(:, 2);
d1 = cumsum(v) * 0.01;
fusedPosition(:, 1) = real(d1);
fusedPosition(:, 2) = imag(d1);
 figure;
 plot(real(d1), imag(d1), 'r-', 'LineWidth', 1.5);  % Change from scatter plot to red solid line plot
hold on;  % Hold current plot so we can add more elements
%plot(0, 0, '+', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');  % Plot purple plus at (0, 0)
x = [0, 1.6];
y = [0, 0];
plot(x, y, '+', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'LineWidth', 1.5, 'MarkerSize', 30);  % Plot purple plus at (0, 1.6)
plot(0, 0, 'x', 'MarkerEdgeColor',  [126 47 142]/255, 'LineWidth', 1.5, 'MarkerSize', 30);  % Plot yellow cross at (0, 0)
legend('Line', 'Wi-Fi', 'Microphone array');  % Add legend
axis([-.5,0.5,-.5,0.5]);  % Set axis limits

merged_velocity = zeros(800,2);

length_wifi = length(wifiVelocity);
length_acou = length(acouVelocity);

% 将向量复制到新的矩阵中，只复制到向量的长度处
merged_velocity(1:length_wifi,1) = wifiVelocity;
merged_velocity(1:length_acou,2) = acouVelocity;
save(['./gesture/', filname,'.mat'], 'merged_velocity');



% Known time points
known_times = [0, time/3, time*2/3, time];
	
	
	


% Known x and y coordinates, replace x0, y0, ..., x3, y3 with real values
known_xs = [0, -.7, -0.05, -0.65];
known_ys = [0, -.35, -.53, -0.7];

% Times at which to interpolate

% Perform the interpolation
inter_xs = interp1(known_times, known_xs, timeNew);
inter_ys = interp1(known_times, known_ys, timeNew);

% Combine the interpolated x and y coordinates
inter_coords = [inter_xs; inter_ys]';

error_fusion_no = sqrt(sum((inter_coords - fusedPosition).^2, 2));
% error_acou = sqrt(sum((inter_coords - acouPosition).^2, 2));
% error_acou = error_acou - .4;

figure;
boxplot([error_fusion_no, error_fusion_no, error_fusion_no, error_fusion_no, error_fusion_no, error_fusion_no, error_fusion_no]);
set(gca,'LineWidth',1.2, 'FontSize', 20);
set(gcf,'Position',[200,200,600,400])
grid on;


% figure;
% h1 = cdfplot(error_fusion_no);
% set(h1,'LineWidth',2);
% hold on;
% h1 = cdfplot(error_fusion_no);
% set(h1,'LineWidth',2);
% hold on;
% h1 = cdfplot(error_fusion_no);
% set(h1,'LineWidth',2);
% legend('Sampling rate = 24000 Hz', 'Sampling rate2', 'Sampling rate3', 'Location', 'southeast');
% title('pure acoustic cdf', 'FontSize', 20)
% xlim([0,2]);
% xlabel('Distance (m)')
% set(gca,'LineWidth',1.2, 'FontSize', 20);
% set(gcf,'Position',[200,200,600,400])
% grid on;


end


