clear all; close all;
acou = load("acouResperson2a.mat");
acou = acou.savedData;
acouDistance = acou(:, 1);
acouAngle = acou(:, 2);
acouVelocity = acou(:, 3);
wifi = load("wifiResperson_2a.mat");
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

figure
subplot(2,1,1)
plot(acouVelocityInterp)
title('v_s after interp')
subplot(2,1,2)
plot(wifiVelocityInterp)
title('v_w after interp')

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
figure;
subplot(2,1,1);
plot(timeNew, v_result);
xlabel('Time');
ylabel('Velocity');
legend('vx', 'vy');
title('Velocity versus Time');


% v_t_res = v_result(:,1) + 1i * v_result(:,2);
% v_t_res = v_t_res - mean(v_t_res);
% v_f_res1 = abs(fft(v_t_res));
% v_f_res2 = abs(fft(v_result(:,2)));

tmp = abs(v_result(:,1) + 1i * v_result(:,2));
subplot(2,1,2);
plot(timeNew, tmp);
xlabel('Time');
ylabel('Velocity');

v = v_result(:,1) + 1i * v_result(:, 2);
d1 = cumsum(v) * 0.01;

figure;
plot(real(d1), imag(d1), 'r-', 'LineWidth', 1.5);  % Change from scatter plot to red solid line plot
% hold on;  % Hold current plot so we can add more elements
% %plot(0, 0, '+', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');  % Plot purple plus at (0, 0)
% x = [0, 1.6];
% y = [0, 0];
% plot(x, y, '+', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'LineWidth', 1.5, 'MarkerSize', 30);  % Plot purple plus at (0, 1.6)
% plot(0, 0, 'x', 'MarkerEdgeColor',  [126 47 142]/255, 'LineWidth', 1.5, 'MarkerSize', 30);  % Plot yellow cross at (0, 0)
% legend('Line', 'Wi-Fi', 'Microphone array');  % Add legend
axis([-.5,0.5,-.5,0.5]);  % Set axis limits
