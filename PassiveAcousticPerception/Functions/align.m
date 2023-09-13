function [Signal, BGSignal] = align(SignalParam,Signal)

% raww = Signal;
% if size(raww, 2) > 1
%     raww = mean(raww, 2); % Average the two channels.
% end
% 
% stftMatrix2 = stft(raww, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);
% stftMatrix2 = mag2db(abs(stftMatrix2));
% t1 = (0:size(stftMatrix2,2)-1)*(1024-512)/48000; % 1024是窗口大小，512是重叠长度
% f1 = (-(size(stftMatrix2,1)-1)/2:size(stftMatrix2,1)/2)*48000/size(stftMatrix2,1);
% % Convert the amplitude spectrogram to dB scale
% figure;
% imagesc(t1,f1,stftMatrix2)
% set(gca,'YDir','normal')
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% title('Spectrogram of audio signal after first alignment')
% colorbar;
% 
% xlim([2,4]);
% hold on;
% 
% energy_threshold = 180;  % 能量阈值
% prev_energy = -inf;  % 上一个能量值
% peak_freq = 100000;  % 频率值初始化为 NaN
% peak_time = 0;  % 时间值初始化为 NaN
% found = false;  % 是否找到了峰值
% 
% % 遍历时间
% for t = 1:15
%     % 获取当前时间点的能量值
%     energy_values = stftMatrix2(:, t);
%     %fprintf('Frequency blcok is from %d to %d\n', f1(round(size(stftMatrix2, 1)*2 / 3)), f1(size(stftMatrix2, 1)));
%     % 遍历当前时间点上半部分的能量值
%     for freq_idx = round(size(stftMatrix2, 1)*2/3):round(size(stftMatrix2, 1)-5)
%     %for freq_idx = 1:size(stftMatrix2, 1)/3
%         % 检查频率是否大于 200
%         %fprintf('id:%d, energy_values(id):%d\n', freq_idx, energy_values(freq_idx));
%         if energy_values(freq_idx) > 190
%             % 计算当前点及其上面 4 行点的能量值和的平均值
%             energy_sum = sum(energy_values(freq_idx-4:freq_idx));
%             energy_avg = energy_sum / 5;
%             %fprintf("energy_avg=%f, energy_threshold=%f, f1(freq_idx)=%f, peak_freq=%f\n ", energy_avg, energy_threshold, f1(freq_idx), peak_freq);
%             % 检查能量和的平均值是否满足阈值要求，并且频率是否大于之前找到的最大频率
%             if energy_avg >= energy_threshold && f1(freq_idx) < peak_freq
%                 %fprintf('im in!!!!   ');
%                 %fprintf("energy_avg=%f, energy_threshold=%f, f1(freq_idx)=%f, peak_freq=%f\n ", energy_avg, energy_threshold, f1(freq_idx), peak_freq);
%                 peak_time = t1(t);
%                 tmp = t;% 峰值对应的时间
%                 peak_freq = f1(freq_idx);  % 峰值的频率
%                 peak_energy = energy_values(freq_idx);  % 峰值的能量
% 
%                 found = true;  % 设置找到标志
%             end
%         end
%     end
% end
% 
% % 如果找到峰值，则将峰值标记在图中
% if found
%     scatter(peak_time, peak_freq, 'r', 'filled');
%     fprintf('The peak time is %f, ', peak_time);
%     fprintf('The peak frequency is %f, ', peak_freq);
%     fprintf('The peak energy is %f\n', peak_energy);
%     fprintf('tmp%d\n', tmp);
% else
%     fprintf('No peak was found\n');
% end
% 
% hold off;


% 找到peak_time在t1中的索引
% [~, peak_time_idx] = min(abs(t1-peak_time));
% fprintf('\n\npeak_time=%f, peak_time_idx=%f\n', peak_time, peak_time_idx);
% % 删除stftMatrix2中前peak_time_idx列的数据
% stftMatrix2(:, 1:(peak_time_idx)) = [];
% 
% % 计算peak_time秒对应的样本数
% samples_to_remove = round(SignalParam.SampleFrequency * peak_time);
% samples_to_remove = 0;
% % 删除Signal前peak_time时间的数据
% fprintf('-- Read size: %d, %d;\n',size(Signal,1), size(Signal,2));
% Signal = Signal(samples_to_remove + 1:end,:);
% fprintf('-- Read size: %d, %d;\n',size(Signal,1), size(Signal,2));
% Signall = Signal;
% if size(Signall, 2) > 1
%     Signall = mean(Signall, 2); % Average the two channels.
% end
% stftMatrix3 = stft(Signall, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);
% 
% stftMatrix3 = mag2db(abs(stftMatrix3));
% t1 = (0:size(stftMatrix3,2)-1)*(1024-512)/48000; % 1024是窗口大小，512是重叠长度
% f1 = (-(size(stftMatrix3,1)-1):size(stftMatrix3,1)-1)*48000/size(stftMatrix3,1);
% % Convert the amplitude spectrogram to dB scale
% figure;
% imagesc(t1,f1,stftMatrix3)
% set(gca,'YDir','normal')
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% title('Spectrogram of audio signal after second alignment')
% colorbar;
% xlim([2,4]);

fprintf('-- Read size: %d, %d;\n',size(Signal,1), size(Signal,2));

 
%Background Alignment
BGSignal=zeros(SignalParam.ChirpSize,6);
%fprintf('BGSignal size=%d\n', size(BGSignal));
for ii=1:5
        %fprintf('+Signal:%f\n', Signal(1:SignalParam.ChirpSize,:));
        BGSignal=BGSignal+Signal(1:SignalParam.ChirpSize,:);
        %fprintf('BGSignal=%f\n', BGSignal);
        Signal(1:SignalParam.ChirpSize,:)=[];
end % delete 25 chirps i.e. 2s (48000 samples)
    
    % in total, will delete 2.4s ~ 2.44s, i.e., 115200 ~ 117120 samples
    %fprintf('-- raw size new new: %d;\n',size(raw_data(:,mic_index)));
    BGSignal=BGSignal/5;
end