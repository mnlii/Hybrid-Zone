clear;

% 读取第一个音频文件
[audioData1, sampleRate1] = audioread('./chirp4.wav');
if size(audioData1, 2) > 1
    audioData1 = mean(audioData1, 2); % Average the two channels.
end

% 计算第一个音频的STFT
stftMatrix1 = stft(audioData1, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);
stftMatrix1 = mag2db(abs(stftMatrix1));

% 计算第一个音频的时间和频率矢量
t1 = (0:size(stftMatrix1,2)-1)*(1024-512)/sampleRate1; % 1024是窗口大小，512是重叠长度
f1 = (0:size(stftMatrix1,1)-1)*sampleRate1/size(stftMatrix1,1);

% 读取第二个音频文件
[audioData2, sampleRate2] = audioread('./source/1.wav');
if size(audioData2, 2) > 1
    audioData2 = mean(audioData2, 2); % Average the two channels.
end

% 计算第二个音频的STFT
stftMatrix2 = stft(audioData2, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);
stftMatrix2 = mag2db(abs(stftMatrix2));

% 计算第二个音频的时间和频率矢量
t2 = (0:size(stftMatrix2,2)-1)*size(stftMatrix2,2)/sampleRate2;
f2 = (0:size(stftMatrix2,1)-1)*sampleRate2/size(stftMatrix2,1);

% 创建一个新的图形窗口
figure

% 画出第一个音频的STFT矩阵
subplot(2,1,1);
imagesc(t1,f1,stftMatrix1)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of audio signal 1')
colorbar

% 画出第二个音频的STFT矩阵
subplot(2,1,2);
imagesc(t2,f2,stftMatrix2)
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of audio signal 2')
colorbar

% 对齐两个子图的时间轴
linkaxes([subplot(2,1,1), subplot(2,1,2)], 'x');
xlim([1 4]);  % 设置时间轴范围为0-3s
