clear;
% %[audioData, sampleRate] = audioread('./chirp4.wav');
% 
% [audioData, sampleRate] = audioread('./source/1.wav');
% if size(audioData, 2) > 1
%     audioData = mean(audioData, 2); % Average the two channels.
% end
% 
% % 计算音频的STFT
% stftMatrix = stft(audioData, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);
% 
% % 取绝对值并转换为分贝
% stftMatrix = mag2db(abs(stftMatrix));
% 
% % 计算时间和频率矢量
% t = (0:size(stftMatrix,2)-1)*size(stftMatrix,2)/sampleRate;
% f = (0:size(stftMatrix,1)-1)*sampleRate/size(stftMatrix,1);
% 
% % 创建一个新的图形窗口
% figure
% 
% % 画出STFT矩阵
% imagesc(t,f,stftMatrix)
% 
% % 翻转y轴方向使得频率从低到高
% set(gca,'YDir','normal')
% 
% % 添加坐标轴标题和图形标题
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% title('Spectrogram of audio signal')
% colorbar

% [audioData, sampleRate] = audioread('./source/1.wav');
[audioData, sampleRate] = audioread('./chirp.wav');
% If the audio data has two columns, it's stereo sound. We should convert it to mono sound.
if size(audioData, 2) > 1
    audioData = mean(audioData, 2); % Average the two channels.
end

% Then perform the spectrogram calculation and display
[S, F, T] = stft(audioData, sampleRate, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);
S = mag2db(abs(S));
F = F * sampleRate;

figure
imagesc(T, F, S);
set(gca,'YDir','normal');

xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of audio signal')
colorbar
