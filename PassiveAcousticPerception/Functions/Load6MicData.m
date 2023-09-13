function [Signal, mic_index] = Load6MicData(SignalParam,FilterParam,filname)
   
    %% 读取数据

    raw_data = double(audioread(['C:\Users\mengn\Desktop\summer\data\0721\',filname,'.wav'],'native')); % filname='0613';


% 
    if size(raw_data, 2) > 1
    raw = mean(raw_data, 2); % Average the two channels.
    end
    %[S, F, T] = spectrogram(raw, hamming(SignalParam.SampleFrequency), round(0.5 * SignalParam.SampleFrequency), SignalParam.SampleFrequency, SignalParam.SampleFrequency, 'yaxis');
stftMatrix1 = stft(raw, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);

stftMatrix1 = mag2db(abs(stftMatrix1));
t1 = (0:size(stftMatrix1,2)-1)*(1024-512)/48000; % 1024是窗口大小，512是重叠长度
f1 = (-size(stftMatrix1,1)+1:size(stftMatrix1,1)-1)*48000/size(stftMatrix1,1);
% Convert the amplitude spectrogram to dB scale
% figure;
% imagesc(t1,f1,stftMatrix1)
% set(gca,'YDir','normal')
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% title('Spectrogram of raw_data')
% colorbar;

    %raw_data(1:round(SignalParam.SampleFrequency * time_diff),:) = []; %delete 0.07s
    pre = raw_data(1:SignalParam.SampleFrequency*0.5,:);
    fprintf("%d", length(raw_data));
    raw_data(1:SignalParam.SampleFrequency,:) = []; %delete 1s
    stdNum = sum(pre.*pre,1);
    [~, idx] = mink(stdNum,2);
    startIdx=mod(max(idx),8)+1;
    if min(idx)==1 && max(idx)==8
        startIdx=2;
    end
    mic_index=ceil(mod(startIdx:startIdx+5,8.1));

    % figure;
    % for ii=1:8
    %     subplot(8,1,ii);
    %     plot(pre(:,ii));
    % end
    % sgtitle(num2str(mic_index))
    % 
    % Sr=raw_data(1:SignalParam.ChirpSize*3,mic_index); 
    % figure;
    % for ii=1:6
    %     subplot(6,3,3*ii-2);
    %     plot(Sr(:,ii));
    % end

    %% 对齐数据

    %粗略对齐
    chirp=GenerateFMCWSignal(SignalParam);
    base=real(chirp);
    xcorrNum=15;
    Sr=raw_data(1:xcorrNum*SignalParam.ChirpSize,mic_index);
    raw_data(1:xcorrNum*SignalParam.ChirpSize,:)=[]; % delete 10 chirps i.e. 0.8s
    AlignIdx1=zeros(6,1);
    for ii=1:6
        % subplot(6,3,3*ii-1);
        corrIdx=xcorr(Sr(:,ii),base);
        % plot(corrIdx)
        [~,idx]=max(corrIdx);
        AlignIdx1(ii)=mod(idx-size(Sr,1),SignalParam.ChirpSize)+1;
    end
    
    Sr=raw_data((1:SignalParam.ChirpSize)+round(mean(rmoutliers(AlignIdx1))),mic_index); 
    AlignIdx2=zeros(6,1);
    for ii=1:6
        % subplot(6,3,3*ii);
        filtered = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp.*Sr(:,ii));
        s=filtered;
        fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
        fft_ = fftshift(fft_,1);
        fft_abs=abs(fft_);
        f_=-300:300;
        range = f_+SignalParam.FFTLength/2+1;
        % plot(f_,fft_abs(range,:));
        [~,idx]=max(fft_abs);
        AlignIdx2(ii)=(idx-SignalParam.FFTLength/2-1)/SignalParam.ChirpBW*SignalParam.ChirpT*SignalParam.SampleFrequency;
    end
    %sgtitle(['StartMic: ',num2str(startIdx),'  AlignDiff: ',num2str(max(AlignIdx1)-min(AlignIdx1)),'  AlignDiff: ',num2str(max(AlignIdx2)-min(AlignIdx2))])
    
    % already deleted 1.4s, i.e. 67200 samples
    fprintf('-- raw size old: %d;\n',size(raw_data));
    raw_data(1:mod(round(mean(rmoutliers(AlignIdx1))+mean(rmoutliers(AlignIdx2))),SignalParam.ChirpSize),:)=[]; %minor deletion, delete < 1 chirp, i.e. 0.08s
   

    fprintf('-- raw size new: %d;\n',size(raw_data));
   

    fprintf('-- raw size new: %d;\n',size(raw_data(:,mic_index)));
    % 背景减法
    

    Signal=raw_data(:, mic_index);
    fprintf('-- Read size: %d, %d;\n',size(Signal,1), size(Signal,2));


