function [P_up, P_down, v, d, idx_shift, x] = Estimation(SignalParam,FilterParam,MicroParam,data,fig, idx_shift)
%     chirp=GenerateFMCWSignal(SignalParam);
%     data=data./max(data);
%     datafiltered = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp.*data);
% 
%     s=datafiltered;
%     fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
%     fft_ = fftshift(fft_,1);
%     data_fft_abs=abs(fft_);
%     [~,idx]=max(data_fft_abs);
%     r=(idx-SignalParam.SampleFrequency/2-1)/SignalParam.ChirpBW*SignalParam.ChirpT*2*SignalParam.Speed;
%     TF = ~isoutlier(r).*(1:6);
%     TF = TF(TF~=0);
%     meanr=mean(r(TF));
%     diffr=r-meanr;
%     Energy=zeros(121,1);
%     for theta=-60:60
%         Energy(theta+61)=sum(abs(cos((MicroParam.Angles(TF)-theta)/180*pi)*MicroParam.toEdge+diffr(TF)));
%     end
%     [~,idx]=min(Energy);
%     theta=idx-61;
%     P=[meanr/2,theta];
%     if fig
%         f_=-1000:1000;
%         range = f_+SignalParam.FFTLength/2+1;
%         figure
%         for ii=1:6
%             subplot(6,1,ii)
%             plot(f_,data_fft_abs(range,ii));
%         end
%     end
%     midpoint = round(size(chirp, 1) / 2);
%     P_up = P;
%     P_down = 0;
% end



chirp=GenerateFMCWSignal(SignalParam);
midpoint = round(size(chirp, 1) / 2);
chirp_up = chirp(1:midpoint, :);
%chirp_up = circshift(chirp_up, [idx_shift, 0]);
chirp_down = chirp(midpoint+1:end, :);
%chirp_down = circshift(chirp_down, [idx_shift, 0]);
data=data./max(data);
data_up = data(1:midpoint, :);
%data_up = data_up - mean(data_up);
data_down = data(midpoint+1:end, :);

% 处理上部分数据
datafiltered_up = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_up.*data_up);
% raw_data = datafiltered_up;
%     if size(raw_data, 2) > 1
%     raw = mean(raw_data, 2); % Average the two channels.
%     end
%     %[S, F, T] = spectrogram(raw, hamming(SignalParam.SampleFrequency), round(0.5 * SignalParam.SampleFrequency), SignalParam.SampleFrequency, SignalParam.SampleFrequency, 'yaxis');
% stftMatrix1 = stft(raw, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', 1024);
% 
% stftMatrix1 = mag2db(abs(stftMatrix1));
% t1 = (0:size(stftMatrix1,2)-1)*(1024-512)/48000; % 1024是窗口大小，512是重叠长度
% f1 = (-size(stftMatrix1,1)+1:size(stftMatrix1,1)-1)*48000/size(stftMatrix1,1);
% % Convert the amplitude spectrogram to dB scale
% figure;
% imagesc(t1,f1,stftMatrix1)
% set(gca,'YDir','normal')
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% title('Spectrogram of raw_data')
% colorbar;

s_up = datafiltered_up;
s_up = s_up - mean(s_up);
fft_up = fft(s_up, SignalParam.FFTLength) / size(s_up,1); 
fft_up = fftshift(fft_up,1);
data_fft_abs_up = abs(fft_up);
[~,idx_up] = max(data_fft_abs_up);


frequency_up = idx_up - SignalParam.SampleFrequency/2 - 1;
%fprintf('max frequency in up chrip: %f\n',frequency_up);
goodvalue_up = ~isoutlier(frequency_up).*(1:6);
%fprintf('max frequency in up chrip: %f\n',g);
goodvalue_up = goodvalue_up(goodvalue_up~=0);
goodvalue_up = mean(frequency_up(goodvalue_up));
%fprintf('mean in up chrip: %f\n',goodvalue_up);

r_up = (idx_up - SignalParam.SampleFrequency/2 - 1) / SignalParam.ChirpBW * SignalParam.ChirpT * SignalParam.Speed;
TF_up = ~isoutlier(r_up).*(1:6);
TF_up = TF_up(TF_up~=0);
meanr_up = mean(r_up(TF_up));
diffr_up = r_up - meanr_up;
%diffr_up = min(max(diffr_up, -0.0475), 0.0475);

Energy_up = zeros(121,1);
for theta=-60:300
    Energy_up(theta+61) = sum(abs(cos((MicroParam.Angles(TF_up)-theta)/180*pi)*MicroParam.toEdge + diffr_up(TF_up)));
end
[~,idx_up] = min(Energy_up);
theta_up = idx_up - 61;
P_up = [meanr_up/2, theta_up];
% 处理下部分数据
datafiltered_down = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_down.*data_down);
s_down = datafiltered_down;
fft_down = fft(s_down, SignalParam.FFTLength) / size(s_down,1); 
fft_down = fftshift(fft_down,1);
data_fft_abs_down = abs(fft_down);
[~,idx_down] = max(data_fft_abs_down);




frequency_down = idx_down - SignalParam.SampleFrequency/2 - 1;
%fprintf('frequency in down chrip: %f\n',frequency_down);
goodvalue_down = ~isoutlier(frequency_down).*(1:6);
goodvalue_down = goodvalue_down(goodvalue_down~=0);
goodvalue_down = mean(frequency_down(goodvalue_down));
%fprintf('mean in down chrip: %f\n',goodvalue_down);

goodvalue_down = abs(goodvalue_down);

fd = (goodvalue_up + goodvalue_down )/2;
fv = (goodvalue_up - goodvalue_down )/2;

%fprintf('fd = %f, fv = %f\n\n',fd, fv);

v = (fv*SignalParam.Speed)/(20000);
d = fd*SignalParam.Speed*SignalParam.ChirpT/2/SignalParam.ChirpBW;


r_down = (idx_down - SignalParam.SampleFrequency/2 - 1) / SignalParam.ChirpBW * SignalParam.ChirpT * SignalParam.Speed;
TF_down = ~isoutlier(r_down).*(1:6);
TF_down = TF_down(TF_down~=0);
meanr_down = mean(r_down(TF_down));
diffr_down = r_down - meanr_down;
x = [diffr_up; abs(diffr_down)];

Energy_down = zeros(121,1);
for theta=-60:300
    Energy_down(theta+61) = sum(abs(cos((MicroParam.Angles(TF_down)-theta)/180*pi)*MicroParam.toEdge + diffr_down(TF_down)));
end
[~,idx_down] = min(Energy_down);
theta_down = idx_down - 61 - 180;
P_down = [meanr_down/2, theta_down];

% 如果需要绘图
if fig
    f_=-2000:2000;
    range = f_ + SignalParam.FFTLength/2 + 1;

    % 绘制上半部分数据的图
    figure
    for ii=1:6
        subplot(6,1,ii)
        plot(f_, data_fft_abs_up(range, ii));
    end

    % 绘制下半部分数据的图
    figure
    for ii=1:6
        subplot(6,1,ii)
        plot(f_, data_fft_abs_down(range, ii));
    end


% Number of shifts
nShifts = 28; 

% Number of iterations
nIterations = 100; 

shiftt = zeros(1, nIterations+1);
% Loop for shifting
for i = 0:nIterations
    chirp_shift = circshift(chirp_up, [i * nShifts, 0]);
    shift = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_shift.*data_up);
    tmp = shift;
    tmp = tmp - mean(tmp);
    fft_shift = fft(tmp, SignalParam.FFTLength) / size(tmp,1); 
    fft_shift = fftshift(fft_shift,1);
    data_shift = abs(fft_shift);
    [idx_shift,value] = max(data_shift);
    fprintf('%d ', value);
%fprintf('frequency in down chrip: %f\n',frequency_down);
    goodvalue_shift = ~isoutlier(idx_shift).*(1:6);
    goodvalue_shift = goodvalue_shift(goodvalue_shift~=0);
    goodvalue_shift = mean(idx_shift(goodvalue_shift));
    shiftt(i+1) = goodvalue_shift;
end
figure;
plot(0:nIterations, shiftt);
xlabel('Iteration');
ylabel('Good Value Shift');
title('Good Value Shift over Iterations');
[~, idx_shift] = max(shiftt);
idx_shift = idx_shift - 1;
idx_shift = 0;
end




















% 
% chirp=GenerateFMCWSignal(SignalParam);
% chirp = circshift(chirp, [i * nShifts, 0]);
% midpoint = round(size(chirp, 1) / 2);
% chirp_up = chirp(1:midpoint, :);
% chirp_down = chirp(midpoint+1:end, :);
% data=data./max(data);
% data_up = data(1:midpoint, :);
% %data_up = data_up - mean(data_up);
% data_down = data(midpoint+1:end, :);
% 
% % 处理上部分数据
% datafiltered_up = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_up.*data_up);
% s_up = datafiltered_up;
% s_up = s_up - mean(s_up);
% fft_up = fft(s_up, SignalParam.FFTLength) / size(s_up,1); 
% fft_up = fftshift(fft_up,1);
% data_fft_abs_up = abs(fft_up);
% [~,idx_up] = max(data_fft_abs_up);
% 
% 
% frequency_up = idx_up - SignalParam.SampleFrequency/2 - 1;
% %fprintf('max frequency in up chrip: %f\n',frequency_up);
% goodvalue_up = ~isoutlier(frequency_up).*(1:6);
% %fprintf('max frequency in up chrip: %f\n',g);
% goodvalue_up = goodvalue_up(goodvalue_up~=0);
% goodvalue_up = mean(frequency_up(goodvalue_up));
% %fprintf('mean in up chrip: %f\n',goodvalue_up);
% 
% r_up = (idx_up - SignalParam.SampleFrequency/2 - 1) / SignalParam.ChirpBW * SignalParam.ChirpT * SignalParam.Speed;
% TF_up = ~isoutlier(r_up).*(1:6);
% TF_up = TF_up(TF_up~=0);
% meanr_up = mean(r_up(TF_up));
% diffr_up = r_up - meanr_up;
% Energy_up = zeros(121,1);
% for theta=-60:60
%     Energy_up(theta+61) = sum(abs(cos((MicroParam.Angles(TF_up)-theta)/180*pi)*MicroParam.toEdge + diffr_up(TF_up)));
% end
% [~,idx_up] = min(Energy_up);
% fprintf('distance = %f', P_up);
% P_up = P_up + [meanr_up/2, theta_up];
% fprintf(', %f\n', meanr_up/2);
% % 处理下部分数据
% datafiltered_down = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_down.*data_down);
% s_down = datafiltered_down;
% fft_down = fft(s_down, SignalParam.FFTLength) / size(s_down,1); 
% fft_down = fftshift(fft_down,1);
% data_fft_abs_down = abs(fft_down);
% [~,idx_down] = max(data_fft_abs_down);
% P_down = P_down + [meanr_down/2, theta_down];
% 
% 
% 
% 
% frequency_down = idx_down - SignalParam.SampleFrequency/2 - 1;
% %fprintf('frequency in down chrip: %f\n',frequency_down);
% goodvalue_down = ~isoutlier(frequency_down).*(1:6);
% goodvalue_down = goodvalue_down(goodvalue_down~=0);
% goodvalue_down = mean(frequency_down(goodvalue_down));
% %fprintf('mean in down chrip: %f\n',goodvalue_down);
% 
% goodvalue_down = abs(goodvalue_down);
% fd = (goodvalue_up + goodvalue_down )/2;
% d2 = fd*SignalParam.Speed*SignalParam.ChirpT/2/SignalParam.ChirpBW;
% d = d1 + d2;


end



    % factor=2i*pi*(SignalParam.ChirpF0+SignalParam.ChirpBW/SignalParam.ChirpT*(SignalParam.SearchChirpRange'-1)/SignalParam.SampleFrequency);
    % originalSignal = reshape(datafiltered(SignalParam.SearchChirpRange,:),[],1);
    % remainSignal=originalSignal;
    % 
    % r=0.2:0.01:0.6;
    % theta=-60:60;
    % Energy=zeros(length(r),length(theta));
    % for ridx=1:length(r)
    %     for thetaidx=1:length(theta)
    %         tao=(2*r(ridx)-cos((MicroParam.Angles-theta(thetaidx))/180*pi)*MicroParam.toEdge)/SignalParam.Speed;
    %         rebuild=reshape(exp(-factor*tao),[],1);
    %         Energy(ridx,thetaidx)=abs(remainSignal'*rebuild);
    %     end
    % end
    % 
    % maxEnergy=max(max(Energy));
    % [ridx,thetaidx]=find(Energy==maxEnergy);
    % maxr=r(ridx(1));
    % maxtheta=theta(thetaidx(1));
    % maxamplitude=maxEnergy/length(originalSignal);
    % 
    % tao=(2*maxr-cos((MicroParam.Angles-maxtheta)/180*pi)*MicroParam.toEdge)/SignalParam.Speed;
    % rebuild=maxamplitude*reshape(exp(factor*tao),[],1);
    % remainSignal=remainSignal-rebuild;
    % 
    % P=[maxr,maxtheta];
    % 
    % if fig
    %     figure
    %     imagesc(theta,r,Energy);
    %     title(['距离:',num2str(maxr),'   ','角度:',num2str(maxtheta),'   ','振幅:',num2str(round(maxamplitude,4))])
    % end