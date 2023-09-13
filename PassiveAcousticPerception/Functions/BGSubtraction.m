function [BGSubtractionData,Signal] = BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,fig)
% 
    data = Signal(1:SignalParam.ChirpSize,:);
    Signal(1:SignalParam.ChirpSize,:)=[];
    BGSubtractionData = data - BGSignal * 1;
    %BGSubtractionData = data;

    if fig
        figure;
        for ii=1:6
            subplot(6,3,3*ii-2)
            plot(data(:,ii))
            subplot(6,3,3*ii-1)
            plot(BGSignal(:,ii))
            subplot(6,3,3*ii)
            plot(BGSubtractionData(:,ii))
            ylim([-5e8, 5e8])
        end
        sgtitle('Background Denoise')


        chirp=GenerateFMCWSignal(SignalParam);
        datafiltered = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp.*data);
        s=datafiltered;
        fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
        fft_ = fftshift(fft_,1);
        data_fft_abs=abs(fft_);
        BGSubtractionDatafiltered = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp.*BGSubtractionData);
        s=BGSubtractionDatafiltered;
        fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
        fft_ = fftshift(fft_,1);
        BGSubtractionData_fft_abs=abs(fft_);

        maxplot=max([max(max(data_fft_abs)),max(max(BGSubtractionData_fft_abs))]);
        data_fft_abs=data_fft_abs/maxplot;
        BGSubtractionData_fft_abs=BGSubtractionData_fft_abs/maxplot;

        figure
        for ii=1:6
            subplot(6,2,ii*2-1)
            f_=-1000:1000;
            range = f_+SignalParam.FFTLength/2+1;
            plot(f_,data_fft_abs(range,ii));
            ylim([-1,1]);

            subplot(6,2,ii*2)
            f_=-1000:1000;
            range = f_+SignalParam.FFTLength/2+1;
            plot(f_,BGSubtractionData_fft_abs(range,ii));
            ylim([-1,1]);
        end
    end
end

% data=Signal(1:SignalParam.ChirpSize,:);
% Signal(1:SignalParam.ChirpSize,:)=[];
% 
% % 把每个周期信号分为上下两部分
% midpoint = round(size(data, 1) / 2);
% data_up = data(1:midpoint, :);
% data_down = data(midpoint+1:end, :);
% BGSignal_up = BGSignal(1:midpoint, :);
% BGSignal_down = BGSignal(midpoint+1:end, :);
% 
% % 处理上半部分的数据
% BGSubtractionData_up = data_up - BGSignal_up;
% 
% % 处理下半部分的数据
% BGSubtractionData_down = data_down - BGSignal_down;
% BGSubtractionData = [BGSubtractionData_up; BGSubtractionData_down];
% 
% 
% % 如果需要绘图
% if fig
%     figure;
%     for ii=1:6
%         subplot(6,3,3*ii-2)
%         plot(data_up(:,ii))
%         subplot(6,3,3*ii-1)
%         plot(BGSignal_up(:,ii))
%         subplot(6,3,3*ii)
%         plot(BGSubtractionData_up(:,ii))
%     end
%     sgtitle('背景减法 - 上半部分')
% 
%     figure;
%     for ii=1:6
%         subplot(6,3,3*ii-2)
%         plot(data_down(:,ii))
%         subplot(6,3,3*ii-1)
%         plot(BGSignal_down(:,ii))
%         subplot(6,3,3*ii)
%         plot(BGSubtractionData_down(:,ii))
%     end
%     sgtitle('背景减法 - 下半部分')
% 
%     % Generate two chirp signals to match the sizes of data_up and data_down
%     chirp_up = GenerateFMCWSignal(SignalParam);
%     chirp_down = GenerateFMCWSignal(SignalParam);
% 
%     % Adjust size of chirp signals to match data size
%     chirp_up = chirp_up(1:size(data_up, 1), :);
%     chirp_down = chirp_down(size(data_down, 1)+1:size(data_up, 1)+size(data_down, 1), :);
% 
%     % 上半部分数据的处理
%     datafiltered_up = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_up.*data_up);
%     s = datafiltered_up;
%     fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
%     fft_ = fftshift(fft_,1);
%     data_fft_abs_up = abs(fft_);
%     BGSubtractionDatafiltered_up = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_up.*BGSubtractionData_up);
%     s = BGSubtractionDatafiltered_up;
%     fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
%     fft_ = fftshift(fft_,1);
%     BGSubtractionData_fft_abs_up = abs(fft_);
% 
%     % 下半部分数据的处理
%     datafiltered_down = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_down.*data_down);
%     s = datafiltered_down;
%     fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
%     fft_ = fftshift(fft_,1);
%     data_fft_abs_down = abs(fft_);
%     BGSubtractionDatafiltered_down = filter(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,chirp_down.*BGSubtractionData_down);
%     s = BGSubtractionDatafiltered_down;
%     fft_ = fft(s,SignalParam.FFTLength)/size(s,1); 
%     fft_ = fftshift(fft_,1);
%     BGSubtractionData_fft_abs_down = abs(fft_);
% 
%     maxplot_up = max([max(max(data_fft_abs_up)),max(max(BGSubtractionData_fft_abs_up))]);
%     data_fft_abs_up = data_fft_abs_up/maxplot_up;
%     BGSubtractionData_fft_abs_up = BGSubtractionData_fft_abs_up/maxplot_up;
%     maxplot_down = max([max(max(data_fft_abs_down)),max(max(BGSubtractionData_fft_abs_down))]);
%     data_fft_abs_down = data_fft_abs_down/maxplot_down;
%     BGSubtractionData_fft_abs_down = BGSubtractionData_fft_abs_down/maxplot_down;
% 
%     figure
%     for ii=1:6
%         subplot(12,2,ii*2-1)
%         f_=-1000:1000;
%         range = f_+SignalParam.FFTLength/2+1;
%         plot(f_,data_fft_abs_up(range,ii));
%         ylim([-1,1]);
%         subplot(12,2,ii*2)
%         f_=-1000:1000;
%         range = f_+SignalParam.FFTLength/2+1;
%         plot(f_,BGSubtractionData_fft_abs_up(range,ii));
%         ylim([-1,1]);
%         subplot(12,2,ii*2+11)
%         f_=-1000:1000;
%         range = f_+SignalParam.FFTLength/2+1;
%         plot(f_,data_fft_abs_down(range,ii));
%         ylim([-1,1]);
%         subplot(12,2,ii*2+12)
%         f_=-1000:1000;
%         range = f_+SignalParam.FFTLength/2+1;
%         plot(f_,BGSubtractionData_fft_abs_down(range,ii));
%         ylim([-1,1]);
%     end
% end