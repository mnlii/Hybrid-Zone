tic
clc;clear;
close all;
addpath('.\Functions')

%%Section 1: Define parameters and Read raw data
FigureParam = InitFigParam();
SignalParam = InitSignalParam();
FilterParam = InitFilterParam(SignalParam);
MicroParam = InitMicroParam(SignalParam,FigureParam);

%%Section 2: Generate transmitted signals
chirp = GenerateFMCWSignal(SignalParam);

%%Section 3: Design filter to seprate single tone and FMCW signals
%Band pass filter (Derive single tone)
[BPF_b,BPF_a] = DesignBPF(FilterParam.F_Center,FilterParam.BPF_Filter_Pass,FilterParam.BPF_Filter_Stop,SignalParam.SampleFrequency,FigureParam);
% High pass filter (Derive FMCW waves)
[HPF_b,HPF_a] = DesignHPF(FilterParam.HPF_Filter_Pass,FilterParam.HPF_Filter_Stop,SignalParam.SampleFrequency,FigureParam);
% Low pass filter (Derive FMCW waves)
[LPF_b,LPF_a] = DesignLPF(FilterParam.LPF_Filter_Pass,FilterParam.LPF_Filter_Stop,SignalParam.SampleFrequency,FigureParam);

%%Section 4: Load the voice signal / Generate simulation signal
[Signal,BGSignal] = Load6MicData(SignalParam,FilterParam,'2');

%读取数据
[data,Signal]=BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,true);
P = Estimation(SignalParam,FilterParam,MicroParam,data,true);

%% Section 5: 
trajectory=P;
% figure
% plot(P(1),P(2),'o');
% axis([0.3 0.9 -60 60])
while size(Signal,1)>SignalParam.ChirpSize
    [data,Signal]=BGSubtraction(SignalParam,FilterParam,Signal,BGSignal,false);
    P = Estimation(SignalParam,FilterParam,MicroParam,data,false);
    trajectory=[trajectory;P];
%     hold on;
%     plot(P(1),P(2),'o');
%     axis([0.3 0.9 -60 60]) 
%     pause(SignalParam.ChirpT)
end


figure
subplot(2,1,1)
plot(trajectory(:,1))
ylim([0 1])

subplot(2,1,2)
plot(trajectory(:,2))
ylim([-90 90])

%figure
%subplot(2,1,1)
%plot(hampel(trajectory(:,1),15,0.1))
%ylim([0.1 0.8])

%subplot(2,1,2)
%plot(hampel(trajectory(:,2),10))
%ylim([-90 90])
