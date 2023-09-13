function SignalParam = InitSignalParam()
    SignalParam.ChirpBW = 4e3;% The bandwidth of Chirp (4kHz)
    SignalParam.ChirpT = 4e-2;
    SignalParam.ChirpF0 = 18e3;% The initial frequency of Chirp (18kHz)
    SignalParam.ChirpF1 = SignalParam.ChirpF0+SignalParam.ChirpF0;

    SignalParam.ScaleFactors = 1;% We calculate k chirp at the same time
    SignalParam.ChirpWindows = SignalParam.ChirpT*SignalParam.ScaleFactors;% The duration of Window for calculation
    SignalParam.SampleFrequency = 48e3;%48kHz
    SignalParam.ChannelNum = 6;% The number of microphones
    SignalParam.Speed = 346;% The speed of sound signals
    SignalParam.SampleTime = 1/SignalParam.SampleFrequency;% Sampling time for each sample

    SignalParam.ChirpSize = SignalParam.ChirpT *2/SignalParam.SampleTime; % Samples for each chirp
   
    SignalParam.SearchChirpRange= 300:30:1920;%40个子采样点

    %跟踪目标个数
    SignalParam.targetNum=1;

    SignalParam.FFTLength=48e3;
end