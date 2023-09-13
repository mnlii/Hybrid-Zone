clear;
ChirpBW = 4e3;% The bandwidth of Chirp (4kHz)
    ChirpT = 4e-2;
    ChirpF0 = 18e3;% The initial frequency of Chirp (18kHz)
    ChirpF1 = ChirpF0+ChirpF0;

    ScaleFactors = 1;% We calculate k chirp at the same time
    ChirpWindows = ChirpT*ScaleFactors;% The duration of Window for calculation
    SampleFrequency = 48e3;%48kHz
    ChannelNum = 6;% The number of microphones
    Speed = 346;% The speed of sound signals
    SampleTime = 1/SampleFrequency;% Sampling time for each sample

    ChirpSize = ChirpT*2/SampleTime; % Samples for each chirp
   
    SearchChirpRange= 300:30:1920;%40个子采样点

    %跟踪目标个数
    targetNum=1;

    FFTLength=48e3;




t_up = (0:SampleTime:ChirpT - SampleTime)';
chirp_up = exp(1i*2*pi*(ChirpF0*t_up + ChirpBW/ChirpT/2*t_up.^2));
t_down = (ChirpT:SampleTime:2*ChirpT - SampleTime)';
chirp_down = exp(1i*2*pi*((ChirpF0+ChirpBW*2)*t_down - ChirpBW/ChirpT/2*t_down.^2));  % 注意这里的负号
chirp_combined = [real(chirp_up); real(chirp_down)];
%b = cos(2*pi*ChirpF0*t);
%chirp_combined = real(b) + chirp_combined;
%chirp = chirp_combined;
s=cos(2*pi*600*(1:SampleFrequency)/SampleFrequency)'*0.05;
s(9601:19200)=0;
s(28801:38400)=0;
s = [s; repmat(chirp_combined, round(300/(2*ChirpT)), 1)];  % 注意这里的2*ChirpT，因为我们现在有两个chirp，所以总时间是原来的两倍
s = s / max(abs(s));
%plot(s(1:3000))
audiowrite("Functions\chirp4.wav",s,SampleFrequency); 