function chirp = GenerateFMCWSignal(SignalParam)
    t = (0:SignalParam.SampleTime:SignalParam.ChirpT - SignalParam.SampleTime)';
   %chirp = exp(1i*2*pi*(SignalParam.ChirpF0*t + SignalParam.ChirpBW/SignalParam.ChirpT/2*t.^2));
    %chirp = exp(1i*2*pi*((SignalParam.ChirpF0+SignalParam.ChirpBW)*t - SignalParam.ChirpBW/SignalParam.ChirpT/2*t.^2));
    % 
    % 
    % s=cos(2*pi*600*(1:SignalParam.SampleFrequency)/SignalParam.SampleFrequency)'*0.05;
    % s(9601:19200)=0;
    % s(28801:38400)=0;
    % s=[s;repmat(real(chirp),round(300/SignalParam.ChirpT),1)];
    % audiowrite("Functions\chirp.wav",s,SignalParam.SampleFrequency);


    t_up = (0:SignalParam.SampleTime:SignalParam.ChirpT - SignalParam.SampleTime)';
    chirp_up = exp(1i*2*pi*(SignalParam.ChirpF0*t_up + SignalParam.ChirpBW/SignalParam.ChirpT/2*t_up.^2));
    t_down = (SignalParam.ChirpT:SignalParam.SampleTime:2*SignalParam.ChirpT - SignalParam.SampleTime)';
    chirp_down = exp(1i*2*pi*((SignalParam.ChirpF0+2*SignalParam.ChirpBW)*t_down - SignalParam.ChirpBW/SignalParam.ChirpT/2*t_down.^2));  % 注意这里的负号
    chirp_combined = [chirp_up; chirp_down];
    % 
    % %b = cos(2*pi*SignalParam.ChirpF0*t);
    % %chirp_combined = real(b) + chirp_combined;
    chirp = chirp_combined;
    %fprintf('chirp size: %d x %d\n', size(chirp, 1), size(chirp, 2));


     %s=cos(2*pi*600*(1:SignalParam.SampleFrequency)/SignalParam.SampleFrequency)'*0.05;
     % s(9601:19200)=0;
     % s(28801:38400)=0;
     % s = [s; repmat(chirp_combined, round(300/(2*SignalParam.ChirpT)), 1)];  % 注意这里的2*SignalParam.ChirpT，因为我们现在有两个chirp，所以总时间是原来的两倍
     % s = s / max(abs(s));
    %plot(s(1:3000))
    % audiowrite("Functions\chirp4.wav",s,SignalParam.SampleFrequency); 

