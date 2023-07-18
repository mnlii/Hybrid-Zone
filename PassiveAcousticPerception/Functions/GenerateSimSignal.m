FigureParam = InitFigParam();
SignalParam = InitSignalParam();
FilterParam = InitFilterParam(SignalParam);
MicroParam = InitMicroParam(SignalParam,FigureParam);

chirp = GenerateFMCWSignal(SignalParam);

audioPos=0;

frameNum=300;
v=0.1;
radv=v/0.1;
curpos=0.4;
pos=zeros(frameNum,1);
for idx=1:frameNum
    dirt=pi/2-radv*idx*SignalParam.ChirpT;
    curpos=curpos+exp(1i*dirt)*v*SignalParam.ChirpT;
    pos(idx)=curpos;
end

Signal=zeros(frameNum*SignalParam.ChirpSize+500,8);
for idx=1:300
    dis=abs(audioPos-MicroParam.Positions);
    delay=round(dis/SignalParam.Speed*SignalParam.SampleFrequency);
    for mic=1:6
        Signal(delay(mic)+(1:SignalParam.ChirpSize)+(idx-1)*SignalParam.ChirpSize,mic)=...
            Signal(delay(mic)+(1:SignalParam.ChirpSize)+(idx-1)*SignalParam.ChirpSize,mic)+...
            chirp*0.5;
    end

    if idx>75
        dis=abs(audioPos-pos(idx))+abs(pos(idx)-MicroParam.Positions);
        delay=round(dis/SignalParam.Speed*SignalParam.SampleFrequency);
        for mic=1:6
            Signal(delay(mic)+(1:SignalParam.ChirpSize)+(idx-1)*SignalParam.ChirpSize,mic)=...
                Signal(delay(mic)+(1:SignalParam.ChirpSize)+(idx-1)*SignalParam.ChirpSize,mic)+...
                chirp*0.05;
        end
    end
end

Signal=add_noisem(Signal,SignalParam,50);
audiowrite('.\Functions\source\Sim.wav',real(Signal),SignalParam.SampleFrequency)