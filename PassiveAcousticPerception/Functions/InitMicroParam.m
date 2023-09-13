function MicroParam = InitMicroParam(SignalParam,FigureParam)
    % 以 x正向轴上的mic1为原点

    MicroParam.MicroTotal = SignalParam.ChannelNum;
    % 麦克风任一点到中心的距离（单位：米）
    MicroParam.toEdge = 0.0475;
    % 麦克风设备仿真的中心位置
    MicroParam.Position = 0;
    % 相邻麦克风夹角差值
    MicroParam.Theta = 2 * pi / MicroParam.MicroTotal;
    % 麦克风阵列各个麦克风的位置
    MicroParam.Angles=-MicroParam.Theta*(0:1:MicroParam.MicroTotal-1)/pi*180+210;
    MicroParam.Positions = MicroParam.toEdge * exp(1i*MicroParam.Angles/180*pi) + MicroParam.Position;
    
    if (FigureParam.MicMap)
        scatter(real(MicroParam.Positions),imag(MicroParam.Positions));
        axis([-0.5,0.5,-0.5,0.5]);
        title('MicroPhone Layout');
        grid on;
    end
end