function FigureParam = InitFigParam()
    FigureParam.MicMap = false;% 是否绘制麦克风位置图
    FigureParam.BPFfig = false;% 是否绘制带通滤波器功率图
    FigureParam.HPFfig = false;% 是否绘制高通滤波器功率图
    FigureParam.LPFfig = false;% 是否绘制低通滤波器功率图
    FigureParam.MovePath = false;% 是否绘制轨迹图
    FigureParam.Heatfig = false;% 是否绘制联合估计时的热图
    FigureParam.Directfig = false;% 是否绘制直接路径消除3D图
    FigureParam.Rebuildfig = false;% 是否绘制重建信号时的拟合图
    FigureParam.EstimateMovefig = false;% 是否绘制估计轨迹图
end