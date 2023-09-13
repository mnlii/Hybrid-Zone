function FilterParam = InitFilterParam(SignalParam)
    FilterParam.F_Center = 20000;
    FilterParam.BPF_Filter_Pass = 1000;
    FilterParam.BPF_Filter_Stop = 2000;
    FilterParam.HPF_Filter_Pass = SignalParam.ChirpF0 - 2000;
    FilterParam.HPF_Filter_Stop = SignalParam.ChirpF0 - 3000;
    FilterParam.LPF_Filter_Pass = 1300;
    FilterParam.LPF_Filter_Stop = 1500;
end