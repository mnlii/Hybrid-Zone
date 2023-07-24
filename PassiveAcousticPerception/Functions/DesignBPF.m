function [b_value,a_value] = DesignBPF(F_Center,Filter_Pass,Filter_Stop,fs)
    Wp=2*pi*[F_Center - Filter_Pass/2,F_Center + Filter_Pass/2]/fs;
    Ws=2*pi*[F_Center - Filter_Stop/2,F_Center + Filter_Stop/2]/fs;
    Rp = 1;
    Rs = 30;
    [n,Wn] = ellipord(Wp/pi,Ws/pi,Rp,Rs);
    [b_value,a_value] = ellip(n,Rp,Rs,Wn);
%     if (FigureParam.BPFfig)
%         figure
%         W=-fs:1:fs;
%         [Hb,wb]=freqz(b_value,a_value,W,fs);
%         plot(wb,20*log10(abs(Hb)),'b');
%         xlabel('Hz');
%         ylabel('·ùÖµ/dB');
%     end
end