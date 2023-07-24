function Y = add_noisem(signal,SignalParam,SNR)
% add_noisem add determinated noise to a signal.
% X is signal, and its sample frequency is fs;
% filepath_name is NOISE's path and name, and the SNR is signal to noise ratio in dB.
load('babble.mat');
NOISE=resample(babble,SignalParam.SampleFrequency,19980);
nx=size(signal,1);
ny=size(signal,2);
NOISE=NOISE(1:nx);
NOISE=NOISE-mean(NOISE);
signal_power = 1/nx/ny*sum(sum(signal.*signal));
noise_variance = signal_power / ( 10^(SNR/10) );
NOISE=sqrt(noise_variance)/std(NOISE)*NOISE;
Y=signal+NOISE;