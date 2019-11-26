clear 
close all
clc
Fe =8e3;

signal=load('./projet_signal-master/fcno04fz.mat');
signal=signal.fcno04fz;
n=length(signal);

RSB_dB = 5; %RSB_dB = 10*log(RSB) =  10*log(Es/Eb)
Es =(signal'*signal)/n;  %%signal'*signal;
Eb =  (Es/10^(RSB_dB/10));
sigma = sqrt(Eb);
bruit=sigma*randn(size(signal));

axeTmp = (1:n)/Fe;
signalb = signal + bruit;
figure;
subplot(2,1,1)
plot(axeTmp,signal)
xlim([0 7])
subplot(2,1,2)
spectrogram(signal,150,0,[],Fe,'yaxis')
xlim([0 7])

figure
subplot(2,1,1)
plot(axeTmp,signalb);
xlim([0 7])
subplot(2,1,2)
spectrogram(signalb,150,0,[],Fe,'yaxis')
xlim([0 7])