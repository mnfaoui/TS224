close all
clear
clc
fe = 8e3;
addpath('projet_signal-master/')
addpath('periodogramme/')
s1 = load('fcno03fz.mat');
s2 = load('fcno04fz.mat');

sigma = 1;
n = 10000;
% bruit = @(n,sigma) randn(1,n)*sigma;

bruit = randn(1,n)*sigma;
axe_frq = @(N) (-N/2:N/2-1)*1/N;
%% autocorrelation 
% theorique  
autocor_theo = zeros(1,n);
autocor_theo(round(n/2)) = sigma;

% estimateur biased
autocor_biased = xcorr(bruit,'biased');

% estimateur unbiased
autocor_unbiased = xcorr(bruit,'unbiased');

%% Affichage
% figure,plot(autocor_theo);
% figure,plot(autocor_biased);
% figure,plot(autocor_unbiased);


%% spectre de puissance 
% sdp = fftshift(abs(fft(bruit))).^2/n;
% sdp1=abs(fft(autocor_theo));
% figure,plot(sdp);
% figure,plot(sdp1);

%% Periodogramme
pdaniel = getpDaniel(bruit,n/100);
pbarlett= getpBartlett(bruit,n/10);
pwelch= getWelch(bruit,n/10);

subplot(311)
plot((1:length(bruit))/length(bruit),bruit);
xlabel('temps');
title('Brut Blanc Gaussien Centré')

subplot(312)
plot(axe_frq(length(pbarlett)),pbarlett);
title('Periodogramme de Bartlett')
xlabel('frequence');

subplot(313)
plot(axe_frq(length(pwelch)),pwelch);
title('Periodogramme de Welch')
xlabel('frequence');

%% Correlogramme

correlogramme = getCorrelogramme(bruit);
figure,plot(correlogramme);











