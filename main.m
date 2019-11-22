close all
clear
clc

addpath('projet_signal-master/')
addpath('periodogramme/')
s1 = load('fcno03fz.mat');
s2 = load('fcno04fz.mat');

sigma = 1;
n = 10000;
% bruit = @(n,sigma) randn(1,n)*sigma;

bruit = randn(1,n)*sigma;

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
figure,plot(bruit);
figure,plot(pdaniel);
figure,plot(pbarlett);
figure,plot(pwelch);

   pxx = pwelch(bruit)