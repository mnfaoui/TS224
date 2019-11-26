close all
clear 
clc 
%% Declaration de variable 
Fe=8e3;
signal=load('./projet_signal-master/fcno03fz');
signal=signal.fcno03fz;

lenWindow = 25e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
len = length(signal);
n = floor(len/lenWindow);
signal1 = signal(lenWindow/2:end);
len1 = length(signal1);
n1 = floor(len1/lenWindow);
%% Fenetrage

A = reshape(signal(1:n*lenWindow),lenWindow,n);
B = reshape(signal1(1:n1*lenWindow),lenWindow,n1);

C = zeros(lenWindow,n1+n);
C(:,1:2:end) = A;
C(:,2:2:end) = B; 
trame =C.*repmat(window(@hamming,lenWindow),1,n1+n) ;


