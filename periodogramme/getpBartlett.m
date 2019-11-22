function [ tmp ] = getpBartlett( signal, k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    n = length(signal);
    a=mod(n,k);
    
    signal = signal(1:end-a);
    l = length(signal)/k;
    tmp = reshape(signal,k,l);
    tmp = abs(fftshift(fft(tmp,k,1))).^2;
    tmp = sum(tmp,2).'/(k*l);
end

