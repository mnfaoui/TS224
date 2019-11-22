function [ densite ] = getpDaniel(signal,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    autocor_unbiased = fftshift(abs(fft(signal))).^2/length(signal);
    bis1 = signal(end-n+1:end);
    bis2 = signal(1:n);
    
    densite = conv([bis1 autocor_unbiased bis2], ones(1,n))/n;
    densite = densite(n+n/2+1:end-3*n/2+1);

end

