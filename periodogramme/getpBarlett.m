function [ tmp ] = getpBarlett( signal, k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    n = length(signal);
    a=mod(n,k);
    A=padarray(signal,a);
    l = floor(n/k);
    tmp = reshape(A,k,l);
    tmp = fftshift(fft(tmp,[],1)).^2/k;
    tmp = abs(sum(tmp,2)).'/l;
end

