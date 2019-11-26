function [ corr ] = getCorrelogramme( signal )
%UNTITLED Summary of this function goes here
%   Corrélogramme :transformée de Fourier del'auto covariance
n=length(signal);
Rxx = xcorr(signal,'biased' );
corr = fftshift(abs(fft(Rxx)));
end

