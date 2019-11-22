function [ densite2 ] = getWelch( signal,l)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
%     n=length(signal);
%     l=n/window;
%     a=(1-noverlap); 
%     densite = zeros(1,window);
%     start = 1;
%     
%     
% %     for ii=1:floor(1/a)
% %         ii
% %         densite(ii,:)=[densite;getpBartlett(signal(start:end),window)];
% %         start =(1-noverlap)*window*ii;
% %     end
%     
%     densite = sum(densite,1)/floor(1/a);
        k=length(signal)/l;
        densite = getpBartlett(signal,l);
        signal1 = [signal(l/2:end) zeros(1,l/2)];
        densite1= getpBartlett(signal1,l);
        densite2=(densite1+densite)/2;
            
end

