function [trameSansFenetre,trameAvecFenetre,signal] = getTrame(window,signal)
%UNTITLED Summary of this function goes here
%   exemple : [trame] = getTrame(window(@hamming,lenWindow),signal)

lenWindow = length(window);
len       = length(signal);
n         = floor(len/lenWindow);
signal    = signal(1:n*lenWindow);

signalbis = signal(lenWindow/2+1:end);
lenbis = length(signalbis);
n1 = floor(lenbis/lenWindow);
signalbis = signalbis(1:n1*lenWindow);

A = reshape(signal,lenWindow,n);
B = reshape(signalbis,lenWindow,n1);

C = zeros(lenWindow,n1+n);
C(:,1:2:end) = A;
C(:,2:2:end) = B;
trame =C.*repmat(window,1,n1+n) ;

trameSansFenetre = C;
trameAvecFenetre = trame;
end

