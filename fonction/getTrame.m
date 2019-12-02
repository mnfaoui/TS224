function [trame] = getTrame(window,signal)
%UNTITLED Summary of this function goes here
%   exemple : [trame] = getTrame(window(@hamming,lenWindow),signal)
lenWindow = length(window);

A = reshape(signal(1:n*lenWindow),lenWindow,n);
B = reshape(signal1(1:n1*lenWindow),lenWindow,n1);

C = zeros(lenWindow,n1+n);
C(:,1:2:end) = A;
C(:,2:2:end) = B; 
trame =C.*repmat(window,1,n1+n) ;
end

