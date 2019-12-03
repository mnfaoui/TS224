close all
clear 
clc 
%% Declaration de variable 
Fe=8e3;
signal=load('./projet_signal-master/fcno04fz');
signal=signal.fcno04fz;

lenWindow = 25e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
len = length(signal);
n = floor(len/lenWindow);

signal1 = signal(lenWindow/2+1:end);
len1 = length(signal1);
n1 = floor(len1/lenWindow);


%% Fenetrage

A = reshape(signal(1:n*lenWindow),lenWindow,n);
B = reshape(signal1(1:n1*lenWindow),lenWindow,n1);

C = zeros(lenWindow,n1+n);
C(:,1:2:end) = A;
C(:,2:2:end) = B; 
trame =C.*repmat(window(@hamming,lenWindow),1,n1+n) ;

%% Hankel

for nbTrame=1:size(trame,2)
    N = length(trame(:,1));
    L = round(1/2*N);
    M = N-L+1;
    H = hankel( trame(1:L,nbTrame), trame(1:M,nbTrame));
    [U,S,V] = svd(H);    
    X=S.*(S>S(1,1)*0.5);
    Hbis = U*X*V;
    tramebis=[Hbis(:,1)', Hbis(end,2:end)];
    trame(:,nbTrame) = tramebis';
end

% tramebis=[Hbis(:,1)', Hbis(end,2:end)];


%% Reconstitution du signal 
% s_k(t) = x(t)W(t-kT);
% s_{k+1}(t) = x(t)W(t-(k+1)T);
% x(t) = (s_k(t)+s_{k+1}(t)) / ( W(t-kT)+W(t-(k+1)T) );
trame = trame(:,1:end-1);
new_Signal = trame(:,1:2:end);
new_Signal = new_Signal(:);

new_Signal_Bis = trame(:,2:2:end);
new_Signal_Bis = new_Signal_Bis(:);

new_Signal(size(trame,1)/2+1:end-size(trame,1)/2) ...
    = new_Signal(size(trame,1)/2 +1:end-size(trame,1)/2) + new_Signal_Bis(1:end);

W = window(@hamming,size(trame,1)) ;
W = repmat(W,1,length(new_Signal)/length(W)) ;
W = W(:);
W(size(trame,1)/2+1:end-size(trame,1)/2) ...
    = W(size(trame,1)/2 +1:end-size(trame,1)/2) + W(1:end-200);

new_Signal = new_Signal./W;

figure,plot(new_Signal);
hold on
plot(signal);

