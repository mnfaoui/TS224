close all
clear
clc

%% Constantes
addpath('fonction')
s=load('./projet_signal-master/fcno03fz');
s=s.fcno03fz;
Fe=8e3;

%% Declaration de variable

SNR_dB    = -5;
SNR       = 10^(SNR_dB/10);
lenWindow = 30e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
L         = lenWindow*0.9;%
windows   = window(@hamming,lenWindow);
% windows   = ones(lenWindow,1);
%% Bruit
Vs  = var(s);%(s'*s)/length(s);
Vb  = Vs/SNR;
signal = s + sqrt(Vb)*randn(size(s));

%% Fenetrage
[trameSansFenetre,trameAvecFenetre,signal] = getTrame(windows,signal);

%% Hankel

trame_rcv = trameAvecFenetre;
trame_new = zeros(size(trame_rcv));
nbTrame   = size(trame_rcv,2);
N = size(trame_rcv(:,1),1);

listValeurSing =zeros(1,size(trame_rcv,2));
M = N+1-L;

listU = zeros(L,M,nbTrame);
listS = zeros(M,M,nbTrame);
listV = zeros(M,M,nbTrame);



sigmaNoise = 2*Vb*(windows'*windows);

for noTrame=1:size(trame_rcv,2)
    H = hankel( trame_rcv(1:L,noTrame), trame_rcv(L:end,noTrame));
    [U,S,V] = svd(H,'econ');
    
    listU(:,:,noTrame) = U;
    listS(:,:,noTrame) = S;
    listV(:,:,noTrame) = V;
    
    dia = diag(S);
    listValeurSing(noTrame) = max(S(:));
    
    X = diag( dia.*dia - sigmaNoise);
    X = X.*(X>0);
    X = sqrt(X);
    Hbis = U*X*V';
    
    a     = Hbis;
    [m,n] = size(a);
    idx   = hankel(1:m,m:(n-1)+m);
    out   = accumarray(idx(:),a(:),[],@mean);
   
    trame_new(:,noTrame)= out;
end

%% Reconstitution du signal
% s_k(t) = x(t)W(t-kT);
% s_{k+1}(t) = x(t)W(t-(k+1)T);
% x(t) = (s_k(t)+s_{k+1}(t)) / ( W(t-kT)+W(t-(k+1)T) );

% trame = trame(:,1:end-1);
% trame = trame_new;
trame = trame_new;
new_Signal = trame(:,1:2:end);
new_Signal = new_Signal(:);

new_Signal_Bis = trame(:,2:2:end);
new_Signal_Bis = new_Signal_Bis(:);

new_Signal(size(trame,1)/2+1:end-size(trame,1)/2) ...
    =  new_Signal(size(trame,1)/2 +1:end-size(trame,1)/2)...
     + new_Signal_Bis(1:end);

W = windows;
W = repmat(W,1,length(new_Signal)/length(W)) ;
W = W(:);
W(size(trame,1)/2+1:end-size(trame,1)/2) ...
    = W(size(trame,1)/2 +1:end-size(trame,1)/2) + W(1:end-lenWindow);

signal_new = new_Signal./W;


figure,hold on
plot(signal);
plot(signal_new);

figure,hold on
plot(signal_new);
plot(signal);

figure,hold on
plot(signal_new);
plot(s);
%%
figure
spectrogram(s,windows,lenWindow/2,[],Fe,'yaxis')
figure
spectrogram(new_Signal,windows,lenWindow/2,[],Fe,'yaxis')

%% Spectrogramme
Specto = @(entry)afficheSpectrogramme(entry,windows,Fe);
figure
subplot(3,1,1)
Specto(s);
subplot(3,1,2)
Specto(signal);
subplot(3,1,3)
Specto(new_Signal);
