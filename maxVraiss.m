close all
clear
clc

%% Constantes
addpath(genpath('fonction'))
s=load('./projet_signal-master/fcno03fz');
s=s.fcno03fz;
Fe=8e3;

%% Declaration de variable

SNR_dB    = -9;
SNR       = 10^(SNR_dB/10);
lenWindow = 30e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
L         = lenWindow*0.9;%
windows   = window(@hamming,lenWindow);

s = s(1:lenWindow*floor(length(s)/lenWindow));
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

N = size(trame_rcv(:,1),1); %% length windows
listValeurSing = zeros(1,size(trame_rcv,2));
listk          = zeros(1,size(trame_rcv,2));
M = N+1-L;

listU = zeros(L,M,nbTrame);
listS = zeros(M,M,nbTrame);
listV = zeros(M,M,nbTrame);


threshold = (trame_rcv(:,1)'*trame_rcv(:,1))*2;

for noTrame=1:size(trame_rcv,2)
    H = hankel( trame_rcv(1:L,noTrame), trame_rcv(L:end,noTrame));
    [U,S,V] = svd(H,'econ');
    
    listU(:,:,noTrame) = U;
    listS(:,:,noTrame) = S;
    listV(:,:,noTrame) = V;
    
    listValeurSing(noTrame) = max(S(1,1));
    
    X=S.*S-threshold;
    X=sqrt(X.*(X>0));
    listk(noTrame) = sum(sum(X>0)); 
    Hbis = U*X*V';
    
    a     = Hbis;
    [m,n] = size(a);
    idx   = hankel(1:m,m:(n-1)+m);
    out   = accumarray(idx(:),a(:),[],@mean); %% moyennage des anti diag
   
    trame_new(:,noTrame)= out;
end

%% Reconstitution du signal
signal_new = reconstructionSignal(trame_new,windows);

%% affichage
afficherTemps(s,signal,signal_new)
afficherSpetr(s,signal,signal_new,windows,Fe)

%% calcul SNR
bruit_new  = s-signal_new;
SNR_new    = var(s)/var(bruit_new);
SNR_dB_new = 10*log10(SNR_new);
