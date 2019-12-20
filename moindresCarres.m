close all
clear
clc

%% Constantes
addpath(genpath('fonction'))
s=load('./projet_signal-master/fcno03fz');
s=s.fcno03fz;
Fe=8e3;

%% Declaration de variable

SNR_dB    = 5;
SNR       = 10^(SNR_dB/10);
lenWindow = 30e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
L         = lenWindow*0.5;%
windows   = window(@hamming,lenWindow);
% windows   = ones(lenWindow,1);
%% Bruit
Vs  = var(s);%(s'*s)/length(s);
Vb  = Vs/SNR;
signal = s + sqrt(Vb)*randn(size(s));

%% Fenetrage
[trameSansFenetre,trameAvecFenetre,signal] = getTrame(windows,signal);

%

%% Hankelthreshold = sqrt(2*(windows'*windows)*Vb);

trame_rcv = trameAvecFenetre;
trame_new = zeros(size(trame_rcv));
nbTrame   = size(trame_rcv,2);
N = size(trame_rcv(:,1),1);

listValeurSing =zeros(1,size(trame_rcv,2));
M = N+1-L;

listU = zeros(L,L,nbTrame);
listS = zeros(L,M,nbTrame);
listV = zeros(M,M,nbTrame);

thresholdTrue = sqrt(2*Vb*L);
threshold = sqrt(2*(windows'*windows)*Vb);
% threshold = sqrt(2*var(trame_rcv(:,2)./windows)*L);
for noTrame=1:size(trame_rcv,2)
    H = hankel( trame_rcv(1:L,noTrame), trame_rcv(L:end,noTrame));
    [U,S,V] = svd(H);
    
    listU(:,:,noTrame) = U;
    listS(:,:,noTrame) = S;
    listV(:,:,noTrame) = V;
    
    dia = diag(S);
    listValeurSing(noTrame) = max(S(:));
     
    X=S.*(S>threshold);
    Hbis = U*X*V';
    
    a     = Hbis;
    [m,n] = size(a);
    idx   = hankel(1:m,m:(n-1)+m);
    out   = accumarray(idx(:),a(:),[],@mean);
   
    trame_new(:,noTrame)= out;
%     plot(out);

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

%%
figure,plot(listValeurSing);
hold on
plot(threshold*ones(length(listValeurSing)));
afficher(s,signal,signal_new)

%% Spectrogramme
Specto = @(entry)afficheSpectrogramme(entry,windows,Fe);
figure
subplot(3,1,1)
Specto(s);
subplot(3,1,2)
Specto(signal);
subplot(3,1,3)
Specto(new_Signal);


%%
figure
plot(new_Signal-signal)
figure
plot(signal)