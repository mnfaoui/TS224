close all
clear
clc

%% Constantes
addpath(genpath('fonction'))
addpath(genpath('projet_signal-master'))
s=load('fcno03fz');
s=s.fcno03fz;
Fe=8e3;

%% Declaration de variable

SNR_dB    = 1;
SNR       = 10^(SNR_dB/10);
lenWindow = 20e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
L         = lenWindow*0.70;%
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

thresholdTrue = sqrt(2*(windows'*windows)*Vb);


% threshold = sqrt(2*(windows'*windows)*Vb)*1.4;

threshold = sqrt(2*(trame_rcv(:,3)'*trame_rcv(:,3)))*1.28;
for noTrame=1:size(trame_rcv,2)
    H = hankel( trame_rcv(1:L,noTrame), trame_rcv(L:end,noTrame));
    [U,S,V] = svd(H,'econ');
    
    listU(:,:,noTrame) = U;
    listS(:,:,noTrame) = S;
    listV(:,:,noTrame) = V;
    
    listValeurSing(noTrame) = max(S(1,1));
    
    X=S.*(S>threshold);
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

axe_temp = @(n,F) (0:n-1)/Fe;
figure
subplot(3,1,1);
plot(axe_temp(length(s),Fe),s);
xlim([1 max(axe_temp(length(s)))] )
title('Signal sans bruit')

subplot(3,1,2);
plot(axe_temp(length(s),Fe),signal);
xlim([1 max(axe_temp(length(s)))] )
title('Signal avec bruit')

subplot(3,1,3);
tmp_axe = (1:lenWindow/2:lenWindow/2*length(listValeurSing))/Fe;
plot(tmp_axe,listValeurSing,'DisplayName','valeur singuliere max');
hold on
plot( tmp_axe,threshold*ones(size(listValeurSing)),'DisplayName','seuil');
hold off
title('Evolution de la valeur singuliere max en fonction du temps')
xlim([1 max(axe_temp(length(s)))] )
xlabel('temps (s)')
legend();
%%
afficherTemps(s,signal,signal_new)
afficherSpetr(s,signal,signal_new,windows,Fe)

%% calcul SNR
bruit_new  = s-signal_new;
SNR_new    = var(s)/var(bruit_new);
SNR_dB_new = 10*log10(SNR_new);
