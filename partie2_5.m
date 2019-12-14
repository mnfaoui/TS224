close all
clear
clc

%% Constantes
s=load('./projet_signal-master/fcno04fz');
s=s.fcno04fz;
Fe=8e3;

%% Declaration de variable

SNR       = 5;
lenWindow = 15e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
L         = 100;%125
threshold = 0.8e4;
windows   = window(@hamming,lenWindow);

%% Bruit
Vs  = (s'*s)/length(s);
Vb  = Vs/SNR;
signal = s + sqrt(Vb)*randn(size(s));

%% Fenetrage
[trameSansFenetre,trameAvecFenetre,signal] = getTrame(windows,signal);

%

%% Hankel

trame_rcv = trameAvecFenetre;
trame_new = zeros(size(trame_rcv));
nbTrame   = size(trame_rcv,2);
N = size(trame_rcv(:,1),1);

listValeurSing =zeros(1,size(trame_rcv,2));
M = N+1-L;

listU = zeros(L,L,nbTrame);
listS = zeros(L,M,nbTrame);
listV = zeros(M,M,nbTrame);

for noTrame=1:size(trame_rcv,2)
    H = hankel( trame_rcv(1:L,noTrame), trame_rcv(L:end,noTrame));
    [U,S,V] = svd(H);
    
    listU(:,:,noTrame) = U;
    listS(:,:,noTrame) = S;
    listV(:,:,noTrame) = V;
    
    dia = diag(S);
    listValeurSing(noTrame) = max(S(:));
    
    i = 7;
    if(S(1,1)>threshold)
%         tmp = dia(1:end-1)-dia(2:end);
%         [~,i] = max(dia(3:end)-2*dia(2:end-1)+dia(1:end-2));
%         i = i+2;
          
        X=S.*(S>S(i,i));
    else
        X=S/1000;
    end
    
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



figure
spectrogram(s,150,0,[],Fe,'yaxis')
figure
spectrogram(new_Signal,150,0,[],Fe,'yaxis')