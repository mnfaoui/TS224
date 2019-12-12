close all
clear 
clc 


%% Declaration de variable 
Fe=8e3;
s=load('./projet_signal-master/fcno04fz');
s=s.fcno04fz;
SNR = 0.25;
lenWindow = 25e-3*Fe;  %% un signal de parole peut être modelise comme un signal quaso-stationnaire sur un intervalle de temps de 25 ms
L = 190;%125
threshold  = 0.5e5;
%% Bruit
Vs  = (s'*s)/length(s);
Vb  = Vs/SNR;
signal = s + sqrt(Vb)*randn(size(s));

%% Fenetrage

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
trame =C.*repmat(window(@hamming,lenWindow),1,n1+n) ;

%% Hankel
trame_rcv = trame;
N = length(trame(:,1));
M = N-L+1;
tmptest =zeros(1,size(trame_rcv,2));

for nbTrame=1:size(trame_rcv,2)
        H = hankel( trame(1:L,nbTrame), trame(L:end,nbTrame));
        [U,S,V] = svd(H);
        dia = diag(S);
        
        tmptest(nbTrame) = max(S(:));
%         figure(1)
%         plot(tmp);
%         figure(2)
%         plot(S)
%         [~,i] = max(dia(1:end-1)-dia(2:end));
%         X=S.*(S>S(i,i));
        if(tmptest(nbTrame)>threshold)
            tmp = dia(1:end-1)-dia(2:end);
            [~,i] = max(dia(3:end)-2*dia(2:end-1)+dia(1:end-2));
%             i = i+2;
            X=S.*(S>S(i,i));
        else
            X = 0*S;
        end
        Hbis = U*X*V';
        
        a = Hbis;
        [m,n] = size(a);
        idx = hankel(1:m,m:(n-1)+m);
        out = accumarray(idx(:),a(:))./[(1:min(n,m)),20*ones(1,abs(n-m)),(min(n,m)-1:-1:1)]';
        %%tramebis=[Hbis(:,1)', Hbis(end,2:end)];
%         sum(out == trame(:,nbTrame))
        trame(:,nbTrame)= out;
end

%% Reconstitution du signal 
% s_k(t) = x(t)W(t-kT);
% s_{k+1}(t) = x(t)W(t-(k+1)T);
% x(t) = (s_k(t)+s_{k+1}(t)) / ( W(t-kT)+W(t-(k+1)T) );

% trame = trame(:,1:end-1);
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
    = W(size(trame,1)/2 +1:end-size(trame,1)/2) + W(1:end-lenWindow);

new_Signal = new_Signal./W;

figure,hold on
plot(new_Signal);
plot(s)
plot(new_Signal);

figure,hold on
plot(new_Signal);
plot(signal);
% figure


figure
spectrogram(s,150,0,[],Fe,'yaxis')
figure
spectrogram(new_Signal,150,0,[],Fe,'yaxis')