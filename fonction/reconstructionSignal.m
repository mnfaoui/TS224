function [signal_new] = reconstructionSignal(trame_new,windows)
% s_k(t) = x(t)W(t-kT);
% s_{k+1}(t) = x(t)W(t-(k+1)T);
% x(t) = (s_k(t)+s_{k+1}(t)) / ( W(t-kT)+W(t-(k+1)T) );

% trame = trame(:,1:end-1);
% trame = trame_new;
lenWindow = length(windows);
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
signal_new(signal_new==inf)  = 0;
signal_new(signal_new==-inf) = 0;
end

