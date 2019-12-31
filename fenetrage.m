clc
close all
clear

Fe = 8000;
duree_fenetre = 30e-3;

N = Fe*duree_fenetre; 
NFFT = 2^(nextpow2(N)+5);

w_p     = window(@rectwin,N);
w_t     = window(@triang,N);
w_hann  = window(@hann,N);
w_hamm  = window(@hamming,N);

W_p = fftshift(abs(fft(w_p,NFFT)));
W_t = fftshift(abs(fft(w_t,NFFT)));
W_hann = fftshift(abs(fft(w_hann,NFFT)));
W_hamm = fftshift(abs(fft(w_hamm,NFFT)));
axe_freq = (-NFFT/2:NFFT/2-1)*Fe/NFFT;

afficher = @(s,title)plot((1:length(s))/Fe*1e3,s,'DisplayName',title);

figure(1)
hold on
afficher(w_p,'porte')
afficher(w_t,'triangulaire')
afficher(w_hann,'hann')
afficher(w_hamm,'hamming')
xlabel('Temps (ms)')
grid on
legend()

%%
figure(2)
semilogx(axe_freq,10*log10(W_p))
hold on
semilogx(axe_freq,10*log10(W_t))
semilogx(axe_freq,10*log10(W_hann))
semilogx(axe_freq,10*log10(W_hamm))
legend('porte','triangulaire','hann','hamming','Location','best')
xlim([0 4e2])
ylabel('amplitude(dB)')
xlabel('Fréquence normalisée')
grid on
