function afficher(signalSansBruit,signalAvecBruit,signalDebruite)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure
subplot(3,1,1);
plot(signalSansBruit);
xlim([1 length(signalSansBruit)])
title('signal sans bruit')

subplot(3,1,2);
plot(signalAvecBruit);
xlim([1 length(signalAvecBruit)])
title('signal avec bruit')

subplot(3,1,3);
plot(signalDebruite);
xlim([1 length(signalDebruite)])
title('signal debruiter')

end

