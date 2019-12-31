function afficherSpetr(s,signal,new_Signal,windows,Fe)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Specto = @(entry)afficheSpectrogramme(entry,windows,Fe);
figure
subplot(3,1,1)
Specto(s);
title('Signal sans bruit')

subplot(3,1,2)
Specto(signal);
title('Signal avec bruit')

subplot(3,1,3)
Specto(new_Signal);
title('Signal debruiter')

end

