function afficheSpectrogramme( s,windows,Fe)
    lenWindow = length(windows);
    spectrogram(s,windows,lenWindow/2,2^nextpow2(length(windows)),Fe,'yaxis')
end

