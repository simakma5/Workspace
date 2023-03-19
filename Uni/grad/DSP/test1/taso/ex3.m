frame=loadbin('frame-006.bin');
s_w  = frame.*hamming(length(frame));
S = 10*log10((abs(fft(s_w)).^2)./length(frame));
S(1:8)