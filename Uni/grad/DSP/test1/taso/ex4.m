sig=loadbin('sigar-003.bin');
p = 12;

sig = sig.*hamming(length(sig));
[a_lpc, er] = aryule(sig, p);
a_lpc