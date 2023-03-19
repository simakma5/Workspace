fs = 200;
mypoles = [0.8*exp(j*pi/4) 0.8*exp(-j*pi/4) 0.95*exp(j*5*pi/6) 0.95*exp(-j*5*pi/6) ];
a = poly(mypoles);
b = 1;
wlen = 512;
nfft = 512;

un = loadbin('un-excitation.bin');
nc1 = filter(b, a, un);
psd = pwelch(nc1, hamming(wlen), wlen/2, nfft, fs, 'twosided', 'psd');
10*log10(psd(1:8))