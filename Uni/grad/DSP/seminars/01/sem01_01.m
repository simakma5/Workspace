% init
Npoints = 1000;
B = 300;
fs = 8000;

a1mod = exp(-B*pi/fs);

b = 1;
a = [1 -a1mod];

% frequency response
figure(1)
freqz(b,a,Npoints,fs,'whole')

% white and colored noise
N = 10000;
un = randn(N,1);

nc1 = filter(b,a,un);

figure(2)
subplot(211), plot(un)
subplot(212), plot(nc1)

% PST via Welch's method
wlen = 512;
figure(3)
pwelch(nc1,wlen,wlen/2,wlen,fs,'twosided')

% PST from definition
figure(4)
plot(10*log10(abs(fft(nc1(1:wlen))).^2/wlen))%, xlim([0 wlen/2])
axis tight  % fit figure (doesn't work with xlim)

% spectrogram
figure(5)
spectrogram(nc1,wlen,wlen/2,wlen,fs,'yaxis')
colorbar off, colormap jet

% plot poles and zeros using zplane func
figure(6)
zplane(b,a)