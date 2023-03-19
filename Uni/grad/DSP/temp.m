%% sem1task1
N = 1e4;
fs = 8e3;
B = 3e2;
wlen = min([512 N]);
Rstream = RandStream.getGlobalStream;

%
Rstream.reset
u = randn(N,1);
% u = u.*hamming(length(u));

%
a1 = exp(-B*pi/fs);
b = 1;
a = [1 -a1];
% H(z) = 1/(1+a1mod*z^(-1));

%
figure(111)
freqz(b,a,N)

%
nc1 = filter(b,a,u);

figure(112)
subplot(311);
plot(nc1)
title('The filtered signal nc1')

subplot(312);
pwelch(nc1,wlen,wlen/2,wlen,fs,'twosided')

subplot(313);
Gx = abs(fft(nc1(1:wlen))).^2./wlen;
plot(10*log10(Gx))
axis tight
title('Periodogram (unsmoothened PSD estimate)')

figure(113)
spectrogram(nc1,wlen,wlen/2,wlen,'yaxis')
colormap jet

%% sem1task2
N = 1e4;
fs = 8e3;
B = 3e2;
wlen = min([512 N]);
Rstream = RandStream.getGlobalStream;

%
Rstream.reset
u = randn(N,1);
% u = u.*hamming(length(u));

%
a1 = exp(-B*pi/fs);
b = 1;
a = [1 a1];
% H(z) = 1/(1-a1mod*z^(-1));

%
figure(121)
freqz(b,a,N)

%
nc1 = filter(b,a,u);

figure(122)
subplot(311);
plot(nc1)
title('The filtered signal nc1')

subplot(312);
pwelch(nc1,wlen,wlen/2,wlen,fs,'twosided')

subplot(313);
Gx = abs(fft(nc1(1:wlen))).^2./wlen;
plot(10*log10(Gx))
axis tight
title('Periodogram (unsmoothened PSD estimate)')

figure(123)
spectrogram(nc1,wlen,wlen/2,wlen,'yaxis')
colormap jet

%% sem1task3
N = 1e4;
fs = 8e3;
wlen = min([512 N]);
Rstream = RandStream.getGlobalStream;

%
Rstream.reset
u = randn(N,1);
% u = u.*hamming(length(u));

%
p1 = 0.9*exp(+1i*0.9);
p2 = 0.9*exp(-1i*0.9);
b = 1;
a = conv([1 p1],[1 p2]);
% H(z) = 1/((1+p1*z^(-1))(1+p2*z^(-1)));

%
figure(131)
freqz(b,a,N)

%
nc3 = filter(b,a,u);

figure(132)
subplot(311);
plot(nc3)
title('The filtered signal nc1')

subplot(312);
pwelch(nc3,wlen,wlen/2,wlen,fs,'twosided')

subplot(313);
Gx = abs(fft(nc3(1:wlen))).^2./wlen;
plot(10*log10(Gx))
axis tight
title('Periodogram (unsmoothened PSD estimate)')

figure(133)
spectrogram(nc3,wlen,wlen/2,wlen,'yaxis')
colormap jet

%% sem2task1
N = 1e4;
fs = 8e3;
wlen = 512;
nc4 = loadbin('nc4.bin');

%
figure(211)
subplot(211);
plot(nc4)
subplot(212);
pwelch(nc4,wlen,wlen/2,wlen,fs)
% two peaks -> 4th order

%
[lpc1,Ep1] = lpc(nc1,1);
[lpc2,Ep2] = lpc(nc2,1);
[lpc3,Ep3] = lpc(nc3,2);
[lpc4,Ep4] = lpc(nc4,4);

%
figure(212)
subplot(221);
zplane(1,lpc1)
title('Pole-Zero Plot nc1')

subplot(222);
zplane(1,lpc2)
title('Pole-Zero Plot nc2')

subplot(223);
zplane(1,lpc3)
title('Pole-Zero Plot nc3')

subplot(224);
zplane(1,lpc4)
title('Pole-Zero Plot nc4')

%% sem2task2
fs = 16e3;
x = loadbin('vm0-512.bin');
x = x.*hamming(512);
N = length(x);

%
p = 30;
[a,Ep] = lpc(x,p);

%
Gx = abs(fft(x).^2)./N;
H = freqz(sqrt(Ep),a,N/2);
LPCspec = abs(H).^2;

figure(221)
subplot(311);
plot(real(x))
title('Waveform')
axis tight

subplot(312);
plot(10*log10(Gx(1:N/2)))
title('Periodogram')
axis tight

subplot(313);
plot(10*log10(LPCspec))
title('Frequency characteristic of the LPC spectrum')

%
figure(222)
plot(10*log10(Gx(1:N/2)))
hold on
plot(10*log10(LPCspec),'r')
hold off
