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
x = nc1;

figure(2)
subplot(211), plot(un)
subplot(212), plot(nc1)

% PST via Welch's method
wlen = min([512 N]);
figure(3)
pwelch(nc1,wlen,wlen/2,wlen,fs,'twosided')

% PST from definition
figure(4)
Gx = abs(fft(x(1:wlen)).^2)/wlen;
plot(10*log10(Gx).^2/wlen)  % xlim([0 wlen/2])
axis tight  % fit figure (doesn't work with xlim)

% spectrogram
figure(5)
spectrogram(nc1,wlen,wlen/2,wlen,fs,'yaxis')
colorbar off, colormap jet

% plot poles and zeros using zplane func
figure(6)
zplane(b,a)

%%

p = 1;
a_lpc = lpc(x,p);

r = roots(a_lpc);

figure(7)
hold on
plot(real(r),imag(r), 'rx')
hold off

figure(8)
[H,ff] = freqz(1,a_lpc,Npoints,fs,'half');
hold on
plot(ff,20*log10(abs(H)), 'r')
hold off

%%

x = loadbin('nc4.bin');
fs = 8000;

figure(11)
pwelch(x,wlen,wlen/2,wlen,fs,"onesided")  % from the figure, we think it's a 4th order filter (two peaks)

p = 4;
a_lpc = lpc(x,p);
figure(12)
freqz(1,a_lpc)

figure(13)
zplane(1,a_lpc)

%%

x = loadbin('vm0-512.bin');
fs = 16e3;

x = x(1:512);  % just to be sure we're working with 512 samples
w = hamming(512);
x = x.*w;  % weighting the signal using Hamming's window

figure(21)
plot(x)

figure(22)
wlen = length(x);
Gx = abs(fft(x).^2)/wlen;  % two-sided spectrum computation
plot(10*log10(Gx(1:wlen/2))), axis tight  % select just the right half with (1:wlen/2)

p = 16;
[aa,Ep] = lpc(x,p);

% several commmon visual errors which might occur in the following
% a) the axes are mismatched (scaling is off) - watch out for the number
%   of points in which I calculate each thing (wlen/2 something)
% b) vertical shift - maybe forgotten scaling by energy
% c) different dynamic of each graph - we have to plot |H(z)|^2 for both
%   -> fixed by checking that (has to be 20log10 and not 10log10)
figure(23)
Hz = freqz(sqrt(Ep),aa,wlen/2);  % Ep scaling to adjust the height to fit
hold on
plot(20*log10(abs(Hz)))
hold off