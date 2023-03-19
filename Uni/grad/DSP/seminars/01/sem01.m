% initialize
N = 1e4;
Npoints = N;
% Npoints = 1e3;
B = 300;
fs = 8e3;
wlen = min([512 N]);
% gain control over the current random number stream
stream = RandStream.getGlobalStream;

%% Generate a LF colored noise using an AR model
% generate a white noise signal with normed gaussian distribution
stream.reset;
un = randn(N, 1);  % zero EV and unitary variance

% determine system function - LP AR model
a1_abs = exp(-B*pi/fs);
% a1_abs = 1 - B/2;  % approximated logarithm formula

b = 1;
a = [1 -a1_abs];

% inspect the frequency characteristics of the system
figure(1)
freqz(b, a, Npoints, fs, "whole")

% generate the desired LF colored noise using the created AR model
nc1 = filter(b, a, un);
x = nc1;

% inspect the difference colored vs white noise
figure(2)
subplot(211);
plot(un)
title("White noise")

subplot(212);
plot(x)
title("LF colored noise generated using an AR model")

% inspect the unsmoothened periodogram PSD estimate
figure(3)
Gx = abs(fft(x(1:wlen)).^2)/wlen;
subplot(211);
plot(10*log10(Gx))
axis tight
title("Unsmoothened periodogram PSD estimate")

% and the smoothened Welch's PSD estimate
subplot(212);
pwelch(x, wlen, wlen/2, wlen, fs, 'twosided')
title("Smoothened Welch's PSD estimate")

% inspect the spectrogram of the signal
figure(4)
spectrogram(x, wlen, wlen/2, wlen, fs, "yaxis")
colorbar off, colormap jet

%% Generate a HF colored noise using an AR, resp. MA, model
% repeat the above task for HF noise generation using the same system

% generate a LW, resp. HF, colored noise by filtering white noise
stream.reset;
un = randn(N, 1);
% using an MA model of the same parameters
b1_abs = exp(-B*pi/fs);
b = [1 -b1_abs];
a = 1;

figure(5)
freqz(b, a, Npoints, fs, "whole")

nc2ma = filter(b, a, un);
x = nc2ma;
figure(6)
subplot(211), plot(un)
title("White noise")
subplot(212), plot(x)
title("HF colored noise generated using an MA model")

% inspect the unsmoothened periodogram PSD estimate
figure(7)
Gx = abs(fft(x(1:wlen)).^2)/wlen;
subplot(211);
plot(10*log10(Gx))
axis tight
title("Unsmoothened periodogram PSD estimate")

% and the smoothened Welch's PSD estimate
subplot(212), pwelch(x, wlen, wlen/2, wlen, fs, "twosided")
title("Smoothened Welch's PSD estimate")

% inspect the spectrogram of the signal
figure(8)
spectrogram(x, wlen, wlen/2, wlen, fs, "yaxis")
colorbar off, colormap jet

% missing LF signal nc1ma

%% Generate a bandlimited colored noise using a 2nd order AR model
stream.reset;
un = randn(N, 1);

a1 = [1 0.9*exp(+1j*0.9)];
a2 = [1 0.9*exp(-1j*0.9)];
b = 1;
a = conv(a1, a2);

figure(9)
freqz(b, a, Npoints, fs, "whole")

nc3 = filter(b, a, un);
x = nc3;
figure(10)
subplot(211), plot(un)
title("White noise")
subplot(212), plot(x)
title("Bandlimited colored noise generated using an MA model")

% inspect the unsmoothened periodogram PSD estimate
figure(11)
Gx = abs(fft(x(1:wlen)).^2)/wlen;
subplot(211);
plot(10*log10(Gx))
axis tight
title("Unsmoothened periodogram PSD estimate")

% and the smoothened Welch's PSD estimate
subplot(212), pwelch(x, wlen, wlen/2, wlen, fs, "twosided")
title("Smoothened Welch's PSD estimate")

% inspect the spectrogram of the signal
figure(12)
spectrogram(x, wlen, wlen/2, wlen, fs, "yaxis")
colorbar off, colormap jet
