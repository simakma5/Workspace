close all; clear; clc

%%
R = 3; % range of target
vr = 5; % speed of target
fc = 77e9; % carrier frequency (common value for modern automotive radars)
BW = 4e9; % bandwidth
ksw = 100e12; % slope of chirp (MHz/s)
fs = 10e6; % sampling frequency
Nch = 32; % number of chirps in frame
c0 = 299792458; % speed of light

tau = 2*R/c0; % round-trip delay
Tch = BW/ksw; % length of single chirp
fb = tau*ksw; % beat frequency
fD = 2*vr*fc/c0; % Doppler shift
Ns = ceil(Tch*fs) + 1; % number of samples per chirp (usually a power of 2)
RMax = fs*c0/(4*ksw); % maximal range
dR = c0/(2*BW); % range resolution
RAxis = linspace(-RMax, RMax-dR, Ns); % range axis
dv = c0/(4*Tch*Nch*fc); % speed resolution
vrmax = c0/(4*fc*Tch); % maximal speed
vAxis = linspace(-vrmax, vrmax-dv, Nch); % speed axis
t = linspace(0, (Ns - 1)/fs, Ns).'; % time vector
sRx = exp(1j*(2*pi*fb*t + 2*pi*fD*(0:Nch-1)*Tch)); % Rx signal (some crazy Matlab shenanigans here)

% FFT: can be done in various ways in order not to obtain messy results
% further, FFT target is often windowed to achieve that
% further further, the window in one dimension can be different from the
% one used on the second dimension, again to improve the result purity

% SRx = fftshift(fft(fft(sRx, [], 1), [], 2));
% SRx = fftshift(fft(fft(sRx).').');
SRx = fftshift(fft2(sRx));  % fftshift creates a symmetric spectrum from -fs/2 to fs/2 instead of 0 to fs

figure
imagesc(vAxis, RAxis, log10(abs(SRx)))
ylabel('Range (m)')
xlabel('Speed (m/s)')