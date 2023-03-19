%% init
clear
close all
clc

nc1 = loadbin("nc1.bin");  % LF colored noise
nc2 = loadbin("nc2.bin");  % HW colored noise
nc3 = loadbin("nc3.bin");  % banded colored noise
nc4 = loadbin("nc4.bin");  % unknown signal

N = 1e3;
fs = 8e3;
wlen = 512;

%% determination of the unknown signal
figure(1)
pwelch(nc4,wlen,wlen/2,fs,"onesided")
% the plotted graph has two peaks -> probably a 4th order filter

%% LPC analysis of an unknown noise
x = nc1; p = 1;
% x = nc2; p = 1;
% x = nc3; p = 2;
% x = nc4l p = 4;

% find the parameters using LP analysis
a = lpc(x,p);

% inspect the frequency response of such a filter
figure(2)
freqz(1,a,N,'whole')

% inspect the poles and zeros in the z-plane
figure(3)
zplane(a,p)

%% PSD estimate using LPC vs using DFT
% init
fs = 16e3;
x = loadbin('vm0-512.bin');
x = x.*hamming(512);
N = length(x);
p = 30;

% get the coefficients and the variance (power?) of the prediction error
[a,Ep] = lpc(x,p);

% from theory,  it should hold that |X(k)|^2/N \approx |H(k)|^2
% where the LHS is PSD estimate using DTFT and RHS using the LPC spectrum
Gx = abs(fft(x).^2)./N;
H = freqz(sqrt(Ep),a,N/2);
LPCspec = abs(H).^2;

% several commmon visual errors which might occur in the following
% a) the axes are mismatched (scaling is off) - watch out for the number
%   of points in which I calculate each function
% b) vertical shift - maybe forgotten scaling by energy
% c) different dynamic of each graph - have to plot squared module of both
%   -> has to be 10*log10(abs(.)^2) or 20*log10(abs(.))
figure(4)   
plot(10*log10(Gx(1:N/2)))
hold on
plot(10*log10(LPCspec),'r')
hold off
