%% MODELOVÁNÍ A PSD
close all;
clear all;
% Šířka pásma = B
% Vzorkovací kmitočet = fs
% Nulový bod v ...
% Násobná konstanta rovna jedné
% Vybudit šumem un = loadbin('un-excitation.bin')

% Vyhlazený odhad spektrální výkonové hustoty (PSD) Welchovou metodov v dB - twosided
% Délka segmentu a počet bodů FFT = wlen
% Segmentace s 50% překryvem
% Implicitní váhování Hammingovým oknem

B = 500;
fs = 1600;
Ts = 1/fs;
wlen = 512;
un = loadbin('un-excitation.bin')
PSD = PSD_ARMA_first(un,B,fs,wlen,'AR','NF');

%% Výkon chyby predikce
close all;
clear all;

sig=loadbin('nazev.bin');
sig = sig(:);
p = 12;

figure(1)
plot(sig);

w = hamming(512);

frame_w = sig.*w;
[x Ep] = lpc(frame_w, p);

Ep %-> výsledek 3.3358e-07

%% LPC spektrum AR
close all;
clear all;

a = [ 1.00000 -0.65061 0.27425 0.15223 0.32457 -0.23301 0.22956 ]

E = 1.322e-06 %Výkon chyby predikce

%LPC spektrum v X bodech

spec = freqz(sqrt(E),a,X);
spec_dB = 10*log10(abs(spec).^2);
disp(spec_dB(1:10))

%-> b)

%% Výkonové spektrum v dB
close all;
clear all;

sig=loadbin('nazev.bin');

figure(2)
plot(sig)

w = hamming(512);
sig_w = sig.*w;

SPD = (abs(fft(sig_w))).^2/512;
SPD_dB = 10*log10(SPD);

SPD_dB(1:8) %-> a)

%psdnc=10*log10(((abs(fft(sig.*hamming(length(sig))))).^2)/length(sig));
%psdnc=psdnc(1:8); % -> a)





