close all; clear; clc

%% 1
cp=20;
wlen=512 ;
sig1 = loadbin("SA012S01.CSX");
sig2 = loadbin("SA012S01.CS0");


cepr1 = vrceps(sig1, 1, cp, wlen, wlen/2); %segmented cepstrum samples (matrix: rows -- cepstrums for each SEGMENTS, columns -- cepstrum indexes)
cepr2 = vrceps(sig2, 1, cp, wlen, wlen/2); %segmented cepstrum samples (matrix: rows -- cepstrums for each SEGMENTS, columns -- cepstrum indexes)


%cepstral distance 
cdist = cde(cepr2(:,2:end), cepr1(:,2:end), cp); %%%%%%%%%%%%%%%%%%%%%%%%%nejde????
%(+compute mean )
meancd = mean(cdist);

cdist/sum(abs(cepr1), 'all');

%% 2
sig1 = loadbin('frame-001.bin');
sig2 = loadbin('frame-012.bin');

sig1 = sig1.*hamming(512);
sig2 = sig2.*hamming(512);

ceps1 = vrceps(sig1, 1, 12, 512);
ceps2 = vrceps(sig2, 1, 12, 512);

cp = 12;
d = cde(ceps1, ceps2, cp);
dist = mean(d);
disp(dist)

%% 3
s = loadbin("SA018S01.CS0");
n = loadbin("nc3.bin");

scale = 0.8;
n = scale*n;
% x = s+scale*n;

P_s = mean(s.^2);
P_n = mean(n.^2);
snr = 10*log10(P_s/P_n);
disp(snr)

%% 5
load sigs_2chan_07.mat;
samp = 16*8;
w = hamming(samp);

msc_val = mscohere(sig1, sig2, w, samp/2, samp, 16e3);
msc_mean = mean(msc_val);
disp(msc_mean)

%% 6
x = loadbin("frame-002.bin");
X = dct(x);

y = idct(X(1:60), length(X));

figure(1);
subplot(211);
plot(x);
title('x(t) - před kompresí');
subplot(212);
plot(y);
title('y(t) - po kompresi');

figure(2);
subplot(211);
plot(X);
title('X(f) - před kompresí');
subplot(212);
plot(X(1:65));
title('Y(f) - po kompresi');

Px = mean(x.^2);
Py = mean(y.^2);

pomer = (Py/Px)*100;
disp(pomer)

%% 7
load sig_xy_01.mat;

wlen = 1024;
fs = 16e3;  % načteno z sig_xy_02.mat

[Pxy, FF] = cpsd(x, y, hamming(wlen), wlen/2, wlen, fs);

Pxydb = 10*log10(abs(Pxy));

figure(1);
plot(FF, Pxydb);

%% 8
x = loadbin('frame-023.bin');
w = hamming(length(x));
x = x.*w;

X = dctxc1(x);

X(1:10)