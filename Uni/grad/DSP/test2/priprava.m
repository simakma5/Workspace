%% Nevyhlazené výkonové spektrální hustoty
Sxx = (1/NFFT)*(abs(fft(x))).^2;
Syy = (1/NFFT)*(abs(fft(y))).^2;

Sxy= (1/NFFT)*abs(fft(x).*conj(fft(y)));

%% Vyhlazené výkonové spektrální hustoty

% Sxx = PSD_x
PSD_x = pwelch(x,wlen,wstep,NFFT);
PSD_y = pwelch(y,wlen,wstep,NFFT);
% Sxy nezískám pwelchem

%% MSC
%MSC = mscohere(frame1,frame2,wlen2,Noverlap2,NFFT,fs);

%%
% -------------------------------------------------------------------------
% TESTY:
% -------------------------------------------------------------------------

%% 1) průměrná koherence
load("sigs_2chan_06.mat")

fs = 8000;
wdur = 0.016;
wlen = wdur*fs;
wstep = 0.75*wlen;

[Cxy,F] = mscohere(sig1, sig2,wlen,wstep,128);

av = sum(Cxy(:))/65
sprintf('%f', mean(Cxy))
mean(Cxy)

%% 2) vyhlazený odhad vzájemné spektrální hustoty
load("sig_xy_02.mat")

fs = 16000;

wlen = 1024;        % délka krátkodobého signálu
wstep = 0.5*wlen;   % 50% překryv

[Sxy,FF]=cpsd(x,y,wlen,[],wstep,fs);

figure(1)
subplot(211)
plot(FF,10*log10(abs(Sxy)));    % modul CPSD v dB
subplot(212)
plot(FF,unwrap(angle(Sxy)));    % fáze CPSD v rad


%% 3) Prvních 10 hodnot koherenční funkce
load("sigs_2chan_08.mat")

fs = 8000;
wdur = 0.642;   % délka krátkodobého segmentu, moc dlouhá + nesmyslný číslo
wlen = wdur*fs; 

%wlen = 512; %s tímhle to vychází dobře -> wdur = 0.064
[Cxy,F] = mscohere(sig1,sig2,wlen,0.5*wlen);
Cxy(1:10)

%% 4) Určete kosinou transformaci DCT1 nebo DCT2 (vzorec)

x = loadbin('');
x = x*hamming(length(x));
XC2 = dctxc2(x);        % Pollákova funkce
XC2(1:8)
% XC1 = dctxc1(x);
% XC1(1:8)

%% 5) Komprimujte signál

signal = loadbin("frame-006.bin");
keep_N = 80;

signal_dct = dct(signal);
signal_dct(keep_N+1:end) = 0;
signal_compr = idct(signal_dct);

figure(1);
subplot(2, 1, 1);
plot(signal);
subplot(2, 1, 2);
plot(signal_compr);

P_signal = mean(signal.^2);
P_signal_compr = mean(signal_compr.^2);

fprintf("Zastoupeni: %f %%\n", P_signal_compr/P_signal*100);

%% 6)Zkresleni
clc, clear, close all
sig1 = loadbin('SA015S01.CSX');
sig2 = loadbin('SA015S01.CS0');

w = hamming(length(sig1));
cr_sig1 = (ifft(log(abs(fft(sig1.w)))));
cr_sig2 = (ifft(log(abs(fft(sig2.w)))));


cd0_cr_sig=cde(cr_sig1(2:20),cr_sig2(2:20));

%% 7) Euklidovská kepstrální vzdálenost na bázi LPC

sig1 = loadbin('SA001S01.CS0');
sig2 = loadbin('SA001S01.CS1');

cp = 12;            % počet kepstrálních koeficientů
wlen = 512;         % délka[ms]*fs
wstep = wlen/2;
p = 16;             % řád modelu

%LPC ceps.:
ceps1 = vaceps(sig1,1,p,cp,wlen,wstep);
ceps2 = vaceps(sig2,1,p,cp,wlen,wstep);

% Vzdálenost mezi kepstry, bez c[0]
d= cde(ceps1(2:end),ceps2(2:end));      % (2:end) protože nechci výkon (c[0])

figure
subplot(211)
plot(ceps1)
grid on, title('LPC ceps sig1')
subplot(212)
plot(ceps2)
grid on, title('LPC ceps sig2')

%% 8) Euklidovská keps. vzdálenost na bázi real. ceps

sig1 = loadbin('SA001S01.CS0');
sig2 = loadbin('SA001S01.CS1');

cp = 13;            % počet kapstrálních koeficientů
wlen = 512;         % délka[ms]*fs
wstep = wlen/2;

cr_sig1 = vrceps(sig1,1,cp,wlen,wstep);   %implicitně hammingovo okno
cr_sig2 = vrceps(sig2,1,cp,wlen,wstep);

% Vzdálenost mezi kepstry, s c[0]
d= cde(ceps1,ceps2);      % chci i c[0]

%% 9) Mix signálů
% WTF?

% SO_FFT_sigs_01.mat - spektra v řádcích

slen = length(n);
wlen = 512 ;
wstep = wlen/2 ;
wnum = (slen-wlen)/wstep+1 ;
w = hamming(wlen);
N = [];

for i=1:wnum
   ii = (i-1)wstep+1;
   jj = (i-1)wstep+wlen;

   framen=n(ii:jj).*w;
   yframen = abs(fft(framen));

   N = [N yframen];
end
N = mean(N,2);

%odhad toho sumu pred ola a odecitanim
sig = x;
slen=length(sig);
wlen=512;
wstep=wlen/2 ;
wnum=(slen-wlen)/wstep+1 ;
w=hamming(wlen);

out = zeros(slen,1);
for i=1:wnum
   ii=(i-1)wstep+1;
   jj=(i-1)wstep+wlen;

   frame = fft(sig(ii:jj).w);
   % jednocestne usmerneni
   frameAbs = abs(frame) - N;
%  frameAbs = (frameAbs + abs(frameAbs))/2; tato char _/
%  frameAbs(frameAbs < 0) = 0; % stejny vysledek

   rameAbs = abs(frameAbs); % oboustran. usm

   frame = frameAbs.exp(1i.*angle(frame));

   yframe = ifft(frame);
   out(ii:jj)=out(ii:jj)+yframe ;
end

scale2 = sum(w(1:wstep:wlen)) ;
out = real(out./scale2 );
figure(4)
subplot(211)
plot(s)
title('Clean signal')
subplot(212)
plot(out)
title('Cleaned signal with spectral substraction')

%% 10
% noisy = loadbin("SX011S01.CS0");
% clean = loadbin("SA011S01.CS0");
% % SNR 
% noise = noisy - clean;
% 
% Pnoise = mean(abs(noise).^2);
% Psig = mean(abs(clean).^2);
% 
% SNR_out = 10*log10(Psig/Pnoise) % 3dB