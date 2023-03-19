
format long;
clear, clc, close all;



%%%%%%%%%%%%%%%%%%%
% užijte DCT ke kompresi signálu, zachovejte prvních 55 složek
% určete jaké procento původního výkonu signálu je zahrnuto v signálu
% komprimovaném
% x = loadbin('frame-001.bin');
% 
% X = dct(x);
% 
% figure(1);
% subplot(211);
% plot(x);
% title('x(t) - před kompresí');
% subplot(212);
% plot(X);
% title('X(f) - před kompresí')
% 
% y = idct(Y, length(X));
% 
% figure(2);
% subplot(211);
% plot(y);
% title('y(t) - po kompresi')
% subplot(212);
% plot(Y);
% title('Y(f) - po kompresi');
% 
% Px = mean(x.^2);
% Py = mean(y.^2);
% 
% pomer = (Py/Px)*100




%%%%%%%%%%%%%%%%%%%
% určete reálné kepstrum signálu, váhujte hamming oknem příslušné délky a
% zobrazte celé vypočítané kepstrum.
% x = loadbin('frame-001.bin');
% 
% w = hamming(length(x));
% x = x.*w;
% 
% [XHAT, YHAT] = rceps(x);
% 
% figure(1);
% plot(XHAT);





%%%%%%%%%%%%%%%%%%%
% určete kosinovou transf. DCT-1 signálu. Před výpočtem váhujte hamming
% oknem, uveďte prvních 8 koef DCT spektra
% x = loadbin('frame-000.bin');
% w = hamming(length(x));
% x = x.*w;
% Xdct1 = dctxc1(x);
% 
% Xdct1(1:9)





%%%%%%%%%%%%%%%%%%%
% spočtěte vyhlazený odhad CPSD Welchovou metodou pro sig. x a y uložené v
% [###.mat] (load ###.mat). fs=16kHz. wlen=512, hamming, 50% overlap, fft
% stejný jako wlen, implicitní jednostranný odhad CPSD reálných signálů
% load sig_xy_04.mat;
% fs = 16000;
% 
% wlen = 512;
% 
% [Pxy, FF] = cpsd(x, y, hamming(wlen), wlen/2, wlen, fs);
% 
% phase = unwrap(angle(Pxy));
% 
% figure(1);
% plot(FF, phase);






%%%%%%%%%%%%%%%%%%%
% určete kepstrální vzdálenost na bází LPC kepstra mezi dvěma signály.
% Počítejte LPC kepstrum, řád p=16, váhujte hamming oknem. Vzdálenost
% počítejte z prvních 13 koef. včetně nultého, tj. z koef. c[0]-c[12].
% 
%%% chybí soubor cde.m, nelze použít pro výpočet vzdálenosti :c






%%%%%%%%%%%%%%%%%%%
% pro sig. sig1 a sig2 s fs=8kHz vypočtěte koherenční fci, konkrétně MSC
% (Magnitude Square Coherence), volte parametry wlen=32ms, hamming okno,
% overlap 50%, NFFT=wlen. Určete průměrnou koherenci na 3 platné cifry
% load sigs_2chan_08.mat;
% fs = 8000;
% 
% wlen = 0.064*fs;
% 
% [Cxy, FF] = mscohere(sig1, sig2, hamming(wlen), wlen/2, wlen);
% 
% figure(1);
% plot(FF, Cxy);
% 
% mean(Cxy)






%%%%%%%%%%%%%%%%%%%
% vypočítejte bázi KLT transformace signálu, jednotlivé realizace jsou v
% řádcích dané matice. První bázový vektor volte pro nejvýznamější
% komponentu. Určete prvních 10 komponent KLT spektra pro 1. realizaci
% daného signálu, tj. signálu v 1. řádku matice realizací sigframes





%%%%%%%%%%%%%%%%%%%
% SNR ze zašumělého a čistého signálu
% noisy = loadbin('SX012S01.CS0');
% clean = loadbin('SA012S01.CS0');
% 
% noise = noisy - clean;
% 
% Pclean = mean(clean.^2);
% Pnoise = mean(noise.^2);
% 
% SNRdb = 10*log10(Pclean/Pnoise)
% 
% 
% figure(1);
% subplot(311);
% plot(noisy);
% title('noisy signal');
% subplot(312);
% plot(clean);
% title('clean signal');
% subplot(313);
% plot(noise);
% title('noise');



