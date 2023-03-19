
%%%%%%%%%%%%%
% druhý test, 2022-12-14 9:15
% 
format long;
clear, clc, close all;



%%%%%%%%%%%%%%%%%%%%%%%
% úloha 1:  [OK]
% určete DCT-2 signálu, před výpočtem váhujte hamming oknem. Uveďte prvních
% koeficientů DCT spektra
% 
% x = loadbin('frame-023.bin');
% w = hamming(length(x));
% x = x.*w;
% 
% X = dctxc2(x);
% 
% X(1:10)



%%%%%%%%%%%%%%%%%%%%%%%
% úloha 2:  [OK]
% komprimujte signál pomocí DCT (dct a idct v matlab). Pro kompresi
% použijte prvních 65 komponent DCT spektra. Spočtěte výkony a zjistěte,
% jaké % výkonu původního signálu je obsaženo v sig komprimovaném
% x = loadbin('frame-013.bin');
% 
% X = dct(x);
% 
% y = idct(X(1:65), length(X));
% 
% figure(1);
% subplot(211);
% plot(x);
% title('x(t) - před kompresí');
% subplot(212);
% plot(y);
% title('y(t) - po kompresi');
% 
% figure(2);
% subplot(211);
% plot(X);
% title('X(f) - před kompresí');
% subplot(212);
% plot(X(1:65));
% title('Y(f) - po kompresi');
% 
% Px = mean(x.^2);
% Py = mean(y.^2);
% 
% pomer = (Py/Px)*100


%%%%%%%%%%%%%%%%%%%%%%%
% úloha 3:  [i have no fucking clue]
% určete LPC KEPSTRUM signálu. Řád LPC p=16, váhujte hamming oknem.
% Zobrazte prvních 21 koef včetně nultého, tj. c[0]-c[20]
% 
% x = loadbin('frame-002.bin');
% w = hamming(length(x));
% x = x.*w;
% 
% p = 16;
% A = lpc(x, p);
% 
% žádný výsledek :c no idea what to do






%%%%%%%%%%%%%%%%%%%%%%%
% úloha 4:   [nejspíš ok]
% spočítejte CPSD pro sig. x a y, fs = 16kHz. Volte parametry wlen=1024,
% hamming okno, 50% overlap, NFFT=wlen, implicitní jednostranný odhad
% reálných signálů. Vykreslete odhad modulu CPSD v decibelech!
% load sig_xy_02.mat;
% 
% wlen = 1024;
% % fs = 16kHz načteno z sig_xy_02.mat
% 
% [Pxy, FF] = cpsd(x, y, hamming(wlen), wlen/2, wlen, fs);
% 
% Pxydb = 10*log10(abs(Pxy));
% 
% figure(1);
% plot(FF, Pxydb);



%%%%%%%%%%%%%%%%%%%%%%%
% úloha 5:   [odhad z nejbližší hodnoty]
% určete euklid. kepstr. vzdálenost reálného kepstra mezi dvěma signály.
% počítejte reálné kepstrum, váhujte hamming oknem, vzdálenost počítejte z
% prvních 13 koef. včetně nultého, tj. c[0]-c[12]

x = loadbin('frame-001.bin');
y = loadbin('frame-012.bin');

w = hamming(length(x));
x = x.*w;
y = y.*w;

[X_ceps, YHAT_x] = rceps(x);
[Y_ceps, YHAT_y] = rceps(y);

X_ceps = X_ceps';
Y_ceps = Y_ceps';

D = cde(X_ceps, Y_ceps, 13);
disp(D)


%%%%%%%%%%%%%%%%%%%%%%%
% úloha 6:   [asi ok]
% pro sig1 a sig2, fs=16kHz vypočtěte koherenční funkci, konkrétně MSC,
% volte parametry wlen=8ms, hamming okno, 50% overlap, NFFT=wlen.
% Určete průměrnou koherenci
% load sigs_2chan_01.mat;
% 
% wlen = 0.008*fs;
% 
% [Cxy, FF] = mscohere(sig1, sig2, hamming(wlen), wlen/2, wlen, fs);
% 
% mean(Cxy)




%%%%%%%%%%%%%%%%%%%%%%%
% úloha 7:   [OK]
% odstup SNR
% 
% 
% noisy = loadbin('SX019S01.CS0');
% clean = loadbin('SA019S01.CS0');
% 
% noise = noisy - clean;
% 
% Pclean = mean(clean.^2);
% Pnoise = mean(noise.^2);
% 
% SNRdb = 10*log10(Pclean/Pnoise);
% SNRdb



%%%%%%%%%%%%%%%%%%%%%%%
% úloha 8:   [i have no fucking clue]
% určete zkreslení delšího sig. na bázi kepstrální vzdálenosti a LPC
% kepstra. Počítejte LPC kepstrum pro wlen=512, 50% overlap a uvažujte
% implicitní hamming váhování. Řád LPC p=16, počet kepstrálních koef (bez
% c[0]) volte cp=20 a vzdálenost počítejte na bázi euklid. vzdál. včetně
% nultého koeficientu, tj. z koef. c[0]-c[20]
% 
% zkresleny = loadbin('SA012S01.CSX');
% referencni = loadbin('SA012S01.CS0');
% 
% i have no fucking idea what to do :c






%%%%%%%%%%%%%%%%%%%%%%%
% úloha :   []
% 
% 
% 









%%%%%%%%%%%%%%%%%%%%%%%
% úloha :   []
% 
% 
% 







