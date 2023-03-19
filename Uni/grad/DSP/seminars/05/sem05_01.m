%% init

N = 100;
fs = 8e3;
f = 5e2;
Tau = 0.001;

t = ((0:N-1)/fs)';
s = cos(2*pi*f*t).*exp(-t/Tau);

figure(1)
subplot(311);
plot(s)
title('Casovy prubeh')

subplot(312);
m = 0.7;
M = 15;
% H(z) = 1 + m*z^(-M);
h = zeros(N,1);
h(1) = 1;
h(M+1) = m;
stem(h)
title('Impulsova odezva')

subplot(313);
x = filter(h,1,s);
plot(x)
title('Vysledny signal s odrazem')

%% spektra
figure(2)
subplot(311);
S = abs(fft(s));
plot(20*log10(S))
title('Spektrum casoveho prubehu')

subplot(312);
H = abs(fft(h));
plot(20*log10(H))
title('Spektrum impulsove odezvy')

subplot(313);
X = abs(fft(x));
plot(20*log10(X))
title('Spektrum vysledneho signalu s odrazem')

%% kepstrum
% kepstra z definice
figure(3)
subplot(311);
Cs = ifft((log(abs(fft(x,N)))),N);
plot(Cs)
title('Kepstrum casoveho prubehu')

subplot(312);
Ch = ifft((log(abs(fft(h,N)))),N);
plot(Ch)
title('Kepstrum impulsove odezvy')

subplot(313)
Cx = ifft((log(abs(fft(x,N)))),N);
plot(Cx)
title('Kepstrum vysledneho signalu s odrazem')

%% odstraneni odrazu (konvolucniho sumu)
L = 12;  % chceme zachovat (odhadem) prvnich 12 vzorku - krome nulteho
% lmask = [ones(L+1,1),zeros(N-L-1,1)]; by bylo spatne, protoze kepstrum
% musi byt vzdycky symetricke
lmask = [ones(L+1,1); zeros(N-2*L-1,1); ones(L,1)];  % na konci uz nulty clen neni, takze jen L

figure(3)
hold on
subplot(313);
Cxlift = Cx.*lmask;
plot(Cxlift)
hold off

%% prevod liftrovaneho kepstra zpet
figure(2)
Xlift = exp(fft(Cxlift));
% plot(Xlift)  % bordel: vlivem aproximaci je tam malinka imag. slozka
Xlift = real(Xlift);  % korekce

% tady zbytkem uz si nejsem jistej
subplot(311);
hold on
plot(20*log10(Xlift))
hold off

subplot(313);
hold on
plot(20*log10(Xlift))
hold off

%% realne kepstrum recoveho signalu (hodne zrychlene)
clear, clc, close all;

N = 512;
s = loadbin('vm0.bin');
s = s(end-N+1:end);
s = s.*hamming(N);

S = abs(fft(s));

cs = rceps(s);  % real cepstrum

figure(5)
subplot(311);
plot(s)

subplot(312)
plot(20*log10(S))

subplot(313);
plot(cs)

L = 16;
lmask = [ones(L+1,1); zeros(N-2*L-1,1); ones(L,1)];
cslift = cs.*lmask;

Slift = real(exp(fft(cslift)));

figure(5)
subplot(313);
hold on
plot(cslift)
hold off

subplot(312)
hold on
plot(20*log10(Slift))
hold off
