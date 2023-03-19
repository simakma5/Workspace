%% TEST 1
%% LPC SPEKTRUM
a = [1.00000 -0.65061 0.27425 0.15223 0.32457 -0.23301 0.22956];
Ep = 1.322e-06;

H=freqz(sqrt(Ep),a,64);
LPCspec=10*log10((abs(H).^2));

%% vyhlazený odhad AR PROCES!!!!!!!!!

fs1=16000;
a = [1.00000 -2.54389 3.09003 -2.98049 2.61336 -1.54185 0.43826]; % u MA se prohoda A a B
G = 0.0023968;
b = 1;

un1 = loadbin('un-excitation.bin');

Num=10000
un=G*un1;

% figure(1)    %frekvencni charakteristika systemu
% freqz(b,a,1000,fs1)
% subplot(211) %vytvari osy

nc=filter(b,a,un);

wlen=512;

N = max([wlen 256]);
%vyhlazeny odhad PSD
hodnota = pwelch(nc,hamming(wlen),wlen/2,N,fs1,'twosided','psd');  %signal, okno, prekryv, ?????, frekvence vz.
vysledek = 10*log10(hodnota);
% title('PSD in dB using Welch method');

% (1:wlen).*hamming(wlen)

%% Vykon spektrum v DB

frame = loadbin('frame-001.bin');

wlen=length(frame) ;
w=hamming(wlen) ;
frame=frame.*w ;

FRAME = fft(frame) ;

DFTspec = 10*log10((abs(FRAME(1:8)).^2)./wlen) ;

%% Autocorel koef
sig = loadbin('sigar-002.bin');

p = 12;

wlen=length(sig) ;
w=hamming(wlen) ;
sig=sig.*w ;

a_lpc=lpc(sig,p);




























