%% 1
a = [ 1.00000 -0.25622 0.47595 0.03518 0.03846 0.01822 -0.20387 ];
Ep = 1.2384e-06;

H=freqz(sqrt(Ep),a,64);
LPCspec=10*log10((abs(H).^2));

%% 2
sig = loadbin('sigar-000.bin');
sig = sig.*hamming(length(sig));
p = 12;

a_lpc=lpc(sig,p);

%% 3
fs = 16e3;
a = [ 1.00000 -0.65061 0.27425 0.15223 0.32457 -0.23301 0.22956 ]; % u MA se prohoda A a B
G = 0.0011498;
b = 1;

un1 = loadbin('un-excitation.bin');

un=G*un1;

% figure(1)
% freqz(b,a,1000,fs)

nc=filter(b,a,un);

wlen=512;

N = max([wlen 256]);
welch = pwelch(nc,hamming(wlen),wlen/2,N,fs,'twosided','psd');  %signal, okno, prekryv, ?????, frekvence vz.
welch_dB = 10*log10(welch);

%% 4
x = loadbin('frame-002.bin');
wlen=length(x) ;
x=x.*hamming(wlen) ;

X = fft(x) ;

DFTspec = 10*log10((abs(X(1:8)).^2)./wlen) ;
