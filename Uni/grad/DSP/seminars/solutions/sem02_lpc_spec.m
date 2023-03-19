% LPC spectral analysis - sample of voiced speech
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Analyzed signal

frame = loadbin ( 'vm0-512.bin' ) ;

% frame = loadbin ( 'vm0.bin' ) ;
% frame=frame(1201:1712);
% savebin('../signaly/vm0-512.bin',frame);

% Weighting  -  needed for speech analysis
%   (for some other applications can be omitted)

wlen=length(frame) ;
w=hamming(wlen) ;
frame=frame.*w ;

FRAME = fft(frame) ;

figure(1) ;
plot(frame) ;


% 2) Basic computation of DFT spectrum

figure(2) ;
plot(abs(FRAME)) ;
title('Magnitude DFT spectrum')

% 3) Basic computation of AR model (LPC)
%    Frequency response with numerator of transfer function equal to 1

p=16;

[a, Ep]=lpc(frame,p);

figure(3)
N=1000;
freqz(1,a,1000) ;
title('Frequency response of AR modelu')

% 4) Magnitude DFT spectrum in dBs

DFTspec = 10*log10(abs(FRAME)) ;
figure(4) ;
plot(DFTspec) ;
title('Magnitude DFT spectrum in dBs')

H=freqz(1,a,N) ;
LPCspec=10*log10(abs(H)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particular steps of comparison of DFT- and LPC-based spectrum
%
% 5) DFT and LPC spectrum in dBs of different lengths

figure(5)
plot(DFTspec(1:wlen/2));
hold on
plot(LPCspec,'r') ;
hold off
title('DFT and LPC spectrum in dBs - DIFFERENT LENGTHS')


% 6) DFT and LPC spectrum in dBs of same lengths - bad scaling

H=freqz(1,a,wlen/2) ;
LPCspec=10*log10(abs(H)) ;

figure(6)
plot(DFTspec(1:wlen/2));
hold on
plot(LPCspec,'r') ;
hold off
title('DFT and LPC spectrum in dBs of same lengths - BAD SCALING')


% 7) Power DFT and LPC spectrum in dBs of same lengths - CORRECT SCALIN

DFTspec = 10*log10((abs(FRAME).^2)./wlen) ;
H=freqz(sqrt(Ep),a,wlen/2) ;
LPCspec=10*log10((abs(H).^2)) ;

figure(7)
plot(DFTspec(1:wlen/2));
hold on
plot(LPCspec,'r') ;
hold off
title('Power DFT and LPC spectrum in dBs of same lengths - CORRECT SCALING')


% 8) POWER DFT and MAGNITUDE LPC spectrum in dBs - different dynamics

DFTspec = 10*log10((abs(FRAME).^2)./wlen) ;
H=freqz(sqrt(Ep),a,wlen/2) ;
LPCspec=10*log10((abs(H))) ;

figure(8)
plot(DFTspec(1:wlen/2));
hold on
plot(LPCspec,'r') ;
hold off
title('POWER DFT and MAGNITUDE LPC spectrum in dBs - different dynamics')


% 9) Power DFT and LPC spectrum in dBs of same lengths - VARIOUS ORDERS

% Choice of the order
% p=6;
% p=10;
% p=16;
% p=30;
p=100;

[a, Ep]=lpc(frame,p);


DFTspec = 10*log10((abs(FRAME).^2)./wlen) ;
H=freqz(sqrt(Ep),a,wlen/2) ;
LPCspec=10*log10((abs(H).^2)) ;

figure(9)
plot(DFTspec(1:wlen/2));
hold on
plot(LPCspec,'r') ;
hold off
title(['DFT and LPC spectrum - AR model ORDER p=' num2str(p)])


