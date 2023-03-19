% LPC spectral analysis - sinusoids in noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Analyzed-signal generation


flen=4000 ;

fs=8000 ;

f1=878 ;
A1=0.8 ;
s1=sin(2*pi*f1*(0:flen-1)./fs) ;

f2=2321 ;
A2=0.7 ;
s2=sin(2*pi*f2*(0:flen-1)./fs) ;

noise=0.5*randn(1,flen) ;

frame=s1 ;
frame=s1 + noise ;
% frame=s1 + s2 ;
frame=s1 + s2 + noise;

% Weighting  -  needed for speech analysis
%   (for some other applications can be omitted)

w=hamming(flen) ;
frame=frame(:) ;   % Vector should be a column one
% frame=frame.*w ;



% 7) Power DFT spectrum and LPC spectrum (PSD) in dBs of SAME LENGTH
FRAME = fft(frame) ;
DFTspec = 10*log10((abs(FRAME).^2)./flen) ;

% A choise of the LPC model order
% p=2 ;  
% p=4 ;
p=8 ;


[a,Ep]=lpc(frame,p) ;
H=freqz(sqrt(Ep),a,flen/2) ;
LPCspec=10*log10((abs(H).^2)) ;

figure(1)
plot(DFTspec(1:flen/2));
hold on
plot(LPCspec,'r') ;
hold off
title('DFT- and LPC-based power spectrum (PSD) in dBs')
ylim([-40 40])


