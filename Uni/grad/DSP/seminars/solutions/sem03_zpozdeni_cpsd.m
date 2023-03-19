% Demo - Dealy estimation using CCF and CPSD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=loadbin('SA001S01.CS0');
fs=16000;

delay = 20;

% ideal delay
%--------------------
M=delay; 

sim1_1 = s(M+1:end) ;
sim1_2 = s(1:end-M) ;

slen = length(s) ;

%----------------
% More realistic model of delayed signal 
%    (i.e. correlated variant only with additive noise) 
%-------------------

% Filter for the creation of correlated signal in 2nd channel
%-------------------

M=30;
% M=2*delay; 
% % wide-band low-pass filter
% b=fir1(M,0.8) ;
% % wide-band high-pass filter
b=fir1(M,0.05,'high') ;
% % wide-band band-pass filter
% b=fir1(M,[0.05 0.8]) ;


sim2_1 = s ;
sim2_2 = filter(b,1,s) ;

%------------------------
% Addition ogPridani sumu
%------------------------
% G = 0 ;     % no noise
G = 0.005 ;    % scaling constant for amplification of background white noise
sim2_1=sim2_1+G*randn(slen,1) ;
sim2_2=sim2_2+G*randn(slen,1) ;


xmin=5.9e4 ;
xmax=6.2e4;


figure(1)
subplot(211)
plot([sim1_1 sim1_2])
xlim([xmin xmax])  % Comment to see the signal of full length
title('simulations: ideal delay')
subplot(212)
plot([sim2_1 sim2_2])  % Comment to see the signal of full length
xlim([xmin xmax])
title('simulations: correlated with a delay (+noisy)')

figure(2)
subplot(311)
spectrogram(sim1_1,512,'yaxis'); colorbar off; colormap jet
title('Spectrogram - SIMULATION 1 (sim1\_1)')
subplot(312)
spectrogram(sim2_1,512,'yaxis'); colorbar off; colormap jet
title('Spectrogram - SIMULATION 2 (sim2\_1)')
subplot(313)
spectrogram(sim2_2,512,'yaxis'); colorbar off; colormap jet
title('Spectrogram - SIMULATION 2 (sim2\_2)')

% Possible observation of frequency response of used filter
figure(3); freqz(b,1,1000)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real1_0=loadbin('SA001S01.CS0');
real1_1=loadbin('SA001S01.CS1');
real1_2=loadbin('SA001S01.CS2');
real1_3=loadbin('SA001S01.CS3');

figure(4)
subplot(311)
plot([real1_0 real1_1])
xlim([xmin xmax])
title('real signals:  CS0 - CS1')
subplot(312)
plot([real1_0 real1_2 ])
xlim([xmin xmax])
title('real signals:  CS0 - CS2')
subplot(313)
plot([real1_0 real1_3])
xlim([xmin xmax])
title('real signals:  CS0 - CS3')

figure(5)
subplot(221)
spectrogram(real1_0,512,'yaxis'); colorbar off; colormap jet
title('REAL SIGNALS (channel CS0)')
subplot(222)
spectrogram(real1_1,512,'yaxis'); colorbar off; colormap jet
title('REAL SIGNALS (channel CS1)')
subplot(223)
spectrogram(real1_2,512,'yaxis'); colorbar off; colormap jet
title('REAL SIGNALS (channel CS2)')
subplot(224)
spectrogram(real1_3,512,'yaxis'); colorbar off; colormap jet
title('REAL SIGNALS (channel CS3)')

%% Signal selection and delay computation


wlen=1024;

% x1=sim1_1(xmin+1:xmin+wlen) ;
% x2=sim1_2(xmin+1:xmin+wlen) ;
x1=sim2_1(xmin+1:xmin+wlen) ;
x2=sim2_2(xmin+1:xmin+wlen) ;
% x1=real1_0(xmin+1:xmin+wlen) ;
% x2=real1_1(xmin+1:xmin+wlen) ;
% x2=real1_2(xmin+1:xmin+wlen) ;
% x2=real1_3(xmin+1:xmin+wlen) ;



% 1) Method based on cross-correlation

figure(10);
subplot(211);
plot([x1 x2]);
title('Analyzed frames of delayed signals');

R12 = xcorr ( x1, x2, 'biased' ) ;


[maxR12,idx]=max(R12) ;

framedelay =  idx(1)-wlen ;
   
subplot(212)
kk=-wlen+1:wlen-1 ;
plot(kk,R12)
a=axis;
h=line([framedelay framedelay],[a(3) a(4)],'Color','r','LineStyle','--');
title(['CCF of selected analyzed frame: estimated delay = ' num2str(framedelay)]);
xlabel ('Delay')
grid

% 2) Method based on cross power spectral density

fs=16000;

% [Sxy,FF]=cpsd(x1,x2) ;

wlen2 = 128; 
noverlap=3*wlen2/4;
NFFT=512;
[Sxy,FF]=cpsd(x1,x2,wlen2,noverlap,NFFT,fs) ;
[Sxyrr,RR]=cpsd(x1,x2,wlen2,noverlap,NFFT) ;
% [Sxy,FF]=cpsd(x1,x2,wlen2,noverlap) ;  % Frequency is in radians 

NDFT=2*(length(Sxy)-1);

figure(11)
subplot(211)
plot(FF,20*log10(abs(Sxy)));
title('Magnitude CPSD of selected analyzed frame');
grid
subplot(212)
plot(FF,unwrap(angle(Sxy)))
title('Phase CPSD of selected analyzed frame');
xlabel('f - Frequency [Hz]')
ylabel('phase [rad]')
grid

figure(12)
subplot(211)
plot(RR,20*log10(abs(Sxy)));
title('Magnitude CPSD of selected analyzed frame');
grid
subplot(212)
plot(RR,unwrap(angle(Sxy)))
title('Phase CPSD of selected analyzed frame');
xlabel('\theta - Nyquist frequency [rad\cdot s^{-1}]')
ylabel('phase [rad]')
grid

figure(13)
subplot(211)
plot(20*log10(abs(Sxy)));
title('Magnitude CPSD of selected analyzed frame');
grid
subplot(212)
plot(unwrap(angle(Sxy)))
title('Phase CPSD of selected analyzed frame');
xlabel('k - discrete frequency [-]')
ylabel('phase [rad]')
grid

% %  POINTS (obtained manually for sim1_1 - sim1_2):   
% %      - Figure(11) - frequency in Hz
% A1 = [1562 12.4] ;
% A2 = [4500 35.06] ;
% %      - Figure(12) - Nyquist frequency in radians
% A1RR = [0.1595 3.319] ;
% A2RR = [1.755 34.98] ;
% %      - Figure(13) - Discrete frequency in samples
% A1kk = [19  4.185] ;
% A2kk = [141 34.32] ;

% POINTS (obtained manually for sim2_1 - sim2_2):   
%     - Figure(11) - frequency in Hz
A1 = [375 2.359] ;
A2 = [2906 17.35] ;
%      - Figure(12) - Nyquist frequency in radians
A1RR = [0.09817 1.511] ;
A2RR = [1.325 19.92] ;
%      - Figure(13) - Discrete frequency in samples
A1kk = [9 1.52] ;
A2kk = [107 19.58] ;

% for frequency realted to fs 
delay_Sxy =  (A2(2)-A1(2)) / (A2(1)-A1(1)) /2/pi * fs  ;

% for frequency defined in RADIANS
delay_SxyRR = (A2RR(2)-A1RR(2)) / (A2RR(1)-A1RR(1)) ;

% for discrete frequency defined in samples - k
delay_Sxykk = (A2kk(2)-A1kk(2)) / (A2kk(1)-A1kk(1)) /2/pi * NFFT  ;

figure(11)
subplot(212)
hold on
plot([A1(1) A2(1)],[A1(2) A2(2)],'rx')
hold off
title(['Phase CPSD of selected analyzed frame = ' num2str(delay_Sxy)  ' samples']);

figure(12)
subplot(212)
hold on
plot([A1RR(1) A2RR(1)],[A1RR(2) A2RR(2)],'rx')
hold off
title(['Phase CPSD of selected analyzed frame = ' num2str(delay_SxyRR)  ' samples']);

figure(13)
subplot(212)
hold on
plot([A1kk(1) A2kk(1)],[A1kk(2) A2kk(2)],'rx')
hold off
title(['Phase CPSD of selected analyzed frame = ' num2str(delay_Sxykk)  ' samples']);



%% Group delay

Sxy_phase=unwrap(angle(Sxy));
M=15; Sxy_phase=filter(ones(1,M+1)./(M+1),1,Sxy_phase);

Sxy_gd=diff(Sxy_phase) /2/pi * NFFT;

figure(14)
subplot(211)
plot(Sxy_phase)
title('Phase (smoothed) CPSD of selected analyzed frame');
grid
subplot(212)
plot(Sxy_gd)
title('Delay based on difference of (smoothed) CPSD phase');
grid




