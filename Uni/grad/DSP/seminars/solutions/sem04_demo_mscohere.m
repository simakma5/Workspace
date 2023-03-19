% Demo - Odhad zpozdeni - zcorr a CPSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulovana data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=loadbin('SA001S01.CS0');
fs=16000;

% Filtr simulujici signal v druhem kanalu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delay = 20;

% idealni zpozdeni
%--------------------
M=delay; 

sim1_1 = s(M+1:end) ;
sim1_2 = s(1:end-M) ;

slen = length(s) ;

%----------------
% korelovany signal 
%-------------------
M=30;
% M=2*delay; 
% % sirokopasmova DP
% b=fir1(M,0.8) ;
% % sirokopasmova HP
b=fir1(M,0.05,'high') ;
% % sirokopasmova PP
% b=fir1(M,[0.05 0.8]) ;


sim2_1 = s ;
sim2_2 = filter(b,1,s) ;

%------------------------
% Pridani sumu
%------------------------
% G = 0 ;
G = 0.001 ;
sim2_1=sim2_1+G*randn(slen,1) ;
sim2_2=sim2_2+G*randn(slen,1) ;


figure(1)
subplot(211)
plot([sim1_1 sim1_2])
% xlim([1.8e4 2.1e4])
title('simulations: ideal delay')
subplot(212)
plot([sim2_1 sim2_2])
% xlim([1.8e4 2.1e4])
title('simulations: correlated with a delay (+noisy)')

wlen=512;

figure(2)
subplot(211)
mscohere(sim1_1,sim1_2,wlen,[],[],fs);
title('MSC via Welch: ideal delay')
subplot(212)
mscohere(sim2_1,sim2_2,wlen,[],[],fs);
title('MSC via Welch: correlated + noise')

%%%%%%%%%%%%%%%%%%%%%%%%%%%

flen=2000; 
wlen=512; 

% flen=512; 
% wlen=512; 

% flen=512; 
% wlen=128; 

figure(3)
subplot(311)
plot([sim2_1(2e4+1:2e4+flen) sim2_2(2e4+1:2e4+flen)])
title('signal - voiced frame')
subplot(312)
plot([sim2_1(3.8e4+1:3.8e4+flen) sim2_2(3.8e4+1:3.8e4+flen)])
title('signal - unvoiced frame')
subplot(313)
plot([sim2_1(1:flen) sim2_2(1:flen)])
title('signal - noise frame')

figure(4)
subplot(311)
mscohere(sim2_1(2e4+1:2e4+flen), sim2_2(2e4+1:2e4+flen),wlen)
title(['MSC - voiced frame: slen=' num2str(flen) ',wlen=' num2str(wlen)])
subplot(312)
mscohere(sim2_1(3.8e4+1:3.8e4+flen), sim2_2(3.8e4+1:3.8e4+flen),wlen)
title(['MSC - unvoiced frame: slen=' num2str(flen) ',wlen=' num2str(wlen)])
subplot(313)
mscohere(sim2_1(1:flen), sim2_2(1:flen),wlen)
title(['MSC - noise frame: slen=' num2str(flen) ',wlen=' num2str(wlen)])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  C O H E R E G R A M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flen=2048;
fstep=flen/2;
wlen=flen/4;
wstep=wlen/2;
NFFT=wlen;
fs=16000;

% Simulation I
fs=16000 ;
x=sim2_1 ;
y=sim2_2 ;

% Real signal I
% fs=16000;
% x=loadbin([sigdir 'SA001S01.CS0']);  % Head-set microphone
% y=loadbin([sigdir 'SA001S01.CS1']);  % Nokia hands-free microphone
% y=loadbin([sigdir 'SA001S01.CS2']);  % Middle-distance (1m) microphone
% y=loadbin([sigdir 'SA001S01.CS3']);  % Far-distance (3m) microphone

% Real signal II
% fs=16000 ;
% x=loadbin([sigdir 'SA176S01.CS0']);  % Head-set microphone
% y=loadbin([sigdir 'SA176S01.CS1']);  % Nokia hands-free microphone
% y=loadbin([sigdir 'SA176S01.CS2']);  % Middle-distance (1m) microphone
% y=loadbin([sigdir 'SA176S01.CS3']);  % Far-distance (3m) microphone

% Real signal III
% fs=16000;
% x=loadbin([sigdir 'SA005S01.CS0']);  % Head-set microphone
% y=loadbin([sigdir 'SA005S01.CS1']);  % Nokia hands-free microphone
% y=loadbin([sigdir 'SA005S01.CS2']);  % Middle-distance (1m) microphone
% y=loadbin([sigdir 'SA005S01.CS3']);  % Far-distance (3m) microphone

% Real signal IV - CAR DATA (cleaner) - symmetric microphone arrangement
% fs=8000 ;
% x=loadbin([sigdir 's0001-l.bin']);  %  1st (symmetric) channel
% y=loadbin([sigdir 's0001-r.bin']);  %  2nd (symmetric) channel

% Real signal V - CAR DATA (strong noise) - symmetric microphone arrangement
% fs=8000 ;
% x=loadbin([sigdir 'x0002-l.bin']);  %  1st (symmetric) channel
% y=loadbin([sigdir 'x0002-r.bin']);  %  2nd (symmetric) channel

figure(5)
subplot(311)
spectrogram(x,wlen,[],[],fs,'yaxis'); colorbar off; 
title('Spectrogram of 1st channel signal')
subplot(312)
spectrogram(y,wlen,[],[],fs,'yaxis'); colorbar off; 
title('Spectrogram of 2nd channel signal')
subplot(313)
cohegram(x,y,flen,fstep,wlen,wstep,NFFT,fs);
title('Cohegram')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Csim2,vtime,freq]=cohegram(x,y,flen,fstep,wlen,wstep,NFFT,fs);

figure(6)
subplot(311)
spectrogram(x,wlen,[],[],fs,'yaxis'); colorbar off; 
title('Spectrogram of 1st channel signal')
subplot(312)
cohegram(x,y,flen,fstep,wlen,wstep,NFFT,fs);
title('Cohegram')
subplot(313)
plot(vtime,mean(Csim2))
title('Short-time average coherence')

figure(7)
subplot(311)
cohegram(x,y,flen,fstep,wlen,wstep,NFFT,fs);
title('Cohegram')
subplot(312)
plot(vtime,mean(Csim2))
title('Short-time average coherence (full-band)')
subplot(313)
flow=1000;
fhigh=5000;
klow=floor(flow/fs*NFFT);
khigh=floor(min([fs/2 fhigh])/fs*NFFT);
plot(vtime,mean(Csim2(klow:khigh,:)))
title(['Short-time average coherence (band ' num2str(flow) '-' num2str(fhigh) 'Hz)'])

return

