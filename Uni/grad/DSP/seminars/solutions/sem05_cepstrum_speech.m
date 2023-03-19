% Demo cepstral analysis
%    SPEECH CEPSTRUM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading of speech signal

s=loadbin('vm0.bin');
% s=loadbin('uf5.bin');
% s=loadbin('vf5.bin');

wlen=512 ;
w=hamming(wlen);

% Selection of short-time frame

s = s (end-wlen+1:end).*w;      % From the end of signal (more stationary) 
% s = s (1:wlen).*w;              % From the beginning
% s = s (end-wlen+1:end);         % Without weighting 
% s=[s(211:310); zeros(412,1)] ;  % Selection of one period only and zero padding

N=length(s) ;

figure(10)
plot(s);
title('Speech frame');

pause(2)

Sabslog =log(abs(fft(s)));
cs=rceps(s) ;

figure(11)
plot(Sabslog);
title('Log spectrum of speech');

pause(2)

figure(12)
plot(cs);
title('Real cepstrum of speech');

pause(2)

% Lifter definition

nc = 16 ;
HL = [ones(nc+1,1); zeros(N-2*nc-1,1); ones(nc,1) ];
HH = [zeros(nc+1,1); ones(N-2*nc-1,1); zeros(nc,1) ];

% Liftering

cs1 = cs.*HL ;
cs2 = cs.*HH ;

% log-magnitude spectrum of liftered signal

SS1abslog = real((fft(cs1)));
SS2abslog = real((fft(cs2)));

% Observation of aperiodic part 
% (TYPICAL RESULT for speech analysis)

figure(12)
hold on
plot(cs1,'r');
hold off
title('Real cepstrum: blue - original, red - LK lifter') ;

pause(2)

figure(11)
hold on
plot(SS1abslog,'r');
hold off

title('Log spectrum:  blue - original, red - LK lifter') ;

pause(2)

% Observation of periodic part 
% (usually not used for speech analysis)

figure(13)
plot(cs);
hold on
plot(cs2,'g');
hold off
title('Real cepstrum: blue - original, green - HK lifter') ;

pause(2)

figure(14)
plot(Sabslog);
hold on
plot(SS2abslog,'g');
hold off

title('Log spectrum:  blue - original, green - HK lifter') ;





