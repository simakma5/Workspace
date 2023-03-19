% Demo cepstral analysis
%    REFLECTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of attenuated cosine
f=500 ;
fs=8000 ;
Tau=0.001 ;

% Generation of time and cosine (result should be COLUMN VECTORS !!!)
N=100 ;
t=(0:N-1)'./fs ;
s=cos(2*pi*f*t).*exp(-t/Tau) ;

% Parameters and impulse response of filter for reflection modeling
M=15; 
m=0.7 ;

h=zeros(N,1); 
h(1)=1;
h(M+1)=m;

% Creation of reflection

ss=filter(h,1,s);

% Particular signals in time domain

figure(1)
subplot(311)
plot(s);
title('Attenuated cosine')
subplot(312)
hh=impz(h,1,N);
% stem(hh);
plot(hh,'.-');
title('Impulse response' )
subplot(313)
plot(ss)
title('Attenuated cosine with a reflection')

pause(2)

% Particular signals in frequency domain

figure(2)
subplot(311)
Sabs =abs(fft(s));
plot(Sabs);
title('Spectrum of attenuated cosine ')
subplot(312)
Habs =abs(fft(hh));
plot(Habs);
title('Frequency repsonse of reflection filter')
subplot(313)
SSabs =abs(fft(ss));
plot(SSabs);
title('Spectrum of attenuated cosine with a reflection')

pause(2)

% Particular signals in cepstral domain

figure(3)
subplot(311)
cs=rceps(s) ;
% cs=real(ifft(log(abs(fft(s))))) ;
plot(cs);
title('Cepstrum of attenuated cosine ')
subplot(312)
ch=real(ifft(log(abs(fft(h))))) ;
plot(ch);
subplot(313)
css=rceps(ss) ;
plot(css);
title('Cepstrum of attenuated cosine with a reflection')

% return

pause(2)


% Lifter definition

nc = 12 ;    % Tray also nc = 14, 20, 10, 6

% Lifter for aperiodic component
HL = [ones(nc+1,1); zeros(N-2*nc-1,1); ones(nc,1) ];

% Lifter for periodic component
HL2 = [zeros(nc+1,1); ones(N-2*nc-1,1); zeros(nc,1) ];

% Liftering

css1 = css.*HL ;
css2 = css.*HL2 ;

% Observation of aperiodic component

figure(3)
subplot(311)
hold on
plot(css1,'r');
hold off
title('Cepstrum of attenuated cosine (blue - orig, red - lifter)')
subplot(313)
hold on
plot(css1,'r');
hold off
title('Cepstrum of attenuated cosine with a reflection (blue - orig, red - lifter)')

pause(2)

SS1abs = real(exp(fft(css1)));
SS2abs = real(exp(fft(css2)));

figure(2)
subplot(311)
hold on
plot(SS1abs,'r');
hold off
title('Spectrum of attenuated cosine (blue - orig, red - lifter)')

subplot(313)
hold on
plot(SS1abs,'r');
hold off
title('Spectrum of attenuated cosine with a reflection (blue - orig, red - lifter)')

pause(2)

% Observation of periodic component (impulse response)

figure(3)
subplot(312)
hold on
plot(css2,'g');
hold off
title('Cepstrum of impulse response (blue - orig, green - lifter)')

pause(2)

figure(2)
subplot(312)
hold on
plot(SS2abs,'g');
hold off



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

