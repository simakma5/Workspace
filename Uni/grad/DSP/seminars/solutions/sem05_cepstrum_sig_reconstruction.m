% Demo cepstral analysis
%    REFLECTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Parameters of attenuated cosine
f=500 ;
fs=8000 ;
Tau=0.001 ;

% Generation of time and cosine (result should be COLUMN VECTORS !!!)
N=100 ;
t=(0:N-1)'./fs ;
s=cos(2*pi*f*t).*exp(-t/Tau) ;

% Parameters and impulse response of filter for refelction modeling
M=15; 
m=0.7 ;

h=zeros(N,1); 
h(1)=1;
h(M+1)=m;

% Creation of reflection

ss=filter(h,1,s);


figure(21)
subplot(211)
plot(s);
title('Signal: blue - original (attenuated cosine)')
subplot(212)
plot(ss)
title('Signal: blue - distorted (att. cosine with reflection)')

pause(2)

Sabs =abs(fft(s));
SSabs =abs(fft(ss));

figure(22)
subplot(211)
plot(Sabs);
title('Spectrum: blue - original')
subplot(212)
plot(SSabs);
title('Spectrum: blue - distorted')

pause(2)

cs=rceps(s) ;
css=rceps(ss) ;

figure(23)
subplot(211)
plot(cs);
title('RCEPS: blue - original')
subplot(212)
plot(css);
title('RCEPS: blue - distorted')


pause(2)

% Lifter
%%%%%%%%

nc = 12 ;
HL = [ones(nc+1,1); zeros(N-2*nc-1,1); ones(nc,1) ];

% Liftrace
%%%%%%%%%%

css1 = css.*HL ;

figure(23)
subplot(211)
hold on
plot(css1,'r');
hold off
title('RCEPS: blue - original, red - liftered distorted')
subplot(212)
hold on
plot(css1,'r');
hold off
title('RCEPS: blue - distorted, red - liftered distorted')

pause(2)

SS1abs = real(exp(fft(css1)));

figure(22)
subplot(212)
hold on
plot(SS1abs,'r');
hold off
title('Spectrum: blue - original, red - liftered distorted')
subplot(211)
hold on
plot(SS1abs,'r');
hold off
title('Spectrum: blue - distorted, red - liftered distorted')

pause(2)

% Signal reconstruction using phase of distorted signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% phase of distorted signal (bad estimation)
SSangle = angle(fft(ss)) ;

% phase of clean signal (ideal estimetion - but not available)
SSangleid = angle(fft(s)) ;


s_out = ifft ( SS1abs.*exp(j*SSangle) ) ;
s_out2 = ifft ( SS1abs.*exp(j*SSangleid) ) ;

figure(21)
subplot(211)
hold on
plot(s_out,'r');
title('Signal: blue - original, red - reconstuct. distorted phase')
hold off
subplot(212)
hold on
plot(s_out,'r');
title('Signal: blue - distorted, red - reconstuct. distorted phase')
hold off


pause(2)

figure(21)
subplot(211)
hold on
plot(s_out2,'g');
title('Signal: blue - original, red - reconstuct. distorted phase, green - ideal phase')
hold off
subplot(212)
hold on
plot(s_out2,'g');
title('Signal: blue - distorted, red - reconstuct. distorted phase, green - ideal phase')
hold off


pause(2)

%% Signal reconstruction using complex cepstrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% complex spectrum is not symmetric (contains phase information)
cc = cceps (ss) ;

% Lifter is not symmetric for complex cepstrum
HLc = [ones(nc+1,1); zeros(N-nc-1,1) ];
cc1 = cc.*HLc ;

% reconstructed signal from complex cepstrum
s_out = ifft ( exp(fft(cc1)) ) ;

% Results observation

figure(24)
subplot(211)
plot(s);
title('Signal: blue - original')
subplot(212)
plot(ss);
title('Signal: blue - distorted')


pause(2)

figure(25);
plot(cc); 
title('CCEPS: blue - distorted')

pause(2)

figure(25);

hold on
plot(cc1,'r');
title('CCEPS: blue - distorted, red - liftered')
hold off

pause(2)

figure(24)

subplot(211)
hold on
plot(s_out,'r');
title('Signal: blue - original, red - reconstructed from cceps')
hold off
subplot(212)
hold on
plot(s_out,'r');
hold off
title('Signal: blue - distorted, red - reconstructed from cceps')






