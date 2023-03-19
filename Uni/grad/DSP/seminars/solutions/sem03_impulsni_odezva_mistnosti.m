% demo - impulsni odersva mistnosti

x0=loadbin('SA001S04.CS0'); x0=x0(2e4:7e4-1);
x1=loadbin('SA001S04.CS1'); x1=x1(2e4:7e4-1);
x2=loadbin('SA001S04.CS2'); x2=x2(2e4:7e4-1);
x3=loadbin('SA001S04.CS3'); x3=x3(2e4:7e4-1);
fs=16000;

figure(20)
subplot(311)
plot([x0 x1])
title('Signals CS0-CS1')

subplot(312)
plot([x0 x2])
title('Signals CS0-CS2')

subplot(313)
plot([x0 x3])
title('Signals CS0-CS3')

figure(21)
subplot(221)
spectrogram(x0,512,[],[],fs,'yaxis'); colorbar off; colormap jet
title('Spectrogram CS0')
subplot(222)
spectrogram(x1,512,[],[],fs,'yaxis'); colorbar off; colormap jet
title('Spectrogram CS1')
subplot(223)
spectrogram(x2,512,[],[],fs,'yaxis'); colorbar off; colormap jet
title('Spectrogram CS2')
subplot(224)
spectrogram(x3,512,[],[],fs,'yaxis'); colorbar off; colormap jet
title('Spectrogram CS3')

wlen=512;

[Hx01,FF] = tfestimate(x0,x1,wlen,[],[],fs,'twosided') ;
impx01 = (ifft(Hx01)) ;
Hx02 = tfestimate(x0,x2,wlen,[],[],fs,'twosided') ;
impx02 = (ifft(Hx02)) ;
Hx03 = tfestimate(x0,x3,wlen,[],[],fs,'twosided') ;
impx03 = (ifft(Hx03)) ;


figure(22)
subplot(211)
plot(FF,20*log10(abs(Hx01)))
title('Mag. freq. response in dBs: CS0-CS1')
subplot(212)
plot(impx01)
title('Ipulse response: CS0-CS1')

figure(23)
subplot(211)
plot(FF,20*log10(abs(Hx02)))
title('Mag. freq. response in dBs: CS0-CS2')
subplot(212)
plot(impx02)
title('Ipulse response: CS0-CS2')

figure(24)
subplot(211)
plot(FF,20*log10(abs(Hx03)))
title('Mag. freq. response in dBs: CS0-CS3')
subplot(212)
plot(impx03)
title('Ipulse response: CS0-CS3')

%% Modelling of distortion

figure(30)

plot(x0,'k:')
hold on
plot(x1,'b')
plot(filter(impx01,1,x0),'r')
hold off
legend('source x0','original x1','modelled x1')

figure(31)

plot(x0,'k:')
hold on
plot(x2,'b')
plot(filter(impx02,1,x0),'r')
hold off
legend('source x0','original x2','modelled x2')

figure(32)

plot(x0,'k:')
hold on
plot(x3,'b')
plot(filter(impx03,1,x0),'r')
hold off
legend('source x0','original x3','modelled x3')



