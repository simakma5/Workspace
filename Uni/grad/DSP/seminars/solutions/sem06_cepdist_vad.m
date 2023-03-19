% Demo cepstral distance

s=loadbin('SA001S04.CS0');
% s=loadbin('SA001S04.CS1');
% s=loadbin('SA001S04.CS2');
% s=loadbin('SA001S04.CS3');

fs=16000;

tt=(0:length(s)-1)./fs;


wlen=512;

channels =1;
cp=16;
wstep=wlen/2;
p=16;

cr=vrceps(s,channels,cp,wlen,wstep);
ca=vaceps(s,channels,p,cp,wlen,wstep);

cravg=mean(cr(1:20,:));
caavg=mean(ca(1:20,:));

wnum=length(cr(:,1));

    
cd_real0 = cde ( cr, cravg ) ;
cd_real1 = cde ( cr(:,2:end), cravg(2:end) ) ;
    
cd_lpc0 = cde ( ca, caavg ) ;
cd_lpc1 = cde ( ca(:,2:end), caavg(2:end) ) ;
    


figure(21)

subplot(211)
plot(tt,s);
xlim([0 max(tt)])
title('Signal waveform')
subplot(212)
spectrogram(s,wlen,[],[],fs,'yaxis'); colorbar off;
title('Signal spectrogram')

figure(22)

subplot(311)
plot(cr); 
xlim([0 wnum])
title('REAL cepstrum of signal')
xlabel('i - frame No. [-]')
subplot(312)
plot(cd_real0); 
xlim([0 wnum])
title('CDE RCEPS from background')
xlabel('i - frame No. [-]')
subplot(313)
plot(cd_real1); 
xlim([0 wnum])
title('CDE RCEPS from background (without c[0])')
xlabel('i - frame No. [-]')

figure(23)

subplot(311)
plot(ca); 
xlim([0 wnum])
title('LPC cepstrum of signal')
xlabel('i - frame No. [-]')
subplot(312)
plot(cd_lpc0); 
xlim([0 wnum])
title('LPC cepstrum of distorted signal')
xlabel('i - frame No. [-]')
subplot(313)
plot(cd_lpc1); 
xlim([0 wnum])
title('LPC cepstrum of distorted signal (without c[0])')
xlabel('i - frame No. [-]')

%%

figure(24)

subplot(211)
plot(cd_real0)
hold on
plot(speechpart,cd_real0(speechpart),'r')
hold off
xlim([0 wnum])
title(['RCEPS - Euclidan distance: ' num2str(cd_real0_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

subplot(212)
plot(cd_real1)
hold on
plot(speechpart,cd_real1(speechpart),'r')
hold off
xlim([0 wnum])
title(['RCEPS - Eucl. dist. without c[0]: ' num2str(cd_real1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

figure(25)

subplot(211)
plot(cd_lpc0)
hold on
plot(speechpart,cd_lpc0(speechpart),'r')
hold off
xlim([0 wnum])
title(['LPCC - Euclidan distance: ' num2str(cd_lpc0_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

subplot(212)
plot(cd_lpc1)
hold on
plot(speechpart,cd_lpc1(speechpart),'r')
hold off
xlim([0 wnum])
title(['LPCC - Eucl. dist. without c[0]: ' num2str(cd_lpc1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

