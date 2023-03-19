% Demo cepstral distance

s0=loadbin('SA001S04.CS0');
s1=loadbin('SA001S04.CS1');
s2=loadbin('SA001S04.CS2');
s3=loadbin('SA001S04.CS3');

% b=fir1(50,0.6,'low');
b=fir1(50,0.1,'high');

s0f=filter(b,1,s0);

fs=16000;
wlen=512;

channels =1;
cp=16;
wstep=wlen/2;
p=16;

cr0=vrceps(s0,channels,cp,wlen,wstep);
ca0=vaceps(s0,channels,p,cp,wlen,wstep);

cr0f=vrceps(s0f,channels,cp,wlen,wstep);
ca0f=vaceps(s0f,channels,p,cp,wlen,wstep);

wnum=length(cr0(:,1));


figure(1)

subplot(211)
spectrogram(s0,wlen,[],[],fs,'yaxis'); colorbar off;
title('Original signal (reference, channel CS0 ')
subplot(212)
spectrogram(s0f,wlen,[],[],fs,'yaxis'); colorbar off;
title('Distorted signal')

figure(2)

subplot(211)
plot(cr0); 
xlim([0 wnum])
title('REAL cepstrum of original signal')
xlabel('i - frame No. [-]')
subplot(212)
plot(cr0f); 
xlim([0 wnum])
title('REAL cepstrum of distorted signal')
xlabel('i - frame No. [-]')

figure(3)

subplot(211)
plot(ca0); 
xlim([0 wnum])
title('LPC cepstrum of original signal')
xlabel('i - frame No. [-]')
subplot(212)
plot(ca0f); 
xlim([0 wnum])
title('LPC cepstrum of distorted signal')
xlabel('i - frame No. [-]')

figure(4)
i=105;
subplot(211)
plot(0:cp,cr0(i,:),'.-b');
title(['Real cepstrum of frame No. ' num2str(i)])
subplot(212)
plot(0:cp,ca0(i,:),'.-r');
title(['LPC cepstrum of frame No. ' num2str(i)])

 
cd_real0 = cde ( cr0, cr0f ) ;
cd_real1 = cde ( cr0(:,2:end), cr0f(:,2:end) ) ;
    
cd_lpc0 = cde ( ca0, ca0f ) ;
cd_lpc1 = cde ( ca0(:,2:end), ca0f(:,2:end) ) ;
    


% Speech part is approx in frames No. 90:250
speechpart=85:260;
cd_real0_avg = mean ( cd_real0(speechpart) ) ;
cd_real1_avg = mean ( cd_real1(speechpart) ) ;
cd_lpc0_avg = mean ( cd_lpc0(speechpart) ) ;
cd_lpc1_avg = mean ( cd_lpc1(speechpart) ) ;

figure(5)

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

figure(6)

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



%% Distortion of real signals from various microphones

cr1=vrceps(s1,channels,cp,wlen,wstep);
ca1=vaceps(s1,channels,p,cp,wlen,wstep);
cr2=vrceps(s2,channels,cp,wlen,wstep);
ca2=vaceps(s2,channels,p,cp,wlen,wstep);
cr3=vrceps(s3,channels,cp,wlen,wstep);
ca3=vaceps(s3,channels,p,cp,wlen,wstep);

  
cd_cs1_real1 = cde ( cr0(:,2:end), cr1(:,2:end) ) ;
cd_cs1_lpc1 = cde ( ca0(:,2:end), ca1(:,2:end) ) ;
    
cd_cs2_real1 = cde ( cr0(:,2:end), cr2(:,2:end) ) ;
cd_cs2_lpc1 = cde ( ca0(:,2:end), ca2(:,2:end) ) ;
    
cd_cs3_real1 = cde ( cr0(:,2:end), cr3(:,2:end) ) ;
cd_cs3_lpc1 = cde ( ca0(:,2:end), ca3(:,2:end) ) ;
    
% Speech part is approx in frames No. 90:250
speechpart=85:260;

cd_cs1_real1_avg = mean ( cd_cs1_real1(speechpart) ) ;
cd_cs1_lpc1_avg = mean ( cd_cs1_lpc1(speechpart) ) ;
cd_cs2_real1_avg = mean ( cd_cs2_real1(speechpart) ) ;
cd_cs2_lpc1_avg = mean ( cd_cs2_lpc1(speechpart) ) ;
cd_cs3_real1_avg = mean ( cd_cs3_real1(speechpart) ) ;
cd_cs3_lpc1_avg = mean ( cd_cs3_lpc1(speechpart) ) ;

figure(10);

subplot(211)
spectrogram(s0,wlen,[],[],fs,'yaxis'); colorbar off;
c=caxis;
title('Reference channel - CS0 (head-set)')
subplot(212)
spectrogram(s1,wlen,[],[],fs,'yaxis'); colorbar off;
caxis(c);
title('Nokia hands-free - CS1 (close-talk, omnidirectional)')

figure(11);

subplot(211)
spectrogram(s0,wlen,[],[],fs,'yaxis'); colorbar off;
c=caxis;
title('Reference channel - CS0 (head-set)')
subplot(212)
spectrogram(s2,wlen,[],[],fs,'yaxis'); colorbar off;
caxis(c);
title('Sehnheiser - CS2 (middle-talk, directional)')

figure(12);

subplot(211)
spectrogram(s0,wlen,[],[],fs,'yaxis'); colorbar off;
c=caxis;
title('Reference channel - CS0 (head-set)')
subplot(212)
spectrogram(s3,wlen,[],[],fs,'yaxis'); colorbar off;
caxis(c);
title('Sehnheiser - CS3 (far-talk, directional)')

figure(13)

subplot(211)
plot(cd_cs1_real1)
hold on
plot(speechpart,cd_cs1_real1(speechpart),'r')
hold off
xlim([0 wnum])
title(['CS1 - RCEPS CDE without c[0]: ' num2str(cd_cs1_real1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

subplot(212)
plot(cd_cs1_lpc1)
hold on
plot(speechpart,cd_cs1_lpc1(speechpart),'r')
hold off
xlim([0 wnum])
title(['CS1 - LPCC CDE without c[0]: ' num2str(cd_cs1_lpc1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

figure(14)

subplot(211)
plot(cd_cs2_real1)
hold on
plot(speechpart,cd_cs2_real1(speechpart),'r')
hold off
xlim([0 wnum])
title(['CS2 - RCEPS CDE without c[0]: ' num2str(cd_cs2_real1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

subplot(212)
plot(cd_cs2_lpc1)
hold on
plot(speechpart,cd_cs2_lpc1(speechpart),'r')
hold off
xlim([0 wnum])
title(['CS2 - LPCC CDE without c[0]: ' num2str(cd_cs2_lpc1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

figure(15)

subplot(211)
plot(cd_cs3_real1)
hold on
plot(speechpart,cd_cs3_real1(speechpart),'r')
hold off
xlim([0 wnum])
title(['CS3 - RCEPS CDE without c[0]: ' num2str(cd_cs3_real1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')

subplot(212)
plot(cd_cs3_lpc1)
hold on
plot(speechpart,cd_cs3_lpc1(speechpart),'r')
hold off
xlim([0 wnum])
title(['CS3 - LPCC CDE without c[0]: ' num2str(cd_cs3_lpc1_avg) ' - average for speech part'])
xlabel('i - frame No. [-]')



