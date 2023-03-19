% Demo cepstral distance

% Preparing data

% s=[];
% s=[s; audioread('fletna_d_dur_1.wav')];
% s=[s; audioread('fletna_d_dur_2.wav')];
% s=[s; audioread('fletna_d_dur_3.wav')];
% s=[s; audioread('fletna_d_dur_4.wav')];
% s=[s; audioread('fletna_d_dur_5.wav')];
% s=[s; audioread('fletna_d_dur_6.wav')];
% s=[s; audioread('fletna_d_dur_7.wav')];
% s=[s; audioread('fletna_d_dur_8.wav')];
% 
% s=decimate(s,3);
% audiowrite('/myhome/pollak/www/html/vyu/b2m31dsp/signaly/fletna_stupnice.wav',s,fs,'BitsPerSample',16);

% t=(1:fs)./fs;
% s=[sin(2*pi*1000*t) sin(2*pi*2000*t) sin(2*pi*3000*t) sin(2*pi*4000*t) sin(2*pi*5000*t)]';
% s=s+0.01*randn(length(s),1);
% audiowrite('/myhome/pollak/www/html/vyu/b2m31dsp/signaly/sinusovky.wav',s,fs,'BitsPerSample',16);

% Loading prepared data

% [s,fs]=audioread('sinusovky.wav');
[s,fs]=audioread('fletna_stupnice.wav');
% s=loadbin('T10100P0.CS0'); fs= 16000;

% RANDOM
wlen = 512;

figure(31)
subplot(211)
tt=(0:length(s)-1)./fs;
plot(tt,s);
title('Analyzed signal waveform')
xlabel('Time [s]')
xlim([0 max(tt)])
ylim([-max(abs(s)) max(abs(s))])
subplot(212)
spectrogram(s,wlen,[],[],fs,'yaxis');
colorbar off
title('Signal spectrogram')
xlabel('Time [s]')

% Cepstrum computation
% 1)realcepstrum

channels =1;
cp=16;
wstep=wlen;  % For segmentation purposes we need to have successive frames
p=16;

cr=vrceps(s,channels,cp,wlen,wstep);
ca=vaceps(s,channels,p,cp,wlen,wstep);

wnum=length(cr(:,1));
ttt=(0:wnum-1)*wstep./fs;

figure(32)

subplot(211)
% plot(ttt,cr)
plot(ttt,cr(:,2:end))
xlim([0 max(ttt)])
title('Real cepstrum')
xlabel('Time [s]')

subplot(212)
% plot(ttt,ca)
plot(ttt,ca(:,2:end))
xlim([0 max(ttt)])
title('LPC cepstrum')
xlabel('Time [s]')



cdist_seg_1_0=zeros(wnum-1,1);
cdist_seg_2_0=zeros(wnum-1,1);
cdist_seg_1_1=zeros(wnum-1,1);
cdist_seg_2_1=zeros(wnum-1,1);

for ii=1:wnum-1,
    cdist_seg_1_0(ii)=cde(cr(ii,:),cr(ii+1,:));
    cdist_seg_2_0(ii)=cde(ca(ii,:),ca(ii+1,:));

    cdist_seg_1_1(ii)=cde(cr(ii,2:end),cr(ii+1,2:end));
    cdist_seg_2_1(ii)=cde(ca(ii,2:end),ca(ii+1,2:end));
    
end

tttt=(1:wnum-1)*wstep./fs;

figure(33)

subplot(311)
plot(tt,s);
xlim([0 max(tt)])
title('Signal waveform')

subplot(312)
spectrogram(s,wlen,[],[],fs,'yaxis'); colorbar off
title('Signal spectrogram')

subplot(313)
plot(tttt,cdist_seg_1_0,'b')
xlim([0 max(tttt)])
hold on
plot(tttt,cdist_seg_1_1,'r')
xlim([0 max(tttt)])
hold off
title('RCEPS CD of successive frames - with (BLUE) and without (RED) c[0]')


figure(34)

subplot(311)
plot(tt,s);
xlim([0 max(tt)])
title('Signal waveform')

subplot(312)
spectrogram(s,wlen,[],[],fs,'yaxis'); colorbar off
title('Signal spectrogram')

subplot(313)
plot(tttt,cdist_seg_2_0,'b')
xlim([0 max(tttt)])
hold on
plot(tttt,cdist_seg_2_1,'r')
xlim([0 max(tttt)])
hold off
title('LPC CD of successive frames - with (BLUE) and without (RED) c[0]')






