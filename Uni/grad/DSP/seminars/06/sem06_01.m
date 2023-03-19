%%
% depdist
clear; close all; clc

s = loadbin('SA001S04.CS0');
fs = 16e3;
wlen = 0.032*fs;  % segment 32 ms

figure(1)
subplot(211);
plot(s)
axis tight

subplot(212);
spectrogram(s,wlen,[],[],fs,'yaxis')
colorbar off; colormap jet

cp = 12;
cepr = vrceps(s,1,cp,wlen,wlen/2);

figure(2)
plot(cepr(:,2:end))  % koeficienty c[0] zkresluji, tak je zahodim
% obecne kdyz nechci, aby intenzita zpracovavaneho signalu zasadne
% ovlivnovala waveformy, tak c[0] zahodim
% nekdy je to vsak potreba zase - zalezi na aplikaci
axis tight

M = 30;
Wc = 0.9;
b = fir1(M,Wc,'low');

figure(3)
freqz(b,1,1e3)

sf = filter(b,1,s);

figure(4)
subplot(211);
spectrogram(sf,wlen,[],[],fs,'yaxis')
colorbar off; colormap jet

ceprf = vrceps(sf,1,cp,wlen,wlen/2);

cdist = cde(cepr(:,2:end),ceprf);
cmean = mean(cdist);

subplot(212);
plot(cdist)
axis tight
title(['Mean CDE: ' num2str(cmean)])

%%
% depdist
clear; close all; clc

s = audioread('sinusovky.wav');
fs = 16e3;
wlen = 0.032*fs;  % segment 32 ms

figure(1)
subplot(211);
plot(s)
axis tight

subplot(212);
spectrogram(s,wlen,[],[],fs,'yaxis')
colorbar off; colormap jet

cp = 12;
cepr = vrceps(s,1,cp,wlen,wlen/2);

figure(2)
plot(cepr(:,2:end))  % koeficienty c[0] zkresluji, tak je zahodim
% obecne kdyz nechci, aby intenzita zpracovavaneho signalu zasadne
% ovlivnovala waveformy, tak c[0] zahodim
% nekdy je to vsak potreba zase - zalezi na aplikaci
axis tight

M = 30;
Wc = 0.9;
b = fir1(M,Wc,'low');

figure(3)
freqz(b,1,1e3)

sf = filter(b,1,s);

figure(4)
subplot(211);
spectrogram(sf,wlen,[],[],fs,'yaxis')
colorbar off; colormap jet

ceprf = vrceps(sf,1,cp,wlen,wlen/2);

cdist = cde(cepr(:,2:end),ceprf);
cmean = mean(cdist);

subplot(212);
plot(cdist)
axis tight
title(['Mean CDE: ' num2str(cmean)])

%% Vimr's take on the last task
s = audioread('sinusovky.wav');
fs = 16e3;
wlen = floor(0.032*fs);  % segment 32 ms

channels = 1;
cp = 16;
% p = 16;
wstep = wlen;
% cr = vaceps(s,channels,p,cp,wlen,wstep);
cr = vrceps(s,channels,cp,wlen,wstep);

wnum = length(cr(:,1));
cdist = zeros(wnum-1,1);

for ii=1:wnum-1
    cdist(ii) = cde(cr(ii,2:end),cr(ii+1,2:end));
end

% tt = wstep*(0:wnum-2)/fs;

figure(31)
subplot(211);
spectrogram(s,wlen,0,wlen,fs,'yaxis');
colormap jet

subplot(212);
plot(cdist)
% xlim([0 max(tt)])  % malinko je to off, ale tohle nefunguje
