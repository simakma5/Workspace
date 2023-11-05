close all; clear; clc; addpath(genpath(fullfile([pwd '\Uni' '\grad' '\RAD' '\sem07'])))

%% Initialization
load DataCv7a.mat
% load DataCv7b.mat
Tp = 25e-6;

%% Processing
% Pulse compression (hamming windowing + matched filtering)
sSumComp = filter(flipud(conj(s0.*hamming(length(s0))))', 1, sSum, [], 1);

% MTI processing (double suppresion / three pulse canceller)
sSumCompMTI = filter([1 -2 1]/4, 1, sSumComp, [], 2);

% MTD processing (two types -- uncomment to compare)
NMTD = 32;  % number of pulses processed
% non-windowed (rectangle) vs. Hamming-windowed MTD filter
% hMTD = repmat(ones(NMTD,1)', NMTD, 1).*exp(1i*2*pi*(0:NMTD-1)'*(0:NMTD-1)/NMTD);
hMTD = repmat(hamming(NMTD)', NMTD, 1).*exp(1i*2*pi*(0:NMTD-1)'*(0:NMTD-1)/NMTD);
sSumCompMTD = cell(NMTD);
for pulse = 1:NMTD
    sSumCompMTD{pulse} = filter(hMTD(pulse,:), 1, sSumComp, [], 2);
end

% Display ampl input/kompress/MTI filterred
figure('name','Signal2D')
subplot(3,1,1)
imagesc(20*log10(abs(sSum)),[0,100])
colormap('jet')
hb=colorbar();
ylabel(hb,'A [dBLSB]')
grid on
title ('Input signal')
xlabel('Azim quanta')
ylabel('Range quanta')
subplot(3,1,2)
imagesc(20*log10(abs(sSumComp)),[0,100])
colormap('jet')
hb=colorbar();
ylabel(hb,'A [dBLSB]')
grid on
title ('Compressed signal')
xlabel('Azim quanta')
ylabel('Range quanta')
subplot(3,1,3)
% imagesc(20*log10(abs(filter([1 -2 1],1,filter(flipud(conj(s0.*hamming(length(s0)))),1,sSum,[],1),[],2))))
imagesc(20*log10(abs(sSumCompMTI)),[0,100])
colormap('jet')
hb=colorbar();
ylabel(hb,'A [dBLSB]')
grid on
title ('MTI filtered signal')
xlabel('Azim quanta')
ylabel('Range quanta')

% Ascope
figure('name','AScopeDkv')
subplot(3,1,1)
plot(20*log10(abs(sSum(:,10))))
title ('A-Scope, on range quantas (in)')
xlabel('R [r. quanta]')
ylabel('A [dBLSB]')
ylim([20,120])
grid on
subplot(3,1,2)
plot(20*log10(abs(sSumComp(:,10))))
title('Compressed signal')
xlabel('R [r. quanta]')
ylabel('A [dBLSB]')
ylim([20,120])
grid on
subplot(3,1,3)
plot(20*log10(abs(sSumCompMTI(:,10))))
xlabel('R [r. quanta]')
ylabel('A [dBLSB]')
title('MTI filtered signal')
ylim([20,120])
grid on

%% Ascope
r=Rmin+(0:size(sSum,1)-1)*1.5e8/fs-Tp*1.5e8;
figure('name','AScope')
subplot(3,1,1)
plot(r/1e3,20*log10(abs(sSum(:,10))))
title ('A-Scope, on range (in)')
xlabel('R [km]]')
ylabel('A [dBLSB]')
ylim([20,120])
grid on
subplot(3,1,2)
plot(r/1e3,20*log10(abs(sSumComp(:,10))))
title('Compressed signal')
xlabel('R [km]]')
ylabel('A [dBLSB]')
ylim([20,120])
grid on
subplot(3,1,3)
plot(r/1e3,20*log10(abs(sSumCompMTI(:,10))))
xlabel('R [km]]')
ylabel('A [dBLSB]')
title('MTI filtered signal')
ylim([20,120])
grid on

figure('name','MTD')
for ii=1:NMTD
  subplot(2,NMTD/2,ii)
  sx=20*log10(abs(sSumCompMTD{ii}))-60;
  sx(sx<0)=0;
  imagesc(sx(:,16:end),[0,90]);
  xlabel([num2str((ii-1)/NMTD*PRF)])
%   xlabel('Azim Q.');
%   ylabel('R. Q.');
  colormap('jet')
end
%%
figure(10)
plot(real(fft([hMTD zeros(NMTD,1024-NMTD)]')))
figure(11)
plot(real(fft([hMTD(10,:) zeros(1,1024-NMTD)]')))
hold on
plot(imag(fft([hMTD(10,:) zeros(1,1024-NMTD)]')))
hold off

figure(12)
plot(abs(fft([1 -1 zeros(1,1024-2)])),'k')
hold on
plot(abs(fft([1 -2 1 zeros(1,1024-3)])),'b')
plot(abs(fft([1 -3 3 -1 zeros(1,1024-4)])),'r')
hold off
