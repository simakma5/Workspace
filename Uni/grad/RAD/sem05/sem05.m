close all; clear; clc; addpath(genpath(fullfile([pwd '\Uni' '\grad' '\RAD' '\sem05'])))

c=3e8;
data = importdata('Cv3_SrecBN.mat');
sSamp = data.s;                     % samples of the Tx signal model's complex envelope [-]
sEnv = data.s0;                     % complex envelope of the Tx signal's replica [-]
fSamp = data.fs;                    % sampling frequency [Hz]
PRF = data.fop;                     % pulse repetition frequency [Hz]
PW = data.Tp;                       % pulse width [s]
time = (1:length(sEnv))/fSamp;

%%
figure(1)
tiles = tiledlayout(2,1);

nexttile
stem(time, abs(sEnv))
xlabel('Time (s)')
ylabel('Amplitude (-)')

nexttile
stem(time, unwrap(angle(sEnv)))
xlabel('Time (s)')
ylabel('\psi (rad)')

title(tiles, 'Received signal')
%%
sEnvHamm = sEnv.*hamming(length(sEnv));
sEnvAutocorr_dB = 20*log10(abs(xcorr(sEnv))/max(abs(xcorr(sEnv))));
sEnvXcorr1_dB = 20*log10(abs(xcorr(sEnv,sEnvHamm))/max(abs(xcorr(sEnv,sEnvHamm))));
sEnvXcorr2_dB = 20*log10(abs(xcorr(sEnv,sEnvHamm))/max(abs(xcorr(sEnv))));

figure(2)
plot(sEnvAutocorr_dB);
hold on
plot(sEnvXcorr1_dB, '-r');
plot(sEnvXcorr2_dB, '-g');
hold off
xlabel('Samples (-)')
ylabel('Amplitude (dB)')
ylim([-100 0])
title('Normalized correlation functions')

%%
sSampXcorr = abs(xcorr(sSamp,sEnv));

figure(3);
distance = (-length(sSamp)+1:length(sSamp)-1)*c/2/fSamp;
plot(distance*1e-3, sSampXcorr)
xlabel('[km]')
title('Cross-correlation of signal samples with its envelope')

%%
sFiltered = filter(conj(flipud(sEnv)), 1, sSamp);

figure(4)
tiles = tiledlayout(2, 1);

nexttile
distance = (0:length(sSamp)-1)*c/2/fSamp;
plot(distance*1e-3, abs(sSamp));
xlabel('[km]')

nexttile
distance = (-length(sEnv)+1:length(sFiltered)-length(sEnv))*c/2/fSamp;
plot(distance*1e-3, abs(sFiltered));
xlabel('[km]')

title(tiles, 'Range detection')

%%
Df=1e6;
% Df=1/Tp;
DR=c/2/Df; %rozlisovaci schopnost [m]/ range resolution

DR_v=round((fSamp/Df+2)); %vzorky/samples
stest=[sEnv; zeros(4*length(sEnv),1)]+[zeros(DR_v,1); sEnv; zeros(4*length(sEnv)-DR_v,1)];

figure(5)
sftest=filter(conj(flipud(sEnv)),1,stest);
R2=(-length(sEnv)+1:length(sftest)-length(sEnv))'/fSamp*c/2;
plot(R2,abs(sftest));
