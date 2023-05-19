close all; clear; clc

%%
% nacteni/load
load RAD_Sem05_SrecB_0p2.mat;
load DataDoppfilt.mat

t = (1:length(s0))'/fs;

%%
figure(1)
plot(abs(sSum(:,1)))

%%
figure(2)
subplot(211)
plot(abs(s0))
subplot(212)
plot(unwrap(angle(s0)))

figure(3)
subplot(211)
plot(real(s0))
subplot(212)
plot(imag(s0))

%%
Rs0 = xcorr(s0,s0);
figure(4)
plot(abs(Rs0))


%%
CPs = xcorr(s,s0);
figure(5)
plot(abs(CPs))

%%
stest_komp = filter(flipud(conj(s0))',1,sSum,[],1); % komprese
DFsCP = filter([1 -2 1],1,stest_komp,[],2);         % Doppler filtrace

figure(6)
plot(abs(DFsCP))
