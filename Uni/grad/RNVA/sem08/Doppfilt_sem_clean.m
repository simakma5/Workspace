clear variables
clear global
clc
close all

Tp = 25e-6;

load DataCv7a01



% Display ampl input/kompress/MTI filterred
figure('name','Signal2D')
imagesc(20*log10(abs(sSum)),[0,100])
colormap('jet')
hb=colorbar();
ylabel(hb,'A [dBLSB]')
grid on
title ('Input signal')
xlabel('Azim quanta')
ylabel('Range quanta')
