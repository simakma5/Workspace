clear all; close all; clc

f   = 440e6;
GTX = 2;
GRX = 2;
PW  = 0.4;
PTX = 10*log10(PW/1e-3)
d = linspace(1e3, 10000e3, 100);

FSLdB = 20.*log10(d)+20*log10(f)-147.55;

figure;
plot(d./1e3,FSLdB);
xlabel('R [km]');
ylabel('FSL [dB]')

Pr = PTX + GRX +GTX - FSLdB;

figure;
plot(d./1e3,Pr);
xlabel('R [km]');
ylabel('Received Power [dBm]')