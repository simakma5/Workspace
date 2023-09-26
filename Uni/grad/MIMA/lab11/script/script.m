close all; clear; clc
c = 3e8;

%% Primkovy rezonator
L = 34e-3;
f = [3.1e9 6.2e9 9.3e9 12.4e9 15.5e9 18.5e9 21.6e9];
lambda = [2*L/1 2*L/2 2*L/3 L/2 2*L/5 L/3 2*L/7];
epsilon_eff = (c./(f.*lambda)).^2;
figure(1)
plot(epsilon_eff,'-o')
axis tight
grid on
grid minor
ylim([1.96 2.1])
xlabel('Frequency [GHz]')
ylabel('\epsilon_{eff} [-]')

%% Prstencovy rezonator
R = 34e-3/2;
w = 0.7e-3;
C = 2*pi*(R-w);
f = [1.72e9 3.46e9 5.2e9 6.93e9 8.63e9 10.33e9 12.04e9 13.72e9 15.41e9 17.09e9 18.75e9 20.41e9 22.06e9 23.71e9];
lambda = [C C/2 C/3 C/4 C/5 C/6 C/7 C/8 C/9 C/10 C/11 C/12 C/13 C/14];
epsilon_eff = (c./(f.*lambda)).^2;
figure(2)
plot(f*1e-9,epsilon_eff,'-o')
axis tight
grid on
grid minor
ylim([2.8 3.05])
xlabel('Frequency [GHz]')
ylabel('\epsilon_{eff} [-]')
